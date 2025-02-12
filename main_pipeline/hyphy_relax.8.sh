##Using hyphy RELAX and BUSTED to determine strength of selection and purifying vs. positive selection, respectively
##Create a central directory for these analyses:
mkdir ~/data/main_pipeline/hyphy_out
mkdir ~/data/main_pipeline/macse_out

## Before running macse, double check that the nucleotide_HOG .ffn files are all unwrapped.
> tmp
ls ~/data/main_pipeline/orthofinder_out/for_downstream/nucleotide_HOGs/*.ffn | \
while read file; do
echo $file >> tmp
a=$(grep -c "^>" $file)
b=$(wc -l < $file)
if (( $(($a*2))!=$b )); then
echo "no bueno!" >> tmp;
fi
done
grep -B 1 "no bueno!" tmp
## All good!
rm tmp


#####################
### Generate MSAs ###
#####################

##Run macse on the single-copy genes.
conda deactivate
conda activate macse
outdir=~/data/main_pipeline/macse_out
mkdir ${outdir}/trimmed_nucleotideHOGs
orthofinder_outdir=~/data/main_pipeline/orthofinder_out/for_downstream
indir=${orthofinder_outdir}/nucleotide_HOGs
cd $outdir
cat ${orthofinder_outdir}/hogs_passed_3_thresh.list | \
while read hog; do
## Subset the .ffn to only those genomes that had exactly 1 locus in this hog
cat ${orthofinder_outdir}/Taz.all_hogs_which_genomes_single_copy.tsv | \
awk -v hog=$hog -F "\t" ' $1==hog { for (i=2;i<=NF;++i) {print $i} } ' > locus_list.tmp
cat ${orthofinder_outdir}/free_living.all_hogs_which_genomes_single_copy.tsv | \
awk -v hog=$hog -F "\t" ' $1==hog { for (i=2;i<=NF;++i) {print $i} } ' >> locus_list.tmp
grep --no-group-separator -A1 -Ff locus_list.tmp ${indir}/${hog}.ffn > ${outdir}/trimmed_nucleotideHOGs/${hog}.ffn
rm locus_list.tmp
macse -prog alignSequences -seq ${outdir}/trimmed_nucleotideHOGs/${hog}.ffn -out_AA ${hog}.macse.faa -out_NT ${hog}.macse.fna
done

## Remove the final stop codon from each of the sequences.
## BUSTED and RELAX don't allow the final stop codon, though nonREV appears to be fine with it.
## But will remove it up front for consistency.
## Replace with 3 dashes to keep alignment correct.
## Also change the header to only the locus tag because hyphy was getting confused by the multiple field header lines
## Also replace the locus_tags with just the species names because I am using the species tree, not gene trees.
mkdir ~/data/main_pipeline/hyphy_out/msas
cd ~/data/main_pipeline/hyphy_out/msas

ls ~/data/main_pipeline/macse_out/*.fna | \
while read msa; do
og=${msa%.macse.fna}
og=${og##*/}
cat $msa | \
awk ' {if ($0 ~ "^>") {printf "\n"$0"\n"} else { printf $0} }; END {printf "\n"} ' | \
tail -n +2 | \
awk '$0 ~ "^>" {print $0"~placeholder~"; next}; {print $0}' | \
tr -d "\n" | \
awk '{print $0}' | \
sed 's/>/\n>/g' | \
tail -n +2 | \
sed -E 's/TAG(-*)$|TGA(-*)$|TAA(-*)$/---\1\2\3/' | \
sed 's/~placeholder~/\n/' | \
awk -F "[[:space:]]" '{print $1}' | \
awk -F "_" '{print $1}' \
> ${og}.noStop.fa
done


##########################
#### Prune gene trees ####
##########################

## I am using the species tree instead of gene trees (orthofinder gives both)
## because I trust the species tree more and we think the potential drawbacks of using the species tree
## (incomplete lineage sorting, HGT, to name 2) are unlikely in this system.
## Not all of these hogs contain a locus from every species
## So must prune the tree to match the hog. 

msa_dir=~/data/main_pipeline/hyphy_out/msas
in_tree=~/data/main_pipeline/orthofinder_out/final_run/Species_Tree/SpeciesTree_rooted_node_labels.txt
tree_outdir=~/data/main_pipeline/hyphy_out/trees
logfile=gotree_prune.log
mkdir $tree_outdir
cd $tree_outdir
conda deactivate
conda activate gotree
## Loop through each N9.HOG (so each line in loci2keep)
## Start at line 2 to skip header
echo "cat all error messages from gotree prune" > $logfile
ls ${msa_dir}/*.noStop.fa | \
while read msa_file; do
hog=${msa_file##*/}
hog=${hog%.noStop.fa}
echo $hog >> $logfile

#Get the list of species to be kept by gotree prune.
grep "^>" $msa_file | sed 's/^>//' | sed 's/[[:space:]].*$//' > species2keep.tmp 

##Prune the corresponding gene tree using gotree
#I am running the trees through gotree prune even if they don't have any leaves to prune to keep everything constant
#Indeed, the input and output files are not identical, although visualizing the trees they look the same to me.
gotree prune -r -i $in_tree -f species2keep.tmp -o ./${hog}.pruned_tree.txt >> $logfile

done
rm *.tmp
#error log contains only the header and the name of each hog. Nothing else printed! success!
rm $logfile

###################
### Label trees ###
###################


## Label the test and reference groups.
## All nodes that contain exclusively T. az. are test except the root of the T. az. zlade which is unlabeled
## All other nodes are labeled reference
conda deactivate
indir=~/data/main_pipeline/hyphy_out/trees
cd $indir
## label_trees_not_ingroup_root.sh checks that the ingroup is monophyletic
## This is not necessary given the species tree is always the starting tree here
## It was put in for gene trees, but should still work fine.
logfile=label-tree.log
echo "label single copy gene trees" > $logfile
ls *.pruned_tree.txt | \
while read in_tree; do
hog=${in_tree%.pruned_tree.txt}
out_tree=${hog}.labeled_tree.txt
echo "$hog" >> $logfile
bash -i ~/scripts/helper/label_trees_not_ingroup_root.sh -i $in_tree -o $out_tree -l $logfile\
 -I ~/Taz.short_prefix.list -O ~/subtree_freeLiving.short_prefix.list
done

## Check one set of erros from helper script
ls *pruned_tree.txt | wc -l
grep --no-group-separator -A1 "Are descendent nodes always before ancestor nodes?" $logfile | grep -v "^Are" | sort | uniq -c
## All yes
## Get the no HOGs
grep -B2 --no-group-separator "^No!" $logfile | awk 'NR %3==1'
## Empty
## logfile looks good
rm $logfile *.pruned_tree.txt


#################
### Run RELAX ###
#################

##Run with test and background branches labeled.

#set variables and make folders
conda deactivate
msa_dir=~/data/main_pipeline/hyphy_out/msas
tree_dir=~/data/main_pipeline/hyphy_out/trees
outdir=~/data/main_pipeline/hyphy_out/relax_out
mkdir $outdir
cd $outdir

ls ${tree_dir}/*.labeled_tree.txt | \
#grep -vFf already_run.tmp | \
## subset to run in parallel
#awk 'NR%2==1' | \
while read tree_file; do
hog=${tree_file##*/}
hog=${hog%.labeled_tree.txt}

msa_file=${msa_dir}/${hog}.noStop.fa
## Only run on MSAs that do not contain !s.
## Will deal with those separately.
if ! grep -q "!" $msa_file; then
outfile=${hog}.relax_out.json
stdoutfile=${hog}.relax_out.txt

## Run RELAX
## --models Minimal just hypothesis tests and doesn't try to fit other models (not sure what this means....see the tutorial)
## --starting-points 10 and srv Yes were recommended by @spond to me in github issue #1661
hyphy relax --alignment $msa_file --tree $tree_file --test "test" --reference "reference"\
 --models Minimal --srv Yes --output $outfile > $stdoutfile

## finish if statement that no exclamation mark
fi
##Finish loop through msa files
done

#############################
### Re-run RELAX for HOGs ### 
### that did not converge ###
#############################

##Identify RELAX results that did not converge
cd ~/data/main_pipeline/hyphy_out/relax_out
> unconverged_HOGs.list
> converged_HOGs.list
for file in N9.HOG*.txt; do
#Take the last line with text on it
cat $file | sed '/^$/d' | tail -n 1 > tmp
#Check if the last line has a message about convergence
if grep -q "converge" tmp; then echo $file >> unconverged_HOGs.list; fi
##Check that the other files all end as I expect
if grep -q "^----" tmp; then echo $file >> converged_HOGs.list; fi
done

echo "Check for error messages:"
grep -i "error" *
## Some runs with errors can just be re-run with the normal call and will be fine.
## Having run ~2000 HOGs several times on the same HOGs, it seems that sometimes the omega
## estimates get too extreme and cause an error. But then this is fixed on re-run.
##Check that I got all files with those two greps
wc -l unconverged_HOGs.list
# 13
wc -l converged_HOGs.list 
# 2009
ls *relax_out.txt | wc -l
# 2022 = 13 + 2009
#Success.
rm converged_HOGs.list

##Also check that all the messages about convergence are identical.
> tmp
cat unconverged_HOGs.list | \
while read file; do
cat $file | sed '/^$/d' | tail -n 1 >> tmp
done
sort -u tmp
##All are identical.
rm tmp

#Rerun the unconverged files using a relax call recommended by spond in a github issue (https://github.com/veg/hyphy/issues/1581)
#Save the original outputs
mkdir unconverged_relax_out
cat unconverged_HOGs.list | \
while read file; do
rm ${file%.txt}.json
mv $file unconverged_relax_out
done
#set variables and make folders
conda deactivate
msa_dir=~/data/main_pipeline/hyphy_out/msas
tree_dir=~/data/main_pipeline/hyphy_out/trees
#Loop and run!
cat unconverged_HOGs.list | \
while read file; do
hog=${file%.relax_out.txt}
outfile=${hog}.relax_out.json
stdoutfile=${hog}.relax_out.txt
## Run RELAX
## --models Minimal just hypothesis tests and doesn't do a couple other things thus saving time
## Got the command options for hyphy from https://github.com/veg/hyphy/issues/1581
hyphy CPU=1 relax --starting-points 100 --grid-size 2000 --models Minimal --srv Yes --alignment ${msa_dir}/${hog}.noStop.fa \
--tree ${tree_dir}/${hog}.labeled_tree.txt --test "test" --reference "reference" --output $outfile > $stdoutfile
done

##Identify RELAX results that did not converge after second attempt
cd ~/data/main_pipeline/hyphy_out/relax_out
> unconverged_HOGs.list
> successfully_run_HOGs.list
for file in N9.HOG*.txt; do
#Take the last line with text on it
cat $file | sed '/^$/d' | tail -n 2 > tmp
#Check if the last line has a message about convergence
if grep -q "converge" tmp; then echo $file >> unconverged_HOGs.list; fi
##Check that the other files all end as I expect
if grep -q "^----" tmp && grep -qi "evidence" tmp; then echo $file >> successfully_run_HOGs.list; fi
done
##Check for files that had other errors or problems than convergence
cat unconverged_HOGs.list successfully_run_HOGs.list > tmp
ls *.relax_out.txt | grep -vFf tmp > other_problems_HOGs.list
rm tmp
 
##Check that I got all files with those two greps
wc -l unconverged_HOGs.list
# 8 
wc -l successfully_run_HOGs.list 
# 2014
wc -l other_problems_HOGs.list
# 0
ls *relax_out.txt | wc -l
# 2022 = 2014 + 8 + 0
#Success.

##Also check that all the messages about convergence are identical.
> tmp
cat unconverged_HOGs.list | \
while read file; do
cat $file | sed '/^$/d' | tail -n 1 >> tmp
done
sort -u tmp
##All are identical.
rm tmp
echo "Check for error messages"
grep -i "error" *
##Move the unconverged outputs
mv unconverged_relax_out unconverged_relax_out_first_attempt
mkdir unconverged_relax_out_second_attempt
cat unconverged_HOGs.list | \
while read file; do
rm ${file%.txt}.json
mv $file unconverged_relax_out_second_attempt
done


######################
### Deal with MSAs ###
## That contain !s ###
######################

bash -i ~/scripts/helper/hyphy_relax_exclamation_marks.sh\
 -m ~/data/main_pipeline -I ~/Taz.short_prefix.list -O ~/subtree_freeLiving.short_prefix.list


##########################
### Parse and Annotate ###
##### Relax Results ######
##########################

## Parse results and put into one table
## Each row is HOG,k-val,p-val,result
echo -e "HOG\tk-val\tp-val" > results.tsv
cat successfully_run_HOGs.list | \
sed 's/\.relax_out.txt//' | \
while read hog; do
json=${hog}.relax_out.json
## Make sure that there is only one line in the json that matches each key I am using
if (( $(grep -c "p-value" $json) == 1  )) && (( $(grep -c "relaxation or intensification parameter" $json) == 1)); then
p=$(grep "p-value" $json | sed 's/.*://' | sed 's/,.*//')
k=$(grep "relaxation or intensification parameter" $json  | sed 's/.*://' | sed 's/,.*//')
echo -e "${hog}\t${k}\t${p}" >> results.tsv
else
echo "$json entries not as expected"
fi
done

## Add a column for the result
cat results.tsv | \
awk '
NR == 1 {print $0"\tresult"; next};
{ printf $0"\t" };
$3 > 0.05 { print "not_significant"; next };
$2 < 1 { print "relaxation"; next };
$2 > 1 { print "intensification"; next };
$2 == 1 { print "equal_intensity"; next};
' \
> tmp
mv tmp results.tsv

##Add a line for each of the hogs that did not converge
ls unconverged_relax_out_second_attempt | sed 's/.relax_out.txt//' | awk '{print $0"\t-\t-\t-"}' > tmp
cat results.tsv tmp > tmp.tmp
mv tmp.tmp results.tsv
rm tmp

## p-values are for each individual test.
## To correct for running thousands of tests,
## calculate FDR q-values!
## add a column to the header!
mv results.tsv results_no_FDR.tsv
head -n 1 results_no_FDR.tsv | \
awk '{print $1"\t"$2"\t"$3"\tq-val"}' \
> results.tsv

## Calculate q-vals
m=$(( $(wc -l < results_no_FDR.tsv) -1 ))
cat results_no_FDR.tsv | \
tail -n +2 | \
awk '$2!="-"' | \
sort -k3,3g \
> tmp
m=$(wc -l < tmp)
cat tmp | \
awk -v m=$m '{print $1"\t"$2"\t"$3"\t"$3*m/NR}' \
>> results.tsv
rm tmp

## add back the lines for hogs that did not converge
cat results_no_FDR.tsv | \
tail -n +2 | \
awk '$2=="-" {print $0}' \
>> results.tsv

## Give a ternary result using q value.
## Run using different q values in ~/data/main_pipeline/hyphy_out/tables_with_different_FDR_cutoffs
## to see if there is any difference in the categorical (KEGG/COG) results
q=0.1 ## arbitrary. 0.1 is used by Jones et al, 2023, Nat. Eco. Evo.
echo "result" > trinary_results.tmp
tail -n +2 results.tsv | \
awk -v q=$q '
$4=="-" {print "-"; next};
$4>q {print "not_significant"; next};
$2<1 {print "relaxation"; next};
$2>1 {print "intensification"; next};
$2==1 {print "equal"; next};
' \
>> trinary_results.tmp
## Did any come back equal? Hopefully 0
awk '$4=="equal"' trinary_results.tmp | wc -l
## Paste
paste results.tsv trinary_results.tmp > paste.tmp
mv paste.tmp results.tsv
rm trinary_results.tmp


########################
### Further analyses ###
########################

## COG and KEGG analysis!
bash -i ~/scripts/helper/hyphy_relax_categorical_analysis.sh

############################
### Annotate results.tsv ###
############################

## Odd order, should probably be before hyphy_relax_categorical_analysis.sh
## But I don't think it will make the code any more or less easy to edit.

cd ~/data/main_pipeline/hyphy_out/relax_out
## Get the hogs in the same order as the annotations in the same as results.tsv
## Conveniently, both headers start with "HOG"
annotation_dir=~/data/main_pipeline/annotations/eggnog-mapper_out/hog_annotations/
> annotation_cols.tmp
cat results.tsv | \
awk -F "\t" '{print $1}' | \
while read hog; do
cat ${annotation_dir}/hog_annotations.tsv | \
awk -v hog=$hog -F "\t" ' $1==hog { printf $2; for (i=3;i<=NF;++i) {printf "\t"$i}; print "" } ' \
>> annotation_cols.tmp
done
## Paste onto results!
paste results.tsv annotation_cols.tmp > annotated_relax_results.tsv
rm annotation_cols.tmp


###########################
## What portion of HOGs ###
#### were included in #####
##### this analysis? ######
###########################

cd ~/data/main_pipeline/hyphy_out/relax_out
## How many HOGs from each genome?
cat ~/data/main_pipeline/orthofinder_out/for_downstream/Taz.n_intact.tsv | \
awk -F "\t" '
{ for (i=2;i<NF;++i) { printf $i"\t" }; print $NF }
' | \
awk -F "\t" '
NR==1 { for (i=1;i<=NF;++i) {A[i]=$i} ; next };
{ for (i=1;i<=NF;++i) { if ($i>0) {B[i]+=1} } };
END { for (i=1;i<=length(A);++i) { print A[i]"\t"B[i] } }
' | \
sort -k1,1 \
> Taz.n_HOGs.tsv

## How many loci from each genome?
cat ~/data/GC_content.tsv | \
grep -Ff ~/Taz.short_prefix.list | \
awk -F "\t" '{print $1"\t"$4}' | \
sort -k1,1 \
> Taz.n_intact.tsv

## How many loci from each genome used?
## Same as how many HOGs because only single copy used.
> loci.tmp
cat successfully_run_HOGs.list | \
sed 's/.relax_out.txt$//' | \
while read hog; do
if grep -q $hog ../msas/msas_with_exclamation_marks.list; then
msa=../msas/msas_with_exclamation_marks_removed/${hog}.noStop.fa
else
msa=../msas/${hog}.noStop.fa
fi
grep "^>" $msa | \
sed 's/^>//' \
>> loci.tmp
done
cat loci.tmp | sort | uniq -c | awk '{print $2"\t"$1}' | grep -Ff ~/Taz.short_prefix.list | sort -k1,1 > Taz.loci_used_counts.tsv
rm loci.tmp

## Paste together and get the percents
echo -e "Genome\tAnalyzed\ttotal_HOGs\ttotal_intact\tp_HOGs_analyzed\tp_intact_analyzed" > Taz.portion_loci_HOGs_analyzed.tsv
paste Taz.loci_used_counts.tsv Taz.n_HOGs.tsv Taz.n_intact.tsv | \
awk '{print $1"\t"$2"\t"$4"\t"$6"\t"($2/$4)"\t"($2/$6)}' \
>> Taz.portion_loci_HOGs_analyzed.tsv
rm Taz.n_HOGs.tsv Taz.n_intact.tsv Taz.loci_used_counts.tsv 
