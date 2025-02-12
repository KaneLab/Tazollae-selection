## Run hyphy RELAX on individual genes, the ingroup and outgroups are variable.
## This is almost identical to hyphy_relax.9.sh, and should replace it if this pipeline were run again

## Set variables
main_dir=~/data/main_pipeline/relax_control
orthofinder_outdir=~/data/main_pipeline/orthofinder_out/for_downstream
ingroup_nickname=Nos2
main_dir=${main_dir}/$ingroup_nickname
##Create a central directory for these analyses:
mkdir ${main_dir}
mkdir ${main_dir}/hyphy_out
mkdir ${main_dir}/macse_out
## Create prefix lists for this analysis
cd ${main_dir}
conda activate gotree
in_tree=~/data/main_pipeline/orthofinder_out/final_run/Species_Tree/SpeciesTree_rooted_node_labels.txt
gotree subtree -i $in_tree -n "^N12$" | \
gotree labels \
> ingroup.short_prefix.list
conda deactivate
grep -vFf ingroup.short_prefix.list ~/subtree.prefix.list > outgroup.short_prefix.list
grep -Ff ingroup.short_prefix.list ~/subtree.shortened.tsv > ingroup.shortened.tsv
grep -Ff outgroup.short_prefix.list ~/subtree.shortened.tsv > outgroup.shortened.tsv

##########################
### Generate hog lists ###
##########################

## Split n_intact.tsv into outgroup and ingroup
infile=${orthofinder_outdir}/n_intact.tsv
echo -e "\
ingroup,ingroup.shortened.tsv
outgroup,outgroup.shortened.tsv\
" | \
sed 's/,/\t/g' | \
while read outfile_prefix prefix_list; do

outfile=${outfile_prefix}.n_intact.tsv
cat $infile | awk -F "\t" '{print $1}' > $outfile

cat $prefix_list | awk '{print $2}' | \
while read prefix; do

cat $infile | \
awk -v prefix=$prefix -F "\t" '
NR==1 { for (i=4;i<=NF;++i) { if ($i==prefix) {c=i} } };
{ print $c }
' \
> this_col.tmp
paste $outfile this_col.tmp > tmp
mv tmp $outfile
rm this_col.tmp

## Finish loop through prefixes
done
## Finish loop through ingroup / outgroup
done

## Get the list of hogs with at least 3 genomes from each group (in/out) having exactly one locus
## And get the list of genomes in those sets.
## First, get a list of genomes with a single copy in each hog.
echo "\
ingroup
outgroup\
" | \
while read file_prefix; do
file=${file_prefix}.n_intact.tsv
## for each hog, get the list of genomes with exactly one locus!
cat $file | \
awk -F "\t" '
NR==1 { for (i=1;i<=NF;++i) {A[i]= $i}; next };
{
        printf $1;
        for (i=2;i<=NF;++i) {if ($i==1) { printf "\t"A[i] } }
        print ""
}
' \
> ${file_prefix}.all_hogs_which_genomes_single_copy.tsv
## Finish loop ingroup vs outgroup
done

## Get all of the HOGs with at least 3 genomes from each of ingroup and outgroup
cat ingroup.all_hogs_which_genomes_single_copy.tsv outgroup.all_hogs_which_genomes_single_copy.tsv | \
awk -F "\t" 'NF>=4 {print $1}' | \
sort | \
uniq -c | \
awk '$1==2 {print $2}' \
> hogs_passed_3_thresh.list

#####################
### Generate MSAs ###
#####################

##Run macse on the single-copy genes.
conda deactivate
conda activate macse
outdir=${main_dir}/macse_out
mkdir ${outdir}/trimmed_nucleotideHOGs
indir=${orthofinder_outdir}/nucleotide_HOGs
cd $outdir
cat ../hogs_passed_3_thresh.list | \
## subset for parralelization
awk 'NR%4==3' | \
while read hog; do
## Subset the .ffn to only those genomes that had exactly 1 locus in this hog
cat ${main_dir}/ingroup.all_hogs_which_genomes_single_copy.tsv | \
awk -v hog=$hog -F "\t" ' $1==hog { for (i=2;i<=NF;++i) {print $i} } ' > ${hog}.locus_list.tmp
cat ${main_dir}/outgroup.all_hogs_which_genomes_single_copy.tsv | \
awk -v hog=$hog -F "\t" ' $1==hog { for (i=2;i<=NF;++i) {print $i} } ' >> ${hog}.locus_list.tmp
grep --no-group-separator -A1 -Ff ${hog}.locus_list.tmp ${indir}/${hog}.ffn > ${outdir}/trimmed_nucleotideHOGs/${hog}.ffn
rm ${hog}.locus_list.tmp
## Some hogs fail to align, and don't align on rerun
## So would probably be good to have macse print to a logfile and parse that to see what missed.
## But there I don't know of anything to be done to improve performance, so no need at the moment.
macse -prog alignSequences -seq ${outdir}/trimmed_nucleotideHOGs/${hog}.ffn -out_AA ${hog}.macse.faa -out_NT ${hog}.macse.fna
done

## Remove the final stop codon from each of the sequences.
## BUSTED and RELAX don't allow the final stop codon, though nonREV appears to be fine with it.
## But will remove it up front for consistency.
## Replace with 3 dashes to keep alignment correct.
## Also change the header to only the locus tag because hyphy was getting confused by the multiple field header lines
## Also replace the locus_tags with just the species names because I am using the species tree, not gene trees.
mkdir ${main_dir}/hyphy_out/msas
cd ${main_dir}/hyphy_out/msas

ls ${main_dir}/macse_out/*.fna | \
while read msa; do
hog=${msa%.macse.fna}
hog=${hog##*/}
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
> ${hog}.noStop.fa
done


##########################
#### Prune gene trees ####
##########################

## I am using the species tree instead of gene trees (orthofinder gives both)
## because I trust the species tree more and we think the potential drawbacks of using the species tree
## (incomplete lineage sorting, HGT, to name 2) are unlikely in this system.
## Not all of these hogs contain a locus from every species
## So must prune the tree to match the hog.

msa_dir=${main_dir}/hyphy_out/msas
in_tree=~/data/main_pipeline/orthofinder_out/final_run/Species_Tree/SpeciesTree_rooted_node_labels.txt
tree_outdir=${main_dir}/hyphy_out/trees
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

#error log contains only the header and the name of each hog. Nothing else printed! success!
rm species2keep.tmp $logfile

###################
### Label trees ###
###################


## Label the test and reference groups.
## All nodes that contain exclusively T. az. are test except the root of the T. az. zlade which is unlabeled
## All other nodes are labeled reference
conda deactivate
indir=${main_dir}/hyphy_out/trees
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
 -I ${main_dir}/ingroup.short_prefix.list -O ${main_dir}/outgroup.short_prefix.list
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

## Run with test and background branches labeled.

## set variables and make folders
## Define main_dir and orthofinder_outdir again to make it easier to parallelize the loop
#main_dir=~/data/main_pipeline/relax_control
#orthofinder_outdir=~/data/main_pipeline/orthofinder_out/for_downstream
#ingroup_nickname=Nos_full
#main_dir=${main_dir}/$ingroup_nickname
conda deactivate
msa_dir=${main_dir}/hyphy_out/msas
tree_dir=${main_dir}/hyphy_out/trees
outdir=${main_dir}/hyphy_out/relax_out
mkdir $outdir
cd $outdir

ls ${tree_dir}/*.labeled_tree.txt | \
#grep -vFf already_run.tmp | \
## subset to run in parallel
awk 'NR%4==0' | \
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

###############################
### Re-run runs with errors ###
###############################

echo "Check for error messages:"
grep -i "error" *
## Some runs with errors can just be re-run with the normal call and will be fine.
## Having run ~2000 HOGs several times on the same HOGs, it seems that sometimes the omega
## estimates get too extreme and cause an error. But then this is fixed on re-run.
grep -i "error" * | \
sed 's/.relax_out.txt:Error:$//' | \
while read hog; do
tree_file=${tree_dir}/${hog}.labeled_tree.txt
msa_file=${msa_dir}/${hog}.noStop.fa
outfile=${hog}.relax_out.json
stdoutfile=${hog}.relax_out.txt
hyphy relax --alignment $msa_file --tree $tree_file --test "test" --reference "reference"\
 --models Minimal --srv Yes --output $outfile > $stdoutfile
done

## Repeat as necessary
echo "Check for error messages:"
grep -i "error" *

#############################
### Re-run RELAX for HOGs ###
### that did not converge ###
#############################

##Identify RELAX results that did not converge
cd ${main_dir}/hyphy_out/relax_out
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

##Check that I got all files with those two greps
wc -l unconverged_HOGs.list
# 54 (Cyl) 60 (DolAna) 91 (Nos) 106 (Nos_full) 84 (CalTol) 116 (Nos2)
wc -l converged_HOGs.list
# 1872 (Cyl) 3157 (DolAna) 4980 (Nos) 5727 (Nos_full) 4552 (CalTol) 5547 (Nos2)
ls *relax_out.txt | wc -l
# 1926 = 54 + 1872 (Cyl) 3217 = 60 + 3157 (DolAna) 5071 = 91 + 4980 (Nos)
# 5833 = 106 + 5727 (Nos_full) 4636 = 84 + 4552 (CalTol) 5663 = 5547 + 116 (Nos2)
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
## Define main_dir again to make parallelization easier
#main_dir=~/data/main_pipeline/relax_control
#ingroup_nickname=Nos_full
#main_dir=${main_dir}/$ingroup_nickname
#cd ${main_dir}/hyphy_out/relax_out
#msa_dir=${main_dir}/hyphy_out/msas
#tree_dir=${main_dir}/hyphy_out/trees
#Loop and run!> tmp
cat unconverged_HOGs.list | \
## Subset for parallelization
awk 'NR%4==2' | \
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

## Check for errors
echo "Check for error messages:"
grep -i "error" *

##Identify RELAX results that did not converge after second attempt
cd ${main_dir}/hyphy_out/relax_out
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
echo "Check for other problems:"
cat other_problems_HOGs.list

##Check that I got all files with those two greps
wc -l unconverged_HOGs.list
# 77 (Cyl) 38 (DolAna) 47 (Nos) 58 (Nos_full) 35 (CalTol) 62 (Nos2)
wc -l successfully_run_HOGs.list
# 2420 (Cyl) 3179 (DolAna) 5024 (Nos) 5775 (Nos_full) 4601 (CalTol) 5601 (Nos2)
wc -l other_problems_HOGs.list
# 0 0 0 0 0 0
ls *relax_out.txt | wc -l
# 2497 = 2420 + 77 + 0 (Cyl) 3217 = 3179 + 38 + 0 (DolAna) 5071 = 5024 + 47 (Nos)
# 5833 = 5775 + 58 (Nos_full) 4636 = 4601 + 35 (CalTol) 5663 = 5601 + 62 (Nos2)
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
 -m $main_dir -I ${main_dir}/ingroup.short_prefix.list -O ${main_dir}/outgroup.short_prefix.list

##########################
### Parse and Annotate ###
##### Relax Results ######
##########################

cd ${main_dir}/hyphy_out/relax_out

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
cat unconverged_HOGs.list | sed 's/.relax_out.txt//' | awk '{print $0"\t-\t-\t-"}' > tmp
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
m=$(wc -l < results_no_FDR.tsv)
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

## Give a trinary result using q value.
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


#########################
### After running all ###
#### of the ingroups ####
#### combine results ####
##### into one table ####
#########################

cd ~/data/main_pipeline/relax_control

## Get the number of hogs that fall into each category for each ingroup
## Not including hogs for which hyphy RELAX did not converge
## Loop 3 q-vals to see if it makes a difference
echo "0.1
0.01
0.001" | \
while read q_thresh; do
outfile=relax_w_diff_ingroups_q_${q_thresh}.tsv
echo "Ingroup (n_genes)
relaxation
not_significant
intensification" > $outfile
## Loop through ingroups
echo "Taz
DolAna
Cylindrospermopsis
Nos2
CalTol
Nos_full" | \
while read ingroup; do

## Remove header and non-convergent results
cat ${ingroup}/hyphy_out/relax_out/results.tsv | \
tail -n +2 | \
awk -F "\t" '$5!="-"' | \
## remove ternary result
awk '{print $1"\t"$2"\t"$3"\t"$4}' \
> results.tmp
## Count results for each ternary category
cat results.tmp |
awk -v ingroup=$ingroup -v q=$q_thresh '
BEGIN {r=0; n=0; i=0};
$4>q {n+=1; next};
$2<1 {r+=1; next};
{i+=1};
END {print ingroup" ("FNR")"; print r; print n; print i}
' > \
this_col.tmp
## paste!
paste $outfile this_col.tmp > paste.tmp
mv paste.tmp $outfile

## End loop through ingroup
done
rm results.tmp
rm this_col.tmp

## End q_thresh loop
done

## Statistical significance of T. az. percents for ternary results differing from those of other ingroups
## Will re-use the same first-col for several tables
echo "-
Cylindrospermopsis
DolAna
Nos
Taz" > first_col.tmp
## Loop 3 q-vals to see if it makes a difference
echo "0.1
0.01
0.001" | \
while read q_thresh; do
infile=relax_w_diff_ingroups_q_${q_thresh}.tsv
## Loop through the ternary results
echo "relaxation
not_significant
intensification" | \
while read ternary; do
## Initiate outfile
outfile=p_val_q_${q_thresh}_${ternary}.tsv
cat first_col.tmp > $outfile
## Loop through ingroups
echo "Cylindrospermopsis
DolAna
Nos
Taz" | \
while read ingroup; do
echo $ingroup > p_val_col.tmp
## Get a 3-col tsv with each row the 3 inputs to binom.test() for ingroup compared to a different ingroup
## Done such that percent comes from the other ingroup, while x and n come from this ingroup
cat $infile | \
awk -v ingroup=$ingroup -v ternary=$ternary -F "\t" '
NR==1 { for (i=2;i<=NF;++i) { if ($i~ingroup) {c=i}; split($i,A,"[[:punct:]]") ; N[i]=A[2] } };
$1==ternary { for (i=2;i<=NF;++i) { print $c"\t"N[c]"\t"$i/N[i] } }
' | \
while read x n p; do
## Call R binom.test()
Rscript --vanilla ~/scripts/helper/binom_test_p_val.R x=$x n=$n p=$p outfile=p_val_col.tmp
done

paste $outfile p_val_col.tmp > paste.tmp
mv paste.tmp $outfile
rm p_val_col.tmp
## End loop through ingroup
done
## End ternary loop
done
## End q_thresh loop
done
rm first_col.tmp

## Now make the tables percents of ternary results as opposed to counts
echo "0.1
0.01
0.001" | \
while read q_thresh; do
infile=relax_w_diff_ingroups_q_${q_thresh}.tsv
cat $infile | \
awk -F "\t" '
NR==1 { print $0; for (i=1;i<=NF;++i) { split($i,A,"[)(]") ; N[i]=A[2] }; next};
{ printf $1; for (i=2;i<=NF;++i) { printf "\t"$i/N[i] }; print "" };
' > percents.tmp
mv percents.tmp $infile
## Finish q_thresh loop
done
