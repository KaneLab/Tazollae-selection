## Run hyphy RELAX and FitMG94 on a Multiple Sequence Alignment (MSA)
## of all single-copy genes concatenated together

## Make an alignment of all the single-copy predicted-intact protein-coding sequences, concatenated
conda deactivate
msa_dir=~/data/main_pipeline/hyphy_out/msas
outdir=~/data/main_pipeline/hyphy_out/concatenated_msas

## Make an alignment of all the single-copy predicted-intact protein-coding sequences, concatenated
## To get a general sense of dN, dS, and dN/dS
mkdir $outdir
cd $outdir

## Gather all of the individual gene MSAs into a single file
msa_file=single_copy_HOGs.fa
> $msa_file
grep -vFf ~/data/main_pipeline/hyphy_out/msas/msas_with_exclamation_marks.list \
~/data/main_pipeline/orthofinder_out/for_downstream/single_copy_genes.list | \
while read hog; do
cat ~/data/main_pipeline/hyphy_out/msas/${hog}.noStop.fa | \
## Remove all but genome name from each header
sed 's/_.*//' | \
## Add the hog_id which will make sure the hogs stay in order in the below processing
sed "s/\(>.*\)/\1\t$hog/" \
>> $msa_file
done

## Make one header and one seq line per genome
cat $msa_file | \
## Put headers and seqs on same line
tr "\n" "\t" | \
sed 's/\t$//' | \
awk '{print $0}' | \
sed 's/>/\n>/g' | \
tail -n +2 | \
## sort to get all of the seqs from the same genome next to one another
sort -k1,1 -k2,2 | \
## put all seqs from the same genome on a single line
## And split header from seq line
awk -F "\t" '
BEGIN {prev=""};
$1==prev {printf $3; next};
{ prev=$1; print ""; print $1; printf $3 }
END {print ""}
' | \
tail -n +2 \
> tmp
## Check that the file looks correct
awk '$0 !~ "^>" {print length($0)}' tmp | sort -u
# Great! All the concatenated sequences are the same length
grep -c "^>" tmp
# 48 great.
grep -c "!" tmp
# 0 great
mv tmp $msa_file

## Prune the species tree to contain only the subtree for N9
in_tree=~/data/main_pipeline/orthofinder_out/final_run/Species_Tree/SpeciesTree_rooted_node_labels.txt
tree_file=single_copy_HOGs.pruned_tree.txt
conda deactivate
conda activate gotree
gotree prune -r -i $in_tree -f ~/subtree.prefix.list -o $tree_file
conda deactivate

## Label the test and reference groups.
## Even though this script is more interested in the "General Descriptive" model, test and reference groups are still required
## All nodes that contain exclusively T. az. are test except the root of the T. az. zlade which is unlabeled
## All other nodes are labeled reference
conda deactivate
logfile=log.tmp
> $logfile
in_tree=$tree_file
out_tree=${tree_file%.pruned_tree.txt}.labeled_tree.txt
bash -i ~/scripts/helper/label_trees_not_ingroup_root.sh -i $in_tree -o $out_tree -l $logfile\
 -I ~/Taz.short_prefix.list -O ~/subtree_freeLiving.short_prefix.list
## Check one set of erros from helper script
grep --no-group-separator -A1 "Are descendent nodes always before ancestor nodes?" $logfile | grep -v "^Are" | sort | uniq -c
## All yes
## Get the no HOGs
grep -B2 --no-group-separator "^No!" $logfile | awk 'NR %3==1'
## Empty
## logfile looks good
## Don't remove the pruned, unlabeled tree. Will be useful for figures.
tree_file=$out_tree

## Run RELAX on the concatenated alignment/tree for all single-copy genes
conda deactivate
outfile=single_copy_HOGs.relax_out.json
stdoutfile=single_copy_HOGs.relax_out.txt
## Run without srv flag for time.
hyphy relax --alignment $msa_file --tree $tree_file --test "test" --reference "reference" --output $outfile > $stdoutfile
mkdir relax_no_flags
mv single_copy_HOGs.relax* relax_no_flags
## Run again with srv
conda deactivate
outfile=single_copy_HOGs.relax_out.srv.json
stdoutfile=single_copy_HOGs.relax_out.srv.txt
hyphy relax --alignment $msa_file --tree $tree_file --test "test" --reference "reference" --srv Yes --starting-points 10 --output $outfile > $stdoutfile

## Make a phylogenetic tree colored by the k-value of the branch
## Similar to Fig 6 from Wertheim et al, 2014
## Get the k-values from the .json
jq -r '(.["branch attributes"]["0"] | to_entries | map ([.key,.value["k (general descriptive)"]])[]) | @tsv ' single_copy_HOGs.relax_out.json | \
sort -k2,2g \
> single_copy_HOGs.relax_k.tsv
## Normalize k-values to the geometric mean. Almost all of the k-values are <1, which is okay because the general descriptive
## model does not have a clearly defined outgroup to set k=1 (https://github.com/veg/hyphy/issues/1507)
k_sum=$(cat single_copy_HOGs.relax_k.tsv | awk 'BEGIN {c=1}; {c=c*$2}; END {print c^(1/FNR)}')
cat single_copy_HOGs.relax_k.tsv | awk -v k_sum=$k_sum '{print $1"\t"$2/k_sum}' > single_copy_HOGs.relax_k_normalized.tsv
## Modify the species tree to contain the k values and a binary yes/no T. az.
## To be used for visualizations on iTOL
cp single_copy_HOGs.pruned_tree.txt tree.tmp
cat single_copy_HOGs.relax_k_normalized.tsv | \
while read node k; do
sed -i "s/${node}:/${node}\[\&k=${k}\]:/" tree.tmp
done
cat ~/Taz.short_prefix.list | \
while read prefix; do
sed -i "s/\(${prefix}\[[^]]*\)\]/\1,endo=TRUE\]/" tree.tmp
done
mv tree.tmp single_copy_HOGs_with_k_tree.txt

###################
### Run FitMG94 ###
###################

outfile=single_copy_HOGs.FitMG94_out.json
stdoutfile=single_copy_HOGs.FitMG94_out.txt
hyphy ~/tools/hyphy/hyphy-analyses/FitMG94/FitMG94.bf --alignment $msa_file --tree $tree_file --type local --output $outfile > $stdoutfile
## Try again with --srv Yes --starting-points 10 
outfile=single_copy_HOGs.FitMG94_out.srv.json
stdoutfile=single_copy_HOGs.FitMG94_out.srv.txt
hyphy ~/tools/hyphy/hyphy-analyses/FitMG94/FitMG94.bf --alignment $msa_file --tree $tree_file --type local  --srv Yes --starting-points 10 --output $outfile > $stdoutfile


## Get the dN/dS values in a 2-col tsv sorted by value
jq -r '(.["branch attributes"]["0"] | to_entries | map ([.key,.value["Confidence Intervals"].MLE])[]) | @tsv ' $outfile | \
sort -k2,2g \
> single_copy_HOGs.FitMG94_omega.txt
## Get the dS values in a 2-col tsv sorted by value
jq -r '(.["branch attributes"]["0"] | to_entries | map ([.key,.value["dS"]])[]) | @tsv ' $outfile | \
sort -k2,2g \
> single_copy_HOGs.FitMG94_ds.txt
## Get the dN values in a 2-col tsv sorted by value
jq -r '(.["branch attributes"]["0"] | to_entries | map ([.key,.value["dN"]])[]) | @tsv ' $outfile | \
sort -k2,2g \
> single_copy_HOGs.FitMG94_dn.txt
## Get the synonymous, nonsynonymous, and combined trees
grep -A1 "^### \*\*Synonymous tree\*\*" $stdoutfile | tail -n 1 > synonymous_tree.txt
grep -A1 "^### \*\*Non-synonymous tree\*\*" $stdoutfile | tail -n 1 > nonsynonymous_tree.txt
grep -A1 "^**Combined tree**" $stdoutfile | tail -n 1 > combined_tree.txt

## Comparing the outputs with and without srv and starting points, they are almost identical
## Will continue with the outputs without those flags because I am not sure what those flags do for this algorithm
## Much less than a 1% change in each value
mkdir FitMG94_with_flags
mv single_copy_HOGs.FitMG94_out.srv.json single_copy_HOGs.FitMG94_out.srv.txt FitMG94_with_flags
## Modify the species tree to contain the dnds values and a binary yes/no T. az.
## To be used for visualizations on iTOL
## Get a list of the nodes with dN/dS = 0 and change to dN/dS = 1e-10 for iTOL visualization purposes.
cat single_copy_HOGs.FitMG94_omega.txt | \
awk '$2==0 {print $1"\t1e-10"; next}; {print $0}' \
> tmp
mv tmp single_copy_HOGs.FitMG94_omega.txt
cp single_copy_HOGs.pruned_tree.txt tree.tmp
cat single_copy_HOGs.FitMG94_omega.txt | \
while read node dnds; do
sed -i "s/${node}:/${node}\[\&dnds=${dnds}\]:/" tree.tmp
done
cat ~/Taz.short_prefix.list | \
while read prefix; do
sed -i "s/\(${prefix}\[[^]]*\)\]/\1,endo=TRUE\]/" tree.tmp
done
mv tree.tmp single_copy_HOGs_with_dnds_tree.txt

## generate tsv for graphing dn, ds, omega against each other
cat single_copy_HOGs.FitMG94_ds.txt | sort -k1,1 > ds.tmp
cat single_copy_HOGs.FitMG94_dn.txt | sort -k1,1 > dn.tmp
cat single_copy_HOGs.FitMG94_omega.txt | sort -k1,1 > omega.tmp
## Check genomes are in the same order!
paste ds.tmp dn.tmp omega.tmp | awk '$1!=$3 || $1!=$5'
## empty, good!
echo -e "clade\tgenome\tdS\tdN\tomega" > ds_dn_omega.tsv
paste ds.tmp dn.tmp omega.tmp | \
awk '{print $1"\t"$2"\t"$4"\t"$6}' \
>> body.tmp
rm ds.tmp dn.tmp omega.tmp
## Add the group to the body
echo "Taz,clade_I,clade_II" | \
sed 's/,/\n/g' | \
while read clade; do
## loop instead of grep -Ff to ensure perfect matches
cat ${clade}.nodes.list | \
while read node; do
cat body.tmp | \
awk -v node=$node -v clade=$clade ' $1==node { print clade"\t"$0 } ' \
>> ds_dn_omega.tsv
## Finish while node loop
done
## Finish while group loop
done
## Clean up
rm body.tmp 


## scp single_copy_HOGs.relax_k.tsv single_copy_HOGs.FitMG94_omega.txt single_copy_HOGs.labeled_tree.txt to local and use iTOL to make figures


####################
### ML + BS tree ###
####################

## Confirm tree topology using the concatenated single copy fasta
## This is basically whato rthofinder does
## But orthofinder casts a wider net, allowing a small portion of genomes to be missing a given gene
cd ~/data/main_pipeline/hyphy_out/concatenated_msas
conda activate orthofinder
## https://isu-molphyl.github.io/EEOB563/computer_labs/lab4/raxml-ng.html
raxml-ng --all --msa single_copy_HOGs.fa --model GTR+G --prefix T15 --seed 2 --threads 8 --bs-metric fbp,tbe
conda deactivate
