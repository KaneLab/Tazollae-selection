main_dir=/home/liam/data/main_pipeline/gene_loss
binary_table=${main_dir}/all.intact_vs_not.tsv
in_tree=/home/liam/data/main_pipeline/hyphy_out/concatenated_msas/single_copy_HOGs.pruned_tree.txt
outdir=${main_dir}/full_tree_single_gene_gain
cd $outdir


## Initiate table that will track the LCA of each HOG (and thus where that HOG was gained)
> hog_lca.tsv
## initiate table that will track likelihood of positive at each node for each HOG
outtable=positive_probabilities.tsv
echo "Node" > $outtable
cat all_nodes.list | sort >> $outtable

conda deactivate
conda activate gotree
## Loop through cols in table, each of which is a HOG
## Could do n_cols=$( cat trimmed_transposed_table.tsv | head -n 1 | awk '{print NF}' )
## But want to check that all lines have the same NF.
n_cols=$( cat trimmed_transposed_table.tsv | awk '{print NF}' | sort -u )
col=2
while (( $col <= $n_cols )); do

echo $col
## Subset table
cat trimmed_transposed_table.tsv | \
awk -v col=$col '{print $1"\t"$col}' \
> this_hog.tsv
## Get the hog
this_hog=$( cat this_hog.tsv | head -n 1 | awk '{print $2}' )
## And remove the header from this_hog.tsv
tail -n +2 this_hog.tsv > tail.tmp
mv tail.tmp this_hog.tsv

## Get list of genomes with this HOG
cat this_hog.tsv | \
awk '$2==1 {print $1}' \
> positive_genomes.list

## If positive_genomes.list is only one genome, then don't need the analysis below
if (( $( wc -l < positive_genomes.list ) == 1 )); then
cat positive_genomes.list | sed 's/$/\t1/' > ape_ace_out.tsv
cat positive_genomes.list | sed "s/^/${this_hog}\t/" >> hog_lca.tsv

else


## Get the subtree 
gotree prune -r -i $in_tree -f positive_genomes.list -o pruned_tree.txt

## Get the LCA of this subtree
LCA=$( cat pruned_tree.txt | sed 's/;$//' | sed 's/.*)//' )
## Add LCA to gained list
echo -e "${this_hog}\t${LCA}" >> hog_lca.tsv
rm pruned_tree.txt

## Get the subtree rooted at the LCA:
gotree subtree -i $in_tree -n "^$LCA" -o subtree.txt

## Get the list of leaves in subtree.txt that are negative
gotree labels -i subtree.txt | \
grep -vFf positive_genomes.list | \
## Create a file of all the subtree leaves with 0 (negative) or 1 (positive)
sed 's/$/\t0/' > subtree_tip_binary.tsv
## Add positive tips
cat positive_genomes.list | \
sed 's/$/\t1/' \
>> subtree_tip_binary.tsv

## If all leaves in the subtree are positive, can't call R:APE:ACE and no need.
n_zeros=$(cat subtree_tip_binary.tsv | awk '$2==0' | wc -l)
if (( $n_zeros == 0 )); then
cp subtree_tip_binary.tsv ape_ace_out.tsv
else
## Run R:APE:ACE
Rscript --vanilla /home/liam/scripts/helper/ape_ace_discrete_wrapper_single_hog.R table_path=subtree_tip_binary.tsv tree_path=subtree.txt outfile=ape_ace_out.tsv
fi

## finish else block of if (( $( wc -l < positive_genomes.list ) == 1 )); then
fi

## Add 0's for each tip/node that isn't in the subtree
cat ape_ace_out.tsv | awk '{print $1}' > in_nodes.tmp
grep -vFf in_nodes.tmp all_nodes.list | sed 's/$/\t0/' >> ape_ace_out.tsv
rm in_nodes.tmp

## Prepare output
echo ${this_hog} > this_out_col.tmp
## Sort nodes so in the same order for each HOG
cat ape_ace_out.tsv | \
sort -k1,1 | \
## Keep only the values
awk '{print $2}' \
>> this_out_col.tmp

paste $outtable this_out_col.tmp > paste.tmp
mv paste.tmp $outtable

## Iterate while loop
col=$(( $col + 1 ))
## Clean up
rm ape_ace_out.tsv subtree.txt this_out_col.tmp subtree_tip_binary.tsv positive_genomes.list this_hog.tsv

## finish while $(( $col <= $n_cols )); do
done
conda deactivate
