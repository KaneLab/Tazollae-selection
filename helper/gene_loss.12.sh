
##########################
### Gain/Loss Analyses ###
##########################

## Run on T. az. subtree with only loss allowed
main_dir=~/data/main_pipeline/gene_loss
binary_table=${main_dir}/Taz.intact_vs_not.tsv
in_tree=~/data/main_pipeline/orthofinder_out/final_run/Species_Tree/SpeciesTree_rooted_node_labels.txt
subtree_list=~/Taz.short_prefix.list
outdir=Taz_subtree_no_gene_gain
approach=no_gain

bash -i ~/scripts/helper/gene_loss_helper.sh -m $main_dir -b $binary_table\
 -t $in_tree -a $approach -l $subtree_list -o $outdir


## Run on whole tree with gain and loss allowed
main_dir=~/data/main_pipeline/gene_loss
binary_table=${main_dir}/all.intact_vs_not.tsv
in_tree=~/data/main_pipeline/orthofinder_out/final_run/Species_Tree/SpeciesTree_rooted_node_labels.txt
subtree_list=~/subtree.prefix.list
outdir=full_subtree_yes_gene_gain
approach=ARD

bash -i ~/scripts/helper/gene_loss_helper.sh -m $main_dir -b $binary_table\
 -t $in_tree -a $approach -l $subtree_list -o $outdir

## Didn't do the following because I don't think it would be helpful:
## Run on Taz-less sub-tree with gain and loss allowed
#main_dir=~/data/main_pipeline/gene_loss
#binary_table=${main_dir}/all.intact_vs_not.tsv
#in_tree=~/data/main_pipeline/orthofinder_out/final_run/Species_Tree/SpeciesTree_rooted_node_labels.txt
#subtree_list=~/subtree.prefix.list
#outdir=full_subtree_yes_gene_gain
#approach=ARD

#bash -i ~/scripts/helper/gene_loss_helper.sh -m $main_dir -b $binary_table\
# -t $in_tree -a $approach -l $subtree_list -o $outdir

############################
### Combine results from ###
#### multiple approaches ###
############################

cd $main_dir

## Replace T. az. values from the larger subtree analysis with those from the Taz-only tree analysis
## Get the edge list
head -n -1 Taz_subtree_no_gene_gain/ancestors.tsv | awk -F "\t" '{print $1}' > Taz_edges.tmp
## Replace values in the cumulative gain/loss files
echo "loss,gain" | \
sed 's/,/\n/' | \
while read change; do
file=cumulative_gene_${change}_probability.tsv
cp Taz_subtree_no_gene_gain/$file ./
grep -vFf Taz_edges.tmp full_subtree_yes_gene_gain/$file >> $file
done
rm Taz_edges.tmp

## Remake Taz_edges.tmp including N35
## Was not included above because I wanted the full_subtree values for gene loss
cat Taz_subtree_no_gene_gain/ancestors.tsv | awk -F "\t" '{print $1}' > Taz_edges.tmp
## But I want the gene count from the Taz_subtree...
## Combine values in the ape_ace_out.tsvs
## The taz_subtree runs did not include HOGs that were unanimously present or absent within the taz.
## And the full_subtree runs did not include HOGs that were unanimously present
## Therefore, add these HOGs.
## How many fields in ape_ace_out.tsv?
echo "Taz_subtree_no_gene_gain
full_subtree_yes_gene_gain" | \
while read subtree; do
tmp_file=${subtree%%_*}.tmp
num_fields=$(cat ${subtree}/ape_ace_out.tsv | awk '{print NF}' | sort -u)
## Check that it is just one value. success.
## Add rows of 1's and 0's for unanimously present and absent HOGs, respectively
cp ${subtree}/ape_ace_out.tsv $tmp_file
cat ${subtree}/unanimous_hogs.list | \
awk -v n=$num_fields ' { printf $1; for (i=2;i<=n;++i) {printf "\t1"}; print "" } ' \
>> $tmp_file
cat ${subtree}/absent_hogs.list | \
awk -v n=$num_fields ' { printf $1; for (i=2;i<=n;++i) {printf "\t0"}; print "" } ' \
>> $tmp_file
## Finish subtree loop
## Put hogs in order
sort -k1,1 $tmp_file > ${tmp_file}.tmp
mv ${tmp_file}.tmp ${tmp_file}
done
## ensured with awk '{print $1}' and cmp that the hog column matches perfectly full.tmp and Taz.tmp

## Now replace the Taz_edges from full.tmp with those from Taz.tmp
## paste each column not from the Taz_subtree onto taz.tmp, one at a time
cat full.tmp | \
head -n 1 | \
sed 's/\t/\n/g' | \
tail -n +2 | \
grep -vFf Taz_edges.tmp | \
while read node; do
cat full.tmp | \
awk -v node=$node '
NR==1 { for (i=1;i<=NF;++i) { if ($i==node) {c=i} } };
{print $c}
' \
> this_col.tmp
paste Taz.tmp this_col.tmp > Taz.tmp.tmp
mv Taz.tmp.tmp Taz.tmp
rm this_col.tmp
## finish node loop
done
mv Taz.tmp combined_ace_ape_out.tsv
rm Taz_edges.tmp full.tmp


## Normalize the cumulative probabilities to the branch lengths
## Need to do this again instead of just using the normalized values from the runs of the helper scirpt
## Because those were normalized separately and thus meaningless to combine as is.
## Put nodes in same order for cumulative p of gene loss and edge length
echo "loss,gain" | \
sed 's/,/\n/' | \
while read change; do
sort -k1,1 cumulative_gene_${change}_probability.tsv > p.tmp
sort -k1,1 full_subtree_yes_gene_gain/edge_lengths.tsv > l.tmp
## Normalize cum_p and branch length so each has an average of 1
echo "p.tmp
l.tmp" | \
while read infile; do
c_ave=$( cat $infile | awk 'BEGIN {c=0}; {c+=$2}; END {print c/FNR} ' )
cat $infile | awk -v c_ave=$c_ave '{print $1"\t"$2/c_ave}' > ${infile}.tmp
mv ${infile}.tmp $infile
done
## Divide cum_p by branch length!
paste p.tmp l.tmp | \
awk ' $1==$3 { print $1"\t"$2/$4 } ' \
> rate_of_cumulative_gene_${change}_probability.tsv
## Log the rate!
cat rate_of_cumulative_gene_${change}_probability.tsv | awk '{print $1"\t"log($2)}' > logged_rate_of_cumulative_gene_${change}_probability.tsv
## Check all rows of input received an output
echo "all the following should be equal:"
wc -l < p.tmp; wc -l < l.tmp; wc -l < rate_of_cumulative_gene_${change}_probability.tsv
## All equal. Success!
rm p.tmp l.tmp
done

## Now also normalize for the number of HOGs predicted at the ancestral node
## Look at the number of HOGs predicted possibly intact in the various genomes and internal nodes
cat combined_ace_ape_out.tsv | \
awk '
NR==1 { for (i=2;i<=NF;++i) { N[i-1]=$i; A[i-1]=0 } };
{ for (i=2;i<=NF;++i) { A[i-1]+=$i } };
END { for (i=1;i<=length(A);++i) {print N[i]"\t"A[i]"\t"} }
' | \
sort -k2,2g \
> n_HOGs_intact.tsv
## Can also use n_HOGs_intact.tsv for
## Using the HOGs (as opposed to ORFs and intact ORFs as in GC_content.tsv)
## to recreate genome size in internal nodes.
## Find the direct ancestral node of each node and then its predicted # intact HOGs
> n_ancestral_HOGs.tmp
echo "loss
gain" | \
while read change; do
cat rate_of_cumulative_gene_${change}_probability.tsv | \
awk '{print $1}' | \
while read node; do
ancestor=$(grep "^$node" full_subtree_yes_gene_gain/ancestors.tsv | sed 's/.*\t//' | sed 's/,.*//')
grep "^$ancestor" n_HOGs_intact.tsv | awk '{print $2}' >> n_ancestral_HOGs.tmp
## finish node loop
done
## Get the average number of HOGs, since the rate_of_cumulative_gene_${change}_probability.tsv  values are already normalized to 1.
n_hog_ave=$(awk ' BEGIN {c=0}; {c+=$0} ; END {print c/FNR}' n_ancestral_HOGs.tmp)
paste rate_of_cumulative_gene_${change}_probability.tsv n_ancestral_HOGs.tmp | \
awk -v a=$n_hog_ave '{print $1"\t"$2/($3/a)}' \
> gene_${change}_norm_by_branch_length_ancestral_HOG_count.tsv
## Log it!
cat gene_${change}_norm_by_branch_length_ancestral_HOG_count.tsv | awk '{print $1"\t"log($2)}' > log_gene_${change}_norm_by_branch_length_ancestral_HOG_count.tsv
rm n_ancestral_HOGs.tmp
## finish change loop
done

## Make metadata doc for iTOL
mkdir iTOL_metadata
echo "loss,gain" | \
sed 's/,/\n/' | \
while read change; do
cat ~/iTol_metadata_template.txt cumulative_gene_${change}_probability.tsv > iTOL_metadata/cum_p_gene_${change}_iTOL_metadata.txt
cat ~/iTol_metadata_template.txt rate_of_cumulative_gene_${change}_probability.tsv > iTOL_metadata/rate_of_cum_p_gene_${change}_iTOL_metadata.txt
cat ~/iTol_metadata_template.txt logged_rate_of_cumulative_gene_${change}_probability.tsv > iTOL_metadata/log_rate_of_cum_p_gene_${change}_iTOL_metadata.txt
cat ~/iTol_metadata_template.txt log_gene_${change}_norm_by_branch_length_ancestral_HOG_count.tsv > iTOL_metadata/log_gene_${change}_norm_by_branch_length_ancestral_HOG_count_iTOL_metadata.txt
done
cat ~/iTol_metadata_template.txt n_HOGs_intact.tsv > iTOL_metadata/genome_size_intact_HOGs_iTOL_metadata.txt

########################
### Jaccard Distance ###
########################

## Jaccard distance based on binary tables
cd ~/data/main_pipeline/gene_loss

## Do both intact/not and present/absent
echo "all.intact_vs_not.tsv,jaccard_distance_intact.tsv
all.presence_absence.tsv,jaccard_distance_present.tsv" | \
sed 's/,/\t/g' | \
while read intable outtable; do

## Initiate table:
cat ~/subtree.prefix.list | tr "\n" "\t" | sed 's/\t$//' | awk '{print $0}' > ${outtable}
## Loop through every pair of genomes
## This will be a full square matrix.
## Could limit it to a triangle, but should run quickly regardless
cat ~/subtree.prefix.list | \
while read prefix_1; do
echo $prefix_1 > this_col.tmp
cat ~/subtree.prefix.list | \
while read prefix_2; do

cat ${intable} | \
awk -v prefix_1=$prefix_1 -v prefix_2=$prefix_2 -F "\t" '
BEGIN { intersection=0; union=0 };
NR==1 { for (i=1;i<=NF;++i) { if ($i==prefix_1) {a=i}; if ($i==prefix_2) {b=i}; }; next };
$a==1 || $b==1 { union+=1 };
$a==1 && $b==1 { intersection+=1 };
END { print 1-intersection/union };
' \
>> this_col.tmp

## Finish prefix_2 loop
done
cat this_col.tmp | tr "\n" "\t" | sed 's/\t$//' | awk '{print $0}' >> ${outtable}
rm this_col.tmp
## Finish prefix_1 loop
done

## finish intable outtable loop
done

## Visualize jaccard distances with trees and PCoAs using script scripts/visualizations/jaccard_visualizations.R
## Paths in that file are for my local computer

################################
### Phylogenetic and Jaccard ###
##### Distance Correlation #####
################################

## Get a matrix of pairwise distances between genomes on the phylogenetic tree
conda activate gotree
in_tree=~/data/main_pipeline/hyphy_out/concatenated_msas/single_copy_HOGs.pruned_tree.txt 
gotree matrix -i $in_tree -o phylogenetic_distance_matrix.tsv
tail -n +2 phylogenetic_distance_matrix.tsv > tail.tmp
mv tail.tmp phylogenetic_distance_matrix.tsv
conda deactivate

binary=intact

## Get phylogenetic distance matrix in same order as jaccard distance matrix
awk '{print $1}' phylogenetic_distance_matrix.tsv > old_order.tmp
awk '{print $1}' jaccard_distance_${binary}.tsv | tail -n +2 > new_order.tmp
## First, get the rows in the same order.
> rows_sorted.tmp
cat new_order.tmp | \
while read prefix; do
cat phylogenetic_distance_matrix.tsv | \
awk -v prefix=$prefix '$1==prefix' \
>> rows_sorted.tmp
done
## Now get the columns in the same order
cp new_order.tmp fully_sorted.tmp
cat new_order.tmp | \
while read prefix; do
n=$( cat old_order.tmp | awk -v prefix=$prefix '$0==prefix {print FNR}' )
cat rows_sorted.tmp | awk -v n=$n '{print $(n+1)}' > this_col.tmp
paste fully_sorted.tmp this_col.tmp > paste.tmp
mv paste.tmp fully_sorted.tmp
rm this_col.tmp
done
## Add header and rename
cat new_order.tmp | tr "\n" "\t" | sed 's/\t$//' | awk '{print $0}' > phylogenetic_distance_matrix.tsv
cat fully_sorted.tmp >> phylogenetic_distance_matrix.tsv
rm fully_sorted.tmp rows_sorted.tmp new_order.tmp old_order.tmp

## Make long format of triangularized jaccard and phylogenetic distance matrices
## Loop inputs
echo "jaccard_distance_intact.tsv
jaccard_distance_present.tsv
phylogenetic_distance_matrix.tsv" | \
while read infile; do

cat $infile | \
awk -F "\t" '
## Get the comparison genome order from the header
NR==1 { for (i=1;i<=NF;++i) { A[i]=$i }; next };
## Loop through genomes that are in an earlier position than the genome for this row
## So dont get converse comparisons nor self comparisons
{ for (i=2;i<FNR;++i) { print $1"\t"A[i-1]"\t"$i } };
' \
> ${infile%.tsv}_long.tsv

## Finish while read infile loop
done

## double check all in the same order
cat jaccard_distance_intact_long.tsv | awk '{print $1"\t"$2}' > tmp
cat jaccard_distance_present_long.tsv | awk '{print $1"\t"$2}' > tmp.tmp
cat phylogenetic_distance_matrix_long.tsv | awk '{print $1"\t"$2}' > tmp.tmp.tmp
cmp tmp tmp.tmp
cmp tmp.tmp tmp.tmp.tmp
rm tmp tmp.tmp tmp.tmp.tmp
## All the same!
## Replace the square matrices which are not useful!
echo "jaccard_distance_intact_long.tsv
jaccard_distance_present_long.tsv
phylogenetic_distance_matrix_long.tsv" | \
while read file; do
mv $file ${file%_long.tsv}.tsv
done
mv phylogenetic_distance_matrix.tsv phylogenetic_distance.tsv

## Combine so can graph jaccard distance against phylogenetic distance
echo "intact
present" | \
while read binary; do
echo -e "Genome_1\tGenome_2\tPhylogenetic\tJaccard" > phylo_jaccard_distances_xy_${binary}.tsv
paste phylogenetic_distance.tsv jaccard_distance_${binary}.tsv | \
awk '{print $1"\t"$2"\t"$3"\t"$6}' \
>> phylo_jaccard_distances_xy_${binary}.tsv
## finish while read binary loop
done

## Label the T. az.
echo "intact
present" | \
while read binary; do
infile=phylo_jaccard_distances_xy_${binary}.tsv
head -n 1 $infile | sed 's/^/Group\t/' > label.tmp
tail -n +2 $infile | grep -Ff ~/Taz.short_prefix.list | sed 's/^/Taz\t/' >> label.tmp
tail -n +2 $infile | grep -vFf ~/Taz.short_prefix.list | sed 's/^/Other\t/' >> label.tmp
mv label.tmp $infile
## finish while read binary loop
done


## Graph in R!
