## Get variables
print_usage() {
  printf "usage requires:\n\
main_dir -m\n\
binary_table -b\n\
in_tree -t\n\
approach -a\n\
subtree_list -l\n\
outdir -o\
"
}

while getopts 'i:l:I:O:o:' flag; do
  case "${flag}" in
    m) main_dir="${OPTARG}" ;;
    b) binary_table="${OPTARG}" ;;
    t) in_tree="${OPTARG}" ;;
    a) approach="${OPTARG}" ;;
    l) subtree_list="${OPTARG}" ;;
    o) outdir="${OPTARG}" ;;
    *) print_usage
       exit 1 ;;
  esac
done

###############################
### R:APE:ACE for gene loss ###
###############################
## Assuming genes cannot be gained
## Or that they can...

outdir=${main_dir}/${outdir}
mkdir $outdir
cd $outdir

## Make a list of HOGs with a positive binary in ALL genomes
## Note that there should not be any rows in the binary table where 0 genomes have the HOG IF the whole subtree is included
## But this won't be true for smlaler subtrees (i.e. the T. az. subtree)
## First, the unanimous (all 1's) HOGs
cat ${binary_table} | \
tail -n +2 | \
awk -F "\t" '
{
for (i=2;i<=NF;++i) { if ($i == 0) { next } }
print $1;
}
' \
> unanimous_hogs.list
## Now the absent (all 0's) HOGs
cat ${binary_table} | \
tail -n +2 | \
awk -F "\t" '
{
for (i=2;i<=NF;++i) { if ($i == 1) { next } }
print $1;
}
' \
> absent_hogs.list

## Transpose and trim the binary table to include only HOGs with at least one locus = 1 AND 1+ locus = 0.
## First trim
head -n 1 ${binary_table} > trimmed_table.tmp
cat ${binary_table} | \
awk -F "\t" '
{
n=0;
p=0;
for (i=2;i<=NF;++i) {
	if ($i == 0) { n=1 }
	else if ($i == 1) {p=1}
}
if (p==1 && n==1) {print $0};
}
' \
>> trimmed_table.tmp



## Now transpose
outfile=trimmed_transposed_table.tsv
> $outfile
i=1
while (( $i <= $( wc -l < trimmed_table.tmp ) )); do
head -n $i trimmed_table.tmp | \
tail -n 1 | \
awk ' { for (j=1;j<=NF;++j) {print $j}; } ' \
> this_col.tmp
paste $outfile this_col.tmp > added_col.tmp
mv added_col.tmp $outfile
i=$(( $i + 1 ))
done
sed -i 's/^\t//' $outfile
rm trimmed_table.tmp this_col.tmp

## Copy and prune to the subtree portion of the species tree
## This is done elsewhere in the larger pipeline, but since it is so quick, I might as well keep it more self-contained here
conda deactivate
conda activate gotree
gotree prune -r -i $in_tree -f ${subtree_list} -o pruned_tree.txt
conda deactivate

## Run ape_ace_discrete_wrapper.R on all appropriate HOGs
table_path=${outdir}/trimmed_transposed_table.tsv
pruned_tree=${outdir}/pruned_tree.txt
R_outfile=${outdir}/ape_ace_out.tsv
Rscript --vanilla ~/scripts/helper/ape_ace_discrete_wrapper.R approach=$approach table_path=$table_path tree_path=$pruned_tree outfile=$R_outfile
## Will return warnings, but they seem fine to me.

## R_outfile contains the order of the nodes from each run of APE ACE so I can check it is consistent
cat $R_outfile | awk 'NR%2==1' | sort -u | wc -l
echo "should print 1:"
## 1. All the same! Remove these lines!
head -n 1 $R_outfile > out.tmp
cat $R_outfile | awk 'NR%2==0' >> out.tmp
mv out.tmp $R_outfile


###########################
### Determine gene loss ###
###### probabilities ######
###########################

## The probability that a gene was lost at a particular node should be equal to the probabilities that it was present in each ancestor and the probability that it is not present in the node of interest, all multiplied together.

## Get a table with a list of ancestor nodes for each node
## Ancestor code is borrowed from label_trees_not_ingroup_root.sh so it does a bit more than necessary for just this purpose.


conda deactivate
conda activate gotree
in_tree=$pruned_tree
##Get list of leaves:
gotree labels -i $in_tree > all_leaves.list

#Get list of internal nodes
#If this call is buggy, just use the fact that all the internal nodes are N followed by a number.
cat $in_tree | tr -c "[[:alnum:]]_" "\n" | grep "^[[:alpha:]]" | grep -vFf all_leaves.list > internal_nodes.list

# descendents.tsv will have col1 be an internal node
# and col 2 will be a csv of descendent nodes/leaves
# and the same thing but only for descendent leaves
# Initiate the column that will list all descendents.
> descendents.tmp
#initiate the column that will list only descendents that are leaves.
> leaf_descendents.tmp
cat internal_nodes.list | \
#loop through each leaf/node 
while read node; do
#get the subtree rooted at that node using gotree
gotree subtree -n "^${node}$" -i $in_tree -o sub_tree.tmp
#Get the list of leaves and nodes in the sub tree
#head -n -1 removes self from list (always last....I hope)
cat sub_tree.tmp | tr -c "[[:alnum:]]_" "\n" | grep "^[[:alpha:]]" | head -n -1 > these_descendents.tmp
#get list of all descendents
cat these_descendents.tmp | tr "\n" "," | awk '{print $0}' | sed 's/,$//' >> descendents.tmp
#get list of only leaf descendents
grep -Ff all_leaves.list these_descendents.tmp | tr "\n" "," | awk '{print $0}' | sed 's/,$//' >> leaf_descendents.tmp
#Finish loop through internal nodes
done
#paste together output files and remove .tmp files
paste internal_nodes.list descendents.tmp > descendents.tsv
paste internal_nodes.list leaf_descendents.tmp > leaf_descendents.tsv

##Generate list of ancestor nodes
#Get list of all nodes, internal and leaf
cat all_leaves.list internal_nodes.list > all_nodes.list
##Initiate the column that will list all ancestor nodes
> ancestors.tmp
#Loop through all internal nodes and leaves
cat all_nodes.list | \
while read prefix; do
cat descendents.tsv | awk -v prefix=$prefix -F "[\t,]" 'BEGIN {f=0}; { for (i=2; i<=NF; ++i) { if ($i==prefix) {printf $1","; f=1} } }; END { if (f==0) {print "-"}; print ""}' | sed 's/,$//' >> ancestors.tmp
done
#paste together output files and remove .tmp files
head -n -1 ancestors.tmp > tmp
paste all_nodes.list tmp > ancestors.tsv
rm leaf_descendents.tmp sub_tree.tmp these_descendents.tmp descendents.tmp ancestors.tmp tmp

##Check that ancestors always come after descendents in the newick tree and thus in all my lists
##This will make narrowing common ancestors down to the LCA easy (just the first CA in the list)
> check_order.tmp
nr=1
while (( $nr <= $(wc -l < descendents.tsv) )); do
echo "row number: $nr" >> check_order.tmp
head -n $nr descendents.tsv > this_block.tmp
## Adding the tab prevents partial overlap hits on the grep -Ff these_descendents.tmp
tail -n 1 this_block.tmp | awk -F "[\t,]" '{for (i=2; i<=NF; ++i) {print $i} }' | grep -vFf all_leaves.list | sed 's/$/\t/' > these_descendents.tmp
a=$(cat this_block.tmp | awk -F "\t" '{print $1"\t"}' | grep -Ff these_descendents.tmp | wc -l) 
echo $a >> check_order.tmp
b=$(wc -l < these_descendents.tmp)
echo $b >> check_order.tmp
if (( $a == $b )); then echo "equal" >> check_order.tmp; else echo "NOT EQUAL!" >> check_order.tmp; fi
nr=$(($nr+1))
done
echo "Are descendent nodes always before ancestor nodes?"
if (( $(grep -c "^equal" check_order.tmp) == $( wc -l < descendents.tsv )  )) && (( $(grep -c "^NOT EQUAL!" check_order.tmp) == 0 )); then echo "yes!"; else echo "No! check check_order.tmp"; fi
rm check_order.tmp this_block.tmp these_descendents.tmp

## Remove extra files generated because this code was borrowed from existing code
rm all_leaves.list internal_nodes.list leaf_descendents.tsv descendents.tsv all_nodes.list

## Initiate output table
outtable=ape_ace_gene_loss_gain_probability.tsv
awk '{print $1}' ${R_outfile} > $outtable

## Loop through the nodes/leaves, not including N35 (the root of the T. az. tree)
cat ancestors.tsv | \
head -n -1 | \
awk '{print $1}' | \
while read node; do

## initiate output table column for this node
echo -e $node"\t"$node > this_output_col.tmp
## Initiate a table of the node of interest and the ancestors
> these_ancestors.tsv
## Get the list of ancestors of this node
## First "ancestor" is the node itself. Second is the immediate ancestor
cat ancestors.tsv | \
awk -v node=$node -F "\t|," '
$1==node { for (i=1;i<=NF;++i) { print $i } }
' | \
## Add the ${R_outfile} for each ancestor to these_ancestors.tsv
while read ancestor; do
cat ${R_outfile} | \
awk -v ancestor=$ancestor '
NR==1 { for (i=1;i<=NF;++i) { if ($i==ancestor) {c=i} } };
{ print $c };
' \
> this_ancestor_col.tmp
paste these_ancestors.tsv this_ancestor_col.tmp > tmp
mv tmp these_ancestors.tsv
## Finish loop to create these_ancestors.tsv
done
sed -i 's/^\t//' these_ancestors.tsv

## calculate the probability that each HOG was lost or gained at the node of interest
cat these_ancestors.tsv | \
tail -n +2 | \
awk ' { print (1-$1)*$2"\t"(1-$2)*$1 } ' \
>> this_output_col.tmp

## add this_output_col to outtable!
paste $outtable this_output_col.tmp > tmp
mv tmp $outtable

## Finish loop through nodes of interest
done
rm these_ancestors.tsv this_ancestor_col.tmp this_output_col.tmp

## Split gene loss and gene gain
cat $outtable | awk -F "\t" ' {printf $1; for (i=2;i<=NF;i=i+2) {printf "\t"$i}; print ""} ' > gene_loss_probability.tsv
cat $outtable | awk -F "\t" ' {printf $1; for (i=3;i<=NF;i=i+2) {printf "\t"$i}; print ""} ' > gene_gain_probability.tsv
rm $outtable

## Remove likely transpon-associated genes.
## And genes with nothing all of the following: kegg_ko, description, and preferred name from eggnog-mapper
## And genes with no annotations! (some annotations.tsv files are empty)
> hogs_that_pass.tmp
cat gene_loss_probability.tsv | \
tail -n +2 | \
awk '{print $1}' | \
while read hog; do
annotation_file=~/data/main_pipeline/annotations/eggnog-mapper_out/hog_annotations/${hog}_annotations.tsv
if ! grep -iq "nuclease\|transposase\|integrase\|helicase" $annotation_file; then
cat $annotation_file | awk -v hog=$hog -F "\t" '$8 != "-" || $9 != "-" || $12 != "-" {print hog; exit} ' >> hogs_that_pass.tmp
fi
done
## Manually inspected for first run with only T. az. in the analysis
## what genes got eliminated by each of the two filters above. Looks good.
echo "loss,gain" | \
sed 's/,/\n/' | \
while read change; do
head -n 1 gene_${change}_probability.tsv > only_annotated_gene_${change}_probability.tsv
grep -Ff hogs_that_pass.tmp gene_${change}_probability.tsv >> only_annotated_gene_${change}_probability.tsv
done

## Count the cumulative probabilities of gene loss or gain for each node
echo "loss,gain" | \
sed 's/,/\n/' | \
while read change; do
cat only_annotated_gene_${change}_probability.tsv | \
awk '
NR==1 {
for (i=2;i<=NF;++i) { names[i-1]=$i };
for (i=2;i<=NF;++i) { A[i-1]=0 };
next;
};
{ for (i=2;i<=NF;++i) { A[i-1]+=$i }; }
END { for (i=1;i<=length(A);++i) { print names[i]"\t"A[i]  } };
' \
> cumulative_gene_${change}_probability.tsv
done

## Normalize the cumulative probabilities to the branch lengths
conda deactivate
conda activate gotree
## Get the lengths of each edge named by the node it leads to
gotree stats edges -i $pruned_tree | \
awk '
NR==1 { for (i=1;i<=NF;++i) { if ($i=="rightname") {a=i} else if ($i=="length") {b=i} } ; next};
{ print $a"\t"$b };
' \
> edge_lengths.tsv
## Put nodes in same order for cumulative p of gene loss and edge length
echo "loss,gain" | \
sed 's/,/\n/' | \
while read change; do
sort -k1,1 cumulative_gene_${change}_probability.tsv > p.tmp
sort -k1,1 edge_lengths.tsv > l.tmp
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

## Make metadata doc for iTOL
mkdir iTOL_metadata
echo "loss,gain" | \
sed 's/,/\n/' | \
while read change; do
cat ~/iTol_metadata_template.txt cumulative_gene_${change}_probability.tsv > iTOL_metadata/cum_p_gene_${change}_iTOL_metadata.txt
cat ~/iTol_metadata_template.txt rate_of_cumulative_gene_${change}_probability.tsv > iTOL_metadata/rate_of_cum_p_gene_${change}_iTOL_metadata.txt
cat ~/iTol_metadata_template.txt logged_rate_of_cumulative_gene_${change}_probability.tsv > iTOL_metadata/log_rate_of_cum_p_gene_${change}_iTOL_metadata.txt
done

