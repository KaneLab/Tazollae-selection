## Label trees.
## hyphy label-tree.bf has a --internal-nodes All option
## That automatically labels all internal nodes that have all descendants labeled, but I don't want the root of the ingroup to be labeled.


## Get variables
print_usage() {
  printf "usage requires:\n\
and in_tree -i\n\
and logfile -l\n\
and ingroup prefix list -I\n\
and outgroup prefix list -O\n\
and out_tree -o\
"
}

while getopts 'i:l:I:O:o:' flag; do
  case "${flag}" in
    i) in_tree="${OPTARG}" ;;
    l) logfile="${OPTARG}" ;;
    I) ingroup_list="${OPTARG}" ;;
    O) outgroup_list="${OPTARG}" ;;
    o) out_tree="${OPTARG}" ;;
    *) print_usage
       exit 1 ;;
  esac
done


#############################
#### Prepare files to be ####
##### used by algorithms ####
#############################

conda deactivate
conda activate gotree
#in_tree=single_copy_HOGs.pruned_tree.txt
#in_msa=single_copy_HOGs.fa

##Get list of leaves:
gotree labels -i $in_tree > all_leaves.list
#Get list of internal nodes
#If this call is buggy, just use the fact that all the internal nodes are N followed by a number.
cat $in_tree | tr -c "[[:alnum:]]_" "\n" | grep "^[[:alpha:]]" | grep -vFf all_leaves.list > internal_nodes.list

# descendents.tsv will have col1 be an internal node
# and col 2 will be a csv of descendent nodes/leaves
# and the same thing but only for descendent leaves
#Initiate the column that will list all descendents.
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
echo "Are descendent nodes always before ancestor nodes?" >> $logfile
if (( $(grep -c "^equal" check_order.tmp) == $( wc -l < descendents.tsv )  )) && (( $(grep -c "^NOT EQUAL!" check_order.tmp) == 0 )); then echo "yes!" >> $logfile; else echo "No! check check_order.tmp" >> $logfile; fi
rm check_order.tmp this_block.tmp these_descendents.tmp

##Add the leaves to the descendents list
paste all_leaves.list all_leaves.list > descendents_incl_leaves.tsv
cat descendents_incl_leaves.tsv descendents.tsv > tmp
mv tmp descendents_incl_leaves.tsv

##############################
#### The actual algorithm ####
##############################

## Get the LCA of the test group

## Get the LCA of a set of leaves
## input is a list of nodes, with one node per line.
## using the order of the nodes (descendents before ancestors) makes this computationally easier, but it might allow some bugs to go unnoticed.

##grep each node in the input list sequentially from descendents.tsv
##which should leave only lines containing all of the nodes.
#These are all the common ancestors.
grep -Ff $ingroup_list all_leaves.list > test_leaves.tmp

cp descendents_incl_leaves.tsv tmp
cat test_leaves.tmp | \
while read node; do
cat tmp | awk -v node=$node -F "[\t,]" '{ for (i=2; i<=NF; ++i) { if ($i==node) {print $0","; } } }' | sed 's/,$//' > tmp.tmp
mv tmp.tmp tmp
done
# Because descendents are always before ancestors, head -n 1 will get the LCA.
# as matches further down the list will be earlier CAs
LCA=$( head -n 1 tmp | awk -F "\t" '{print $1}' )
rm tmp test_leaves.tmp

## Now get all descendants from LCA
cat descendents_incl_leaves.tsv | \
awk -v lca=$LCA -F "\t" '$1==lca {print $1"\n"$2; exit}' | \
tr "," "\n" \
> test.list
## Get all nodes other than LCA and the test nodes
grep -vFf test.list all_nodes.list > reference.list
## Remove LCA from the test.list
tail -n +2 test.list > tmp
mv tmp test.list

## Also get the root of the entire tree and remove this from reference.list
## Again using the fact that descendent nodes come first, so the root node should be last
root_node=$(cat descendents.tsv | awk '{print $1}' | tail -n 1)

## make mapfile for gotree rename
## 2col tsv of old name and new name
cat test.list | awk '{print $1"\t"$1"{test}"}' > mapfile.tmp
cat reference.list | awk '{print $1"\t"$1"{reference}"}' >> mapfile.tmp
echo -e "${LCA}\t${LCA}" >> mapfile.tmp
echo -e "${root_node}\t${root_node}" >> mapfile.tmp

## Label
conda deactivate
conda activate gotree
gotree rename --internal -m mapfile.tmp -i $in_tree -o $out_tree

## Clean up
rm internal_nodes.list leaf_descendents.tsv descendents.tsv all_nodes.list ancestors.tsv \
descendents_incl_leaves.tsv all_leaves.list test.list reference.list mapfile.tmp
