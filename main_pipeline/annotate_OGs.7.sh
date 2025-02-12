## Organize the eggnog-mapper annotations to be more easily useful down the pipeline

## I concatenated all of the intact genes from the orthofinder subtree
## and submitted these to the eggnog-mapper online database with default parameters
## Versions and the exact call to eggnog-mapper can be found in that output.

## These variables were in a loop when I had also used alignments to RefSeq genomes to annotate
## I am leaving them with possibly unnecessary variable naming to prevent bugs due to removing the loop
annotation_file=~/data/main_pipeline/annotations/eggnog-mapper_out/online_intact_annotations/out.emapper.annotations
tool=eggnog-mapper_out

## Split the eggnog-mapper outfile by prefix to speed up search below.
cd ~/data/main_pipeline/annotations/${tool}
mkdir genome_annotations
cat ~/subtree.prefix.list | \
while read prefix; do
grep "^${prefix}" ${annotation_file} > genome_annotations/${prefix}.tsv
done
## Total lines in genome_annotations files is same or 1 less than in annotation_file 
## depending on if annotation_file has a header line

## For each hog, create a file containing all eggNOG-mapper annotations for loci in that hog
mkdir hog_annotations
locus_list_dir=~/data/main_pipeline/orthofinder_out/for_downstream/locus_lists_HOGs
## Loop through hogs
ls ${locus_list_dir}/*_loci.tsv | \
while read locus_list; do
hog=${locus_list##*/}
hog=${hog%_loci.tsv}
echo $hog
> hog_annotations/${hog}_annotations.tsv
## Loop through the loci in each hog
cat $locus_list | \
awk -F "\t|, " ' $2 != "" { for (i=2;i<=NF;++i) {print $1"\t"$i} } ' | \
while read prefix locus; do
cat genome_annotations/${prefix}.tsv | \
awk -v locus=$locus '$1 == locus' \
>> hog_annotations/${hog}_annotations.tsv
## Finish loop through loci in this hog
done
## Finish loop through HOGs
done

###########################
### Get a file for each ###
### type of annotation  ###
###########################

annotation_file=~/data/main_pipeline/annotations/eggnog-mapper_out/online_intact_annotations/out.emapper.annotations
tool=eggnog-mapper_out

cd ~/data/main_pipeline/annotations/${tool}/hog_annotations
## Loop through annotation columns of interest
echo "KO,12
gene_name,9
KEGG_module,14
KEGG_pathway,13
COG_category,7" | \
sed 's/,/\t/g' | \
while read annotation target_col; do
echo $annotation

## Create a tsv with first col the hog, second col a csv of annotations
outfile=search_${annotation}_by_hog.tsv
> $outfile
ls N9*_annotations.tsv | \
while read file; do
hog=${file%_annotations.tsv}
cat $file | \
awk -v target_col=$target_col -F "\t" ' $target_col != "-" {print $target_col} ' | \
tr "," "\n" \
> if.tmp
## the COG_categories are not comma-separated
if [[ $annotation == "COG_category" ]]; then
cat if.tmp | awk -F "" '{for (i=1;i<=NF;++i) {print $i}}' > if.tmp.tmp
mv if.tmp.tmp if.tmp
fi
## The KEGG_pathways contain both ko and map which I have confirmed are redundant
if [[ $annotation == "KEGG_pathway" ]]; then
cat if.tmp | grep "^ko" > if.tmp.tmp
mv if.tmp.tmp if.tmp
fi
cat if.tmp | \
sort -u | \
awk -v hog=$hog '
BEGIN {printf hog"\t"};
{ if (FNR==1) {printf $0} else {printf ","$0} };
END { if (FNR==0) {print "-"} else {print ""} };
' \
>> $outfile
rm if.tmp
## Finish loop through hog annotation files
done

## To improve future searches, also make a table with each row the annotation
## followed by hogs with that annotation
infile=$outfile
outfile=search_hog_by_${annotation}.tsv
> $outfile
cat $infile | \
awk -F "\t" ' { split($2,A,","); for (i=1;i<=length(A);++i) {print A[i]} } ' | \
sort -u | \
awk '$0!="-"' | \
## And find the hogs with that annotation term
while read annotation_term; do
cat $infile | \
awk -v annotation_term=$annotation_term -F "\t" '
BEGIN {printf annotation_term};
{ split($2,A,","); for (i=1;i<=length(A);++i) { if (A[i]==annotation_term) { printf ","$1; next } } };
END {print ""}
' | \
sed 's/,/\t/' \
>> $outfile
## finish annotation_term loop
done

## End annotation target_col loop
done

## The description column is prose and contains punctuation
## So parsing in the above loops is difficult
## So just take the annotation from the closest genome to T. az.
## Get the ordered list of closeness to T. az. (using NosAzo0708 as the center point)
conda deactivate
conda activate gotree
in_tree=~/data/main_pipeline/orthofinder_out/final_run/Species_Tree/SpeciesTree_rooted_node_labels.txt
gotree matrix -i $in_tree -o distance_matrix.tmp
## Get the list of nodes in the order that gotree matrix output them
cat distance_matrix.tmp | awk -F "\t" '{print $1}' | tail -n +2 > prefix_col.tmp
cat distance_matrix.tmp | awk -F "\t" ' $1=="NosAzo0708" {for (i=2;i<=NF;++i) {print $i} } ' > distance_col.tmp
paste prefix_col.tmp distance_col.tmp | \
grep -Ff ~/subtree.prefix.list \
> distance_from_NosAzo0708.tsv
rm prefix_col.tmp distance_col.tmp distance_matrix.tmp
## Sort so shortest distance on top, longest on bottom
cat distance_from_NosAzo0708.tsv | sort -k2,2g > sorted_distance_from_NosAzo0708.tsv
cat sorted_distance_from_NosAzo0708.tsv | awk '{print $1}' > ascending_distance_from_NosAzo0708.tmp
rm distance_from_NosAzo0708.tsv
conda deactivate
## Get one description per hog
annotation=description
target_col=8
outfile=search_${annotation}_by_hog.tsv
> $outfile
## Loop through hogs
ls N9*_annotations.tsv | \
while read file; do
hog=${file%_annotations.tsv}
echo $hog

## if no annotations, assign description of "-" and skip all the code below
if (( $( wc -l < $file ) == 0 )); then
echo -e "${hog}\t-" >> $outfile;
continue;
fi

> this_hog_description.tmp
## i will increment through the genomes to take the annotation from
i=1
i_max=$( wc -l < ascending_distance_from_NosAzo0708.tmp )
while (( $( wc -l < this_hog_description.tmp  ) == 0 )); do
## if have already looped through all prefixes, then assign annotation of "-"
if (( $i > $i_max )); then
	echo -e "${hog}\t-" > this_hog_description.tmp;
else
	prefix=$( cat ascending_distance_from_NosAzo0708.tmp | head -n $i | tail -n 1 )
	cat $file | \
	awk -v target_col=$target_col -v prefix=$prefix -v hog=$hog -F "\t" '
	$1 ~ prefix && $target_col != "-" {print hog"\t"prefix"\t"$target_col; exit}
	' \
	> this_hog_description.tmp
fi
## increment!
i=$(( $i+1 ))
## finish wc -l this_hog_description.tmp while loop
done

cat this_hog_description.tmp >> $outfile
rm this_hog_description.tmp
## Finish loop through hogs
done
rm ascending_distance_from_NosAzo0708.tmp

## Combine the annotations into a single table
## I confirmed that each search_${annotation}_by_hog.tsv has the same hogs in the same order
echo "HOG" > hog_annotations.tsv
cat search_KO_by_hog.tsv | awk -F "\t" '{print $1}' >> hog_annotations.tsv
echo "gene_name
KO
KEGG_module
KEGG_pathway
COG_category
description" | \
while read annotation; do
echo $annotation > this_col.tmp
cat search_${annotation}_by_hog.tsv | awk -F "\t" '{print $NF}' >> this_col.tmp
paste hog_annotations.tsv this_col.tmp > paste.tmp
mv paste.tmp hog_annotations.tsv
rm this_col.tmp
done
#rm search_*_by_hog.tsv



##########################
### Uniquely annotated ###
##########################

## Which HOGs are the only ones to contain a given ko or preffered name?
## Include only HOGs that have at least one T. az. locus.
## So, if ko1234 is in HOG01 and HOG02, but HOG02 contains only free-living,
## then HOG01 counts as the unique HOG for ko1234.
## THis is because I am curious which HOGs are unique genes for the T. az.
## Note the T. az. in that HOG does not have to have that ko/name annotation, but just be present in the HOG
## Get a list of T. az.-containing HOGs.
> HOGs_that_contain_Taz.list
ls ${locus_list_dir}/N9.HOG*_loci.tsv | \
while read file; do
if cat $file | awk -F "\t" '$2!="" {print $2}' | grep -qFf ~/Taz.short_prefix.list; then
hog=${file##*/}
hog=${hog%_loci.tsv}
echo $hog >> HOGs_that_contain_Taz.list
fi
done

## add "ko:" to list so grep -Ff call below keeps ko: lines and removes only hog lines not in HOGs_that_contain_Taz.list
cp HOGs_that_contain_Taz.list grep_list.tmp
echo "ko:" >> grep_list.tmp
## Split search_hog_by_KO.csv into individual lines so can use grep -Ff to remove hogs that don't contain T. az. loci
cat search_hog_by_KO.tsv | \
tr "," "\n" | \
tr "\t" "\n" | \
## keep only ko: and hogs that contain Taz
grep -Ff grep_list.tmp | \
## Put back into 1 line per ko:
awk '$0~"^ko:" {print ""; printf $0; next}; {printf ","$0}' | \
tail -n +2 | \
## Get a list of ko: with exactly 1 HOG
awk -F "," 'NF==2 {print $2"\t"$1}' | \
## Get one line per hog
sort -k1,1 | \
awk '
BEGIN {prev=""};
$1!=prev {print ""; printf $0; prev=$1; next}
{printf ","$2};
END {print ""}
' | \
tail -n +2 \
> HOGs_w_unique_ko.tsv
rm grep_list.tmp
cat HOGs_w_unique_ko.tsv | awk '{print $1}'  > HOGs_w_unique_ko.list

## Now do the same thing, but for preferred_name instead of ko
## First, get the annotations of the HOGs that contain T. az.
> HOGs_w_unique_names.tsv
> HOGs_with_not_unique_names.list
cat search_gene_name_by_hog.tsv | \
grep -Ff HOGs_that_contain_Taz.list | \
## Split each annotation into its own line
awk -F "\t" ' { split($2,A,","); for (i=1;i<=length(A);++i) {print A[i]"\t"$1} } ' | \
## A small number of gene names have punctuation. Deal with that here:
## If running this on different data, should grep "[[:punct:]]" to make sure these make sense.
## some gene names contain "-" and "_". I checked that there are none where the same gene name appears but without punctuation.
awk -F "\t" ' { split($1,A,"/"); for (i=1;i<=length(A);++i) {print A[i]"\t"$2} } ' | \
## sort by annotation
sort -k1,1 | \
## Remove hogs with identical annotation to previous 
## Note that ccmK1 and ccmK2 would count as different, though I don't know if they actually are
## Keep track of HOGs_with_not_unique_names because in theory a HOG could have 2 names,
## one of which is unique and would therefore pass, and one of which is not unique and therefore not pass
## I want that HOG to then not pass.
## This eliminates 3 HOGs!
awk '
BEGIN {prev_name=""; prev_hog=""; p=0};
{
if ($1!=prev_name) {
	if (p==1) {
		print prev_hog"\t"prev_name >> "HOGs_w_unique_names.tsv"
	} else {
		if (NR!=1) { print prev_hog >> "HOGs_with_not_unique_names.list" }
	};
	p=1;
} else {
	print prev_hog >> "HOGs_with_not_unique_names.list"
	p=0
};
prev_name=$1;
prev_hog=$2;
}
'
grep -vFf HOGs_with_not_unique_names.list HOGs_w_unique_names.tsv > grep.tmp
mv grep.tmp HOGs_w_unique_names.tsv
## Get one line per hog
cat HOGs_w_unique_names.tsv | \
sort -k1,1 | \
awk '
BEGIN {prev=""};
$1!=prev {print ""; printf $0; prev=$1; next}
{printf ","$2};
END {print ""}
' | \
tail -n +2 \
> one_line.tmp
mv one_line.tmp HOGs_w_unique_names.tsv
cat HOGs_w_unique_names.tsv | awk '{print $1}' > HOGs_w_unique_names.list

## Finally, get the overlap between HOGs with unique preferred names and unique ko:s
## If the HOG has two lines, that means it was in both lists.
## So combine the annotations and print:
cat HOGs_w_unique_names.tsv HOGs_w_unique_ko.tsv | \
sort -k1,1 | \
awk '
BEGIN {prev_hog=""; prev_ann=""};
{
if ($1==prev_hog) {
	print prev_hog"\t"prev_ann"\t"$2
} else ;
prev_hog=$1
prev_ann=$2;
}
' | \
## put ko and name in consistent order
awk ' $2~"^ko:" {print $0; next}; {print $1"\t"$3"\t"$2} ' \
> HOGs_w_unique_names_AND_KOs.tsv
## Get just the list of HOGs
cat HOGs_w_unique_names_AND_KOs.tsv | awk '{print $1}' > HOGs_w_unique_names_AND_KOs.list


## The KEGG pathways and modules seemed potentially out of sync
## with the current KEGG database, so I got the constituent lists
## and human language definitions from an R package (KEGGREST) on my local computer
## THESE ARE WHAT I WILL USE IN DOWNSTREAM ANALYSES!
## So, keep the eggnog-assigned COGs, gene names, and ko's, but use these pathways and moudules.
cd ~/data/main_pipeline/annotations/eggnog-mapper_out/hog_annotations
KEGG_modules_table=~/KEGG_module_kos.tsv
KEGG_pathways_table=~/KEGG_pathway_kos.tsv

echo "$KEGG_modules_table
$KEGG_pathways_table" | \
while read kegg_table; do

## module_list.tmp while list all modules associated with all KOs for each HOG
> module_list.tmp
## Loop through HOGs to get all associated KO's for each HOG
cat search_KO_by_hog.tsv | \
awk -F "\t" '{print $2}' | \
sed 's/ko://g' | \
while read ko_list; do

## if no associated modules, add "-" to module_list.tmp
if [[ $ko_list == "-" ]];
then
echo "-" >> module_list.tmp
else

## Get modules associated with those KO's
echo $ko_list | tr "," "\n" > ko_list.tmp
grep -Ff ko_list.tmp $kegg_table > these_modules.tmp

## Append list of modules
## Or "-" if none
if (( $( wc -l < these_modules.tmp ) == 0 )); then
echo "-" >> module_list.tmp
else
cat these_modules.tmp | \
awk -F "\t" '{print $1}' | \
tr "\n" "," | \
sed 's/,$//' | \
awk '{print $0}' \
>> module_list.tmp
## Finish if $(( $( wc -l < these_modules.tmp ) == 0 ));
fi
rm ko_list.tmp these_modules.tmp


## Finish if [[ $ko_list == "-" ]];
fi
## Finish loop while read ko_list
done

## Finish loop while read kegg_table
done

## Need to change the table output names
## As it is, I just switched out the corresponding columns in hog_annotations.tsv
