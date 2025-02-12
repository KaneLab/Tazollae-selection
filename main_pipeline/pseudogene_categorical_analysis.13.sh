## Identify HOGs that are present, but not intact

outdir=~/data/main_pipeline/pseudogene_analyses
mkdir $outdir
cd $outdir
indir=~/data/main_pipeline/gene_loss

## Identify HOGs that are present in at least 6 genomes.
## This is parallel to the hyphy RELAX threshold which requires a HOG to be intact in at least 3 T. az. and at least 3 free-living
cat ${indir}/all.presence_absence.tsv | \
tail -n +2 | \
awk '{ c=0; for (i=2;i<=NF;++i) { c+=$i }; if (c>=6) {print $1} }' \
> hogs_present_in_6_genomes.list


## Paste binary tables together
## And use the difference in values to determine hogs in genomes that are present but not intact
## This will therefore ignore multi-copy HOGs where some but not all are pseudogenized
## There are not very many multi-copy HOGs, and the above situation would be less straightforward to interpret

## Loop Taz and free-living binary tables
## Easier than using the all tables and re-splitting
echo "Taz
free_living" | \
while read group; do
mkdir $outdir/$group
cd $outdir/$group

paste ${indir}/${group}.intact_vs_not.tsv ${indir}/${group}.presence_absence.tsv | \
awk -F "\t" '
NR==1 {
	printf $1;
	for (i=2;i<=NF/2;++i) {printf "\t"$i};
	print "";
	next;
};
{
	printf $1;
	for (i=2;i<=NF/2;++i) { printf "\t"; if ($i==0 && $(NF/2+i)==1) {printf 1} else {printf 0} }; 	print ""
}
' \
> pseudo_only.tsv


## Get a list of hogs that are fully pseudogenized in at least one genome
cat pseudo_only.tsv | \
tail -n +2 | \
awk -F "\t" '
{ for (i=2;i<=NF;++i) { if ($i==1) {print $1; next} } }
' \
> pseudo_hogs.list

## Get a list of hogs that are present in at least one genome
## And pass the 6+ genome threshold
cat ${indir}/${group}.presence_absence.tsv | \
tail -n +2 | \
awk -F "\t" '
{ for (i=2;i<=NF;++i) { if ($i==1) {print $1; next} } }
' | \
grep -Ff ../hogs_present_in_6_genomes.list \
> present_hogs.list


## How many genomes are fully pseudogenized for each hog?
## Are most hogs pseudogenized in only 1 genome? in all? Somehwere in between?
cat pseudo_only.tsv | \
tail -n +2 | \
grep -Ff present_hogs.list | \
awk '{c=0; for (i=2;i<=NF;++i) {c+=$i}; print c}' | \
sort -n | \
uniq -c


## Below is adapted from hyphy_relax_categorical_analysis.sh
## binomial test of enrichment of pseudogenization in HOGs of various functional categories

hog_annotation_file=~/data/main_pipeline/annotations/eggnog-mapper_out/hog_annotations/hog_annotations.tsv
## loop COG_category, KEGG_pathway, and KEGG_module
echo "COG_category
KEGG_pathway
KEGG_module" | \
sed 's/,/\t/g' | \
while read ontology; do

## Do the analysis including HOGs without annotations for this ontology and not including those HOGs
echo "inclusive
exclusive" | \
while read inclusive_exclusive; do

## Get the ontology categories for hogs present in this group
> ${ontology}_HOG_annotations.tsv
cat present_hogs.list | \
while read hog; do
cat ${hog_annotation_file} | \
awk -v ontology=$ontology -v hog=$hog -F "\t" '
NR==1 { for (i=1;i<=NF;++i) { if ($i==ontology) {c=i; next} } };
$1==hog {print hog"\t"$c};
' \
>> ${ontology}_HOG_annotations.tsv
done

## Remove hogs with no annotation for this ontology
if [[ $inclusive_exclusive == "exclusive" ]]; then
cat ${ontology}_HOG_annotations.tsv | \
awk -F "\t" '$2!="-"' \
> awk.tmp
mv awk.tmp ${ontology}_HOG_annotations.tsv
fi

## Get the portion of all results that are positive or negative
## Can be used by both COG and KEGG_pathway analyses
#i# This could be done outside of the ontology loop were it not for the inclusive_exclusive loop
global_n=$( wc -l < ${ontology}_HOG_annotations.tsv )
global_number_positive=$( grep -Ff pseudo_hogs.list ${ontology}_HOG_annotations.tsv | wc -l )
global_percent_positive=$( awk -v numer=$global_number_positive -v denom=$global_n 'BEGIN {print numer/denom}' )

## For each category present, get the number of present HOGs, number of pseudogenized HOGs, more/less than global, p-val
## Initiate output columms
echo "${ontology}" > category_col.tmp
echo "present_HOGs (${global_n})" > n_present_col.tmp
echo "pseudo_HOGs (${global_number_positive})" > n_pseudo_col.tmp
echo "relative_to_global_proportion" > more_less_col.tmp
echo "p" > p_val_col.tmp
echo "${ontology}_def" > definitions_col.tmp
## Loop through categories:
cat ${ontology}_HOG_annotations.tsv | \
awk -F "\t" '{print $2}' | \
tr "," "\n" | \
sort -u | \
while read category; do
echo $category >> category_col.tmp
## Get all present hogs with this category as an annotation
cat ${ontology}_HOG_annotations.tsv | \
awk -v category=$category -F "\t|," ' { for (i=2;i<=NF;++i) { if ($i==category) {print $1} } } ' \
> this_category.tsv
n_hogs_present=$( wc -l < this_category.tsv )
echo $n_hogs_present >> n_present_col.tmp
## all hogs in this category with at least one genome fully pseudogenized
n_hogs_pseudo=$( cat this_category.tsv | grep -Ff pseudo_hogs.list | wc -l )
echo $n_hogs_pseudo >> n_pseudo_col.tmp
## Is this more or less than the global proportion of hogs that are fully pseudogenized in at least one genome?
more_less=$(\
awk -v global_p=$global_percent_positive -v this_pseudo=$n_hogs_pseudo -v this_present=$n_hogs_present '
BEGIN { a=this_pseudo/this_present; if (a>global_p) {print "more"} else if (a<global_p) {print "less"} else {print "equal"} }
'\
)
echo $more_less >> more_less_col.tmp
## p-value!
Rscript --vanilla ~/scripts/helper/binom_test_p_val.R x=$n_hogs_pseudo n=$n_hogs_present p=${global_percent_positive} outfile=p_val_col.tmp
## category definition!
cat ~/${ontology}_definitions.tsv | \
awk -v category=$category -F "\t" 'BEGIN {s=0}; $1==category {print $2; s=1}; END { if (s==0) {print "-"} } ' \
>> definitions_col.tmp


## Finish category loop;
done

## Paste output cols together!
paste category_col.tmp n_present_col.tmp n_pseudo_col.tmp more_less_col.tmp p_val_col.tmp definitions_col.tmp > ${ontology}_pseudogene_enrichment.${inclusive_exclusive}.tsv
rm category_col.tmp n_present_col.tmp n_pseudo_col.tmp more_less_col.tmp p_val_col.tmp definitions_col.tmp this_category.tsv ${ontology}_HOG_annotations.tsv
## sort by significance
head -n 1 ${ontology}_pseudogene_enrichment.${inclusive_exclusive}.tsv > sorted.tmp
cat ${ontology}_pseudogene_enrichment.${inclusive_exclusive}.tsv | \
tail -n +2 | \
sort -k5,5g \
>> sorted.tmp
mv sorted.tmp ${ontology}_pseudogene_enrichment.${inclusive_exclusive}.tsv

## Finish inclusive_exclusive loop
done

## Finish ontology loop
done


## End group loop
done
