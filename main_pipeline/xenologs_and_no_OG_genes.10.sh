##Look at the annotations for the loci that were not assigned HOGs by orthofinder
##Are any of them annotated meaningfully?
cd ~/data/main_pipeline/orthofinder_out/for_downstream
#loop through Nostocales
outfile=ungrouped_gene_annotations.tsv
>$outfile
> unassigned_loci.list
cat ~/subtree.shortened.tsv | \
while read long_prefix short_prefix; do
echo $short_prefix
##Get the list of loci that were assigned to a N9.HOG
cat ~/data/main_pipeline/orthofinder_out/for_downstream/N9_noOutgroups.tsv | \
awk -v prefix=$short_prefix -F "\t" '
NR == 1 { for (i=1; i<=NF; ++i) {if ($i==prefix) {c=i} }; next};
{print $c}' | \
tr "," "\n" | \
tr -d " " | \
sed '/^$/d' \
> assigned_loci.tmp

##Get the list of loci that were not assigned to an N9.HOG
##Get the full list of loci
orthofinder_indir=~/data/main_pipeline/orthofinder_out/input_fastas
cat ${orthofinder_indir}/${short_prefix}.fasta | grep "^>" | sed 's/[[:space:]].*//' | sed 's/^>//' > all_loci.tmp
##Find the difference between the two lists
grep -vFf assigned_loci.tmp all_loci.tmp > unassigned_loci.tmp

##Get the annotation lines
annotation_file="~/data/main_pipeline/annotations/eggnog-mapper_out/genome_annotations/${short_prefix}.tsv"
grep -Ff unassigned_loci.tmp $annotation_file >> $outfile

## Save the list of unassigned loci
cat unassigned_loci.tmp >> unassigned_loci.list

done
rm assigned_loci.tmp all_loci.tmp unassigned_loci.tmp


## How many unassigned loci from each genome?
outfile=unassigned_loci_counts.tsv
cat unassigned_loci.list | sed 's/_.*//' | sort | uniq -c | awk '{print $2"\t"$1}' > $outfile
## Out of how many intact loci?
> intact_count.tmp
cat $outfile | \
awk '{print $1}' | \
while read prefix; do
#long_prefix=$( grep ${prefix} ~/subtree.shortened.tsv | awk '{print $1}' )
grep -c "^>" ${orthofinder_indir}/${prefix}.fasta >> intact_count.tmp
done
paste $outfile intact_count.tmp > paste.tmp
mv paste.tmp $outfile
rm intact_count.tmp
## Representing what portion?
cat $outfile | \
awk '
BEGIN {print "genome\tn_unassigned\tn_intact\tportion"};
{print $0"\t"$2/$3}
' \
> awk.tmp
head -n 1 awk.tmp > $outfile
tail -n +2 awk.tmp | sort -k4,4g >> $outfile
rm awk.tmp
