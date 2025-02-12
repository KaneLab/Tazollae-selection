## Generate binary intact/not and presence/absence tables
## by gene_name or KO instead of by HOG.

annotation_dir=~/data/main_pipeline/annotations/eggnog-mapper_out/hog_annotations
core_shell_dir=~/data/main_pipeline/core_vs_shell_analyses
binary_dir=~/data/main_pipeline/gene_loss

main_dir=~/data/main_pipeline/inventory_by_annotation
mkdir $main_dir
cd $main_dir

## Loop KO's and gene names ( two ways to identify genes )
echo "KO
gene_name" | \
while read ontology; do
outdir=${main_dir}/${ontology}
mkdir $outdir
cd $outdir
## Loop intact/not or presence/absence
echo "intact_vs_not
presence_absence" | \
while read binary; do
## Loop through Taz-only and free-living-only binary intables
## Easier than starting with the "all" table and re-splitting
echo "Taz
free_living" | \
while read group; do


## Make a binary table where rows are ${ontology}'s instead of HOGs.
annotation_file=${annotation_dir}/search_hog_by_${ontology}.tsv
binary_intable=${binary_dir}/${group}.${binary}.tsv
binary_outtable=${outdir}/${binary_intable##*/}

## initiate outtable with header from intable
cat $binary_intable | \
head -n 1 | \
awk -F "\t" ' { printf $2; for (j=3;j<=NF;++j) {printf "\t"$j}; print "" } ' \
> $binary_outtable

## Loop through ${ontology}'s
## Use a counter loop instead of just reading the first col of the file
## So can take the whole line of the annotation_file with its inconsistent csv 2nd col
i=1
while (( $i <= $( wc -l < $annotation_file ) )); do
echo $i
## Grab all hogs for the ${ontology} on that line of the annotation file
cat $annotation_file | \
head -n $i | \
tail -n 1 | \
awk -F "\t|," '{ for (j=2;j<=NF;++j) {print $j} }' \
> hog_list.tmp
## Get those lines of the binary intable
grep -Ff hog_list.tmp $binary_intable | \
## And collapse into a single binary row 1 if any rows=1
## start awk for loop at 2 to skip the hog column
awk '
{ for (j=2;j<=NF;++j) { A[j-1]+=$j } };
END { for (j=1;j<=length(A);++j) { if (A[j]>0) {print 1} else if (A[j]==0) {print 0} else {print "ERROR"} } }
' | \
tr '\n' '\t' | \
sed 's/\t$//' | \
awk '{print $0}' \
>> $binary_outtable


## increment counter
i=$(( $i + 1 ))
## finish i while loop
done
rm hog_list.tmp

## Add first col
cat $annotation_file | \
awk -v ontology=$ontology -F "\t" 'BEGIN {print ontology}; {print $1}' \
> first_col.tmp
paste first_col.tmp $binary_outtable > paste.tmp
mv paste.tmp $binary_outtable
rm first_col.tmp

## finish group loop
done

## Finish binary loop
done
## Finish ontology loop
done
