## Run simple t-tests on length-bias 2 col data from prokka vs refseq pipeline
## Run paired t-tests on other 2 col data from prokka vs refseq pipeline
## /Users/liamfriar/Desktop/Taz_stuff/analyses is a copy of
## liam@photobiont.colorado.edu:/home/liam/data/main_pipeline/RefSeq_vs_prokka/analyses

## Run simple t-tests on length-bias 2 col data from prokka vs refseq pipeline
cd /Users/liamfriar/Desktop/Taz_stuff/analyses
outfile=length_ttest_outputs.txt
> $outfile

## Loop through prefixes
cat N13.shortened.tsv | \
while read long_prefix short_prefix; do

## Loop through the 6 input files from each prefix
echo -e "\
prokka_countDiffs_length_ttest.tsv
prokka_presAbs_length_ttest.tsv
prokka_poorOverlap_length_ttest.tsv
refseq_countDiffs_length_ttest.tsv
refseq_presAbs_length_ttest.tsv
refseq_poorOverlap_length_ttest.tsv\
" | \
while read filename; do

infile=${short_prefix}/${filename}
#--vanilla tells R not to look at a user's history, personal settings, etc.
Rscript --vanilla /Users/liamfriar/Documents/Research/Azolla/Genomics/Scripts/Local/Ttests.R infile=$infile outfile=$outfile paired=FALSE header=FALSE

## Finish loop through input filenames
done

## Finish loop through prefixes
done

## Run paired t-tests on other 2 col data from prokka vs refseq pipeline
cd /Users/liamfriar/Desktop/Taz_stuff/analyses
outfile=paired_ttest_outputs.txt
> $outfile

## Loop through 4 input files:
echo -e "\
n_hogs_ttest.tsv
n_intact_ttest.tsv
n_orfs_ttest.tsv
p_intact_ttest.tsv\
" | \
while read filename; do

infile=${filename}
#--vanilla tells R not to look at a user's history, personal settings, etc.
Rscript --vanilla /Users/liamfriar/Documents/Research/Azolla/Genomics/Scripts/Local/Ttests.R infile=$infile outfile=$outfile paired=TRUE header=TRUE

## Finish loop through input filenames
done

## use a csv in the middle instead of tsv because of my local bash
## Arrange the length ttest outputs into a tsv
cat length_ttest_outputs.txt | \
sed 's/\//,/' | \
sed 's/_length_ttest.tsv$//' | \
awk ' NR%3==0 { print $0; next }; { printf $0","} ' \
> length_ttest_outputs.csv

## Split into one table for prokka, one for refseq
for i in prokka refseq; do
grep "$i" length_ttest_outputs.csv | sed "s/${i}_//" > ${i}.tmp
## Rearrange
echo -e "genome\tdiff_counts_ttest\tdiff_counts_pval\tpres_abs_ttest\tpres_abs_pval\tpoor_overlap_ttest\tpoor_overlap_pval" > ${i}_length_ttest_outputs.tsv
cat ${i}.tmp | awk -F "," ' BEGIN { prev="" }; $1 != prev { prev=$1; printf "\n"$1"\t" }; {printf $3"\t"$4"\t"}; END { print "" } ' | \
tail -n +2 | \
sed 's/[[:blank:]]*$//' \
>> ${i}_length_ttest_outputs.tsv
done

rm *tmp
rm length_ttest_outputs.csv length_ttest_outputs.txt

## Rearrange the paired ttest outputs into a table
echo "## negative ttest means refseq values lower, positive means refseq values higher (than prokka values) ##" > tmp
echo -e "dataset\tttest\tpval" >> tmp
cat paired_ttest_outputs.txt | \
sed 's/_ttest.tsv//' | \
awk ' NR%3==0 { print $0; next }; { printf $0"\t"} ' \
>> tmp
mv tmp paired_ttest_outputs.tsv
rm paired_ttest_outputs.txt
 
#scp paired_ttest_outputs.tsv prokka_length_ttest_outputs.tsv refseq_length_ttest_outputs.tsv liam@photobiont.colorado.edu:/home/liam/data/main_pipeline/RefSeq_vs_prokka/analyses
#cd ..
#rm -r analyses
