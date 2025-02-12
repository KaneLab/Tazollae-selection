## Assign pseudogenes to HOGs via diamond blast search.

## Inputs from main pipeline:
#outdir=~/data/main_pipeline/assign_pseudogenes
#N=N9
#prefix_tsv=~/N9.shortened.tsv
#pseudofinder_outdir=~/data/main_pipeline/pseudofinder_out
#orthofinder_downstream_dir=~/data/main_pipeline/orthofinder_out/for_downstream

## Get variables
print_usage() {
  printf "usage requires:\n\
and outdir -o\n\
and orthofinder subtree root (i.e. N9) -n\n\
and subtree prefix list -l\n\
and pseudofinder_outdir -p\n\
and orthofinder_downstream_dir -f\
"
}

while getopts 'o:n:l:p:f:' flag; do
  case "${flag}" in
    o) outdir="${OPTARG}" ;;
    n) N="${OPTARG}" ;;
    l) prefix_tsv="${OPTARG}" ;;
    p) pseudofinder_outdir="${OPTARG}" ;;
    f) orthofinder_downstream_dir="${OPTARG}" ;;
    *) print_usage
       exit 1 ;;
  esac
done


## Make outdir and log file
mkdir $outdir
cd $outdir
log=${outdir}/assign_pseudogenes_2_hogs.log.txt
> $log

## Make .dmnd db of intact protein sequences!
## I am not using the db that I used for pseudofinder because this included sequences that would be predicted pseudogenes
## I debated making a unique db for each genome that does not include NosAzo0708, the MAGs, and the genome itself
## In order to parallel the pseudofinder dbs
## But I don't think this is necessary as it would be valid for a pseudogene to hit a MAG intact gene.
echo "Make diamond db of intact protein sequences" >> $log
db=all_intact.faa
> $db
cat ${prefix_tsv} | \
while read long_prefix short_prefix; do
cat ${pseudofinder_outdir}/${long_prefix}/${long_prefix}_intact.faa >> $db
done
echo "Check that one prefix per header line. All three values should be the same." >> $log
## First get a list of just the short prefixes
cat $prefix_tsv | awk '{print $2}' > short_prefix_list.tmp
grep -c -Ff short_prefix_list.tmp $db >> $log
grep -c "^>" $db >> $log
grep "^>" $db | grep -c -Ff short_prefix_list.tmp >> $log
rm short_prefix_list.tmp
## Make the db!
conda deactivate
conda activate orthofinder
diamond makedb --in $db -d ${db}.dmnd
db=${db}.dmnd

## Confirm that each intact locus appears only once in ${N}_noOutgroups.tsv
## I believe I do this in another script as well.
cd ${orthofinder_downstream_dir}
tail -n +2 ${N}_noOutgroups.tsv | \
awk -F "\t" ' { for (i=4;i<=NF;++i) { split($i,A,", "); { for (j=1;j<=length(A);++j) {print A[j]} } } } ' \
> tmp
echo "The following two numbers should be equal if loci are all unique in ${N}_noOutgroups.tsv:" >> $log
wc -l < tmp >> $log
sort -u tmp | wc -l >> $log
rm tmp

## Combine all the pseudogenes into one fasta for diamond query.
## This will allow a more efficient use of locus_hog_conversion_file.tmp below
## Because an intact locus that hits multiple pseudogenes from multiple genomes will only have to be searched against the full ${N}_noOutgroups.tsv once
echo "Combine all pseudos into one fasta for query" >> $log
cd $outdir
> pseudos.fasta
cat ${prefix_tsv} | \
while read long_prefix short_prefix; do
cat ${pseudofinder_outdir}/${long_prefix}/${long_prefix}_pseudos.shorter_headers.fasta >> pseudos.fasta
done

## Run diamond
echo "Run diamond" >> $log
conda deactivate
conda activate orthofinder
diamond blastx --threads 20 --outfmt 6 --query pseudos.fasta --db $db --out pseudos_vs_intact.dmnd_out.tsv
conda deactivate

## Assign each hit an ${N}.HOG by matching the intact locus id hits from the dmnd search above
## To that locus and the corresponding HOG in locus_hog_conversion_file.tmp
echo "Generate *locus_hog_conversion_file.tmp files" >> $log
infile=pseudos_vs_intact.dmnd_out.tsv
## get a subset of ${N}_noOutgroups.tsv for locus_hog_conversion_file.tmp to speed up search
## Because several intact loci are hit by multiple pseudogene loci
## (I think an average of 5x per intact locus)
## The following subsetting lets us only search the full ${N}_noOutgroups.tsv file 1 time per intact locus
## Get a list of unique intact loci tht were hit 
cat $infile | awk -F "\t" '{print $2}' | sort -u > loci_that_hit.tmp
## Make a separate conversion file for each prefix to aid in later search.
cat ${prefix_tsv} | \
while read long_prefix short_prefix; do

## subset loci_that_hit.tmp for this prefix knowing it is sorted
## grep might still be faster. Both are 0 on the timer.
cat loci_that_hit.tmp | \
awk -v prefix=$short_prefix -F "\n" '
BEGIN {p=0};
$0 ~ "^"prefix {p=1};
p==1 { if ($0 !~ "^"prefix) {exit} else {print $0} } 
' > subset_loci_that_hit.tmp

## subset ${N}_noOutgroups.tsv to only the HOG column and the column for this prefix
## And split each locus into its own row with the HOG (so can be multiple rows with same HOG) 
cat ${orthofinder_downstream_dir}/${N}_noOutgroups.tsv | \
awk -v prefix=$short_prefix -F "\t" '
NR == 1 { for (i=4;i<=NF;++i) { if ($i==prefix) {c=i; next} } };
{ split($c,A,", "); for (i=1;i<=length(A);++i) { print $1"\t"A[i] } }
' > subset.${N}.tsv.tmp

## Get overlap between subset.${N}.tsv.tmp and subset_loci_that_hit.tmp
grep -Ff subset_loci_that_hit.tmp subset.${N}.tsv.tmp > ${short_prefix}.locus_hog_conversion_file.tmp
#Add a row for each locus that isn't in ${N}
cat ${short_prefix}.locus_hog_conversion_file.tmp | awk '{print $2}' > ${N}.loci.tmp
grep -vFf ${N}.loci.tmp subset_loci_that_hit.tmp | awk '{print "-\t"$0}' >> ${short_prefix}.locus_hog_conversion_file.tmp
rm ${N}.loci.tmp

## Finish loop through prefixes
done

## Make sure that each of the conversion files has a unique locus on each line
## This allows the "exit" in the awk call that builds hog.col.tmp
> check_unique.tmp
for file in *locus_hog_conversion_file.tmp; do
a=$(wc -l < $file);
b=$(awk '{print $2}' $file | sort -u | wc -l);
echo $(( $a - $b )) >> check_unique.tmp;
done
echo "the following should be just one 0 if *locus_hog_conversion_file.tmp all have unique loci on each line:" >> $log
sort -u check_unique.tmp >> $log
#Just 0! success
rm check_unique.tmp

## initiate the list of hog matches. Will paste this to the list of pseudos.
echo "match pseudos to HOGs. hog.col.tmp should grow in length (wc -l) to equal length of pseudos_vs_intact.dmnd_out.tsv." >> $log
echo "May take a long time (on the order of 24hrs for 48 input genomes)" >> $log
> hog.col.tmp

## Loop through the diamond output
cat $infile | \
## Get the intact locus
awk -F "\t" '{print $2}' | \
while read locus; do
prefix=${locus%_*}
## get the matching line from the hog table
## Keep both the intact locus and the hog in case helpful to have intact locus downstream
cat ${prefix}.locus_hog_conversion_file.tmp | awk -v locus=$locus ' $2 == locus {print $2"\t"$1; exit}' >> hog.col.tmp
#Finish loop through loci
done

## Paste pseudo.col and hog.col to output file
cat $infile | awk -F "\t" '{print $1}' > pseudo.col.tmp
paste pseudo.col.tmp hog.col.tmp > pseudos_blastHits_HOGs.tsv

## Assign pseudos to HOGs
## pseudos_N9HOGs_greatest_portion_hits.tmp: for each locus, get the hog with the highest portion of hits (as well as that portion)
cat pseudos_blastHits_HOGs.tsv | awk '{print $1"\t"$3}' | \
sort | uniq -c | sort -k2,2 -k1,1nr | \
awk '
BEGIN {
        previous_locus=""; previous_hog=""; c=1; n=1;
};
$2 != previous_locus {
        print previous_locus"\t"previous_hog"\t"n/c;
        previous_locus=$2; previous_hog=$3; n=$1; c=0;
}
{
        c+=$1
}
END {
        print previous_locus"\t"previous_hog"\t"n/c;
}
' | \
tail -n +2 \
> pseudos_N9HOGs_greatest_portion_hits.tmp

## pseudos_N9HOGs_top_hit.tmp: for each locus, get the hog with the top e-val hit (dmnd auto arranges by evalue)
cat pseudos_blastHits_HOGs.tsv | \
awk '
BEGIN { previous_locus="" };
$1 != previous_locus { print $1"\t"$3; previous_locus=$1 };
' | \
## sort to get in the same order as the portion of blasthits tsv
sort -k1,1 \
> pseudos_N9HOGs_top_hit.tmp

## Check that the two files are in the same order with respect to loci
echo "awk will print any lines that aren't correctly matched up for combining the top-hit and portion-hits tables:" >> $log
echo "should be empty (so awk prints nothing)" >> $log
paste pseudos_N9HOGs_greatest_portion_hits.tmp pseudos_N9HOGs_top_hit.tmp | awk '$1!=$4' >> $log
## Nothing returned! good!
## Make a table where each locus has the top hog by portion of hits and by top hit.
echo -e "pseudo\tHighestPortionHits\tPortionHits\tTopHit" > pseudos_best_HOGs.tsv
cat pseudos_N9HOGs_top_hit.tmp | awk '{print $2}' > tmp
paste pseudos_N9HOGs_greatest_portion_hits.tmp tmp >> pseudos_best_HOGs.tsv
rm *tmp

## Get some basic stats
echo -e "\
###################\n\
## Collect stats ##\n\
###################\
" >> $log
echo "How many pseudos?" >> $log
grep -c "^>" pseudos.fasta >> $log
echo "How many got at least one blast hit?" >> $log
cat pseudos_vs_intact.dmnd_out.tsv | awk '{print $1}' | sort -u | wc -l >> $log
echo "How many pseudos assigned had one hog be >50% of assignments?" >> $log
cat pseudos_best_HOGs.tsv | awk ' $3 > 0.5 ' | wc -l >> $log
echo "How about >50% of hits to the same hog as the top hit?" >> $log
cat pseudos_best_HOGs.tsv | awk ' $3 > 0.5 && $2==$4 ' | wc -l >> $log
echo "How about 100% hits to the same hog?" >> $log
cat pseudos_best_HOGs.tsv | awk ' $3 == 1 ' | wc -l >> $log

## Using the portion of hits that go to particular HOG might not be the best choice
## Because many HOGs are smaller than the max number of hits dmnd will return
## So if dmnd returns 10 hits, and the "correct" HOG is only 4 seqs,
## Then even if they are the top 4 hits, it will be inconclusive with a >0.5 cutoff
## So, below I look at how consistently the second hit is the same as the first hit.
## Also note that if e-val=0 for multiple hits then I am not sure how dmnd ranks
cat pseudos_blastHits_HOGs.tsv | awk '
BEGIN { previous_locus="" ; second=0 };
$1 != previous_locus { 
        if ( second==1 ) { print "-"}; #Means previous_locus only had one hit
        printf $1"\t"$3"\t"; previous_locus=$1; second=1; next;
}
second == 1 { second=0; print $3 };
END { if (second==1) {print "--"} }  
' | sort -k1,1 > top2hits.tmp
## Checked using wc -l that same length as pseudos_best_hogs.tsv and first cols appear to match
echo "How many pseudos have the same HOG for the top 2 hits?" >> $log
cat top2hits.tmp | awk '$3 == $2' | wc -l >> $log >> $log
echo "How many have different HOGs for the top 2 hits?" >> $log
cat top2hits.tmp | awk '$3 != $2' | wc -l >> $log
echo "How many have only one hit?" >> $log
cat top2hits.tmp | awk '$3 == "--"' | wc -l >> $log

## Can remove all the files except the output because recreating the files (other than the output) is quick
rm all_intact.faa all_intact.faa.dmnd pseudos.fasta pseudos_vs_intact.dmnd_out.tsv top2hits.tmp

##############################
### Pseudogene count table ###
##############################

##Make a table of hog x genome x number of pseudogenes
##This will give a 3 col tsv with col1=count, col2=prefix, col3=HOG
tail -n +2 pseudos_best_HOGs.tsv | awk -F "\t" '{print $1"\t"$4}' | sed 's/_p_.*\t/\t/' | sort | uniq -c | sort -k3,3 | awk '$3 != "-"' > pseudo_counts.tmp
## Create first col
echo "HOG" > first_col.tmp
cat pseudo_counts.tmp | awk '{print $3}' | sort -u >> first_col.tmp
#initiate the rest of the output file with a header
cat $prefix_tsv | awk '{print $2}' | tr "\n" "\t" | awk '{print $0}' | sed 's/\t$//' > body.tmp

##Loop through HOGs
tail -n +2 first_col.tmp | \
while read hog; do
#Get the lines corresponding to this hog. Could just grep instead.
cat pseudo_counts.tmp | awk -v hog=$hog 'BEGIN {p=0}; $3==hog {p=1; print $0; next}; p==1 {exit}' > this_block.tmp
#Create the output line for this hog
> this_line.tmp
#Loop through prefixes
cat $prefix_tsv | \
while read long_prefix short_prefix; do
cat this_block.tmp | awk -v prefix=$short_prefix '$2 == prefix {print $1; e="no_end"; exit}; END {if (e != "no_end") {print "0"}}' >> this_line.tmp
done
##Make this_line.tmp horizontal and add to the output file
cat this_line.tmp | tr "\n" "\t" | awk '{print $0}' | sed 's/\t$//' >> body.tmp
rm this_line.tmp 
done
#paste the output file together
paste first_col.tmp body.tmp > n_pseudo.tsv
rm *tmp

##Now make a table combining the intact and pseudogene counts for HOG x genome
#Key: #intact/#pseudo
cd $orthofinder_downstream_dir
outfile=n_intact_and_pseudo.tsv
pseudo_infile=${outdir}/n_pseudo.tsv
intact_infile=n_intact.tsv

##Going to take the congruent values from each of the intact and pseudotables
##So want to make sure the columns are in the same order.
cat $intact_infile | awk -F "\t" '{printf $4; for (i=5; i<=NF; ++i) {printf "\t"$i}; print ""}' | head -n 1 > intact_header.tmp
cat $pseudo_infile | awk -F "\t" '{printf $2; for (i=3; i<=NF; ++i) {printf "\t"$i}; print ""}' | head -n 1 > pseudo_header.tmp
## If not, then rearrange so they are.
if ! cmp -s intact_header.tmp pseudo_header.tmp; then
cat $pseudo_infile | awk -F "\t" '{print $1}'> tmp
head -n 1 $intact_infile | awk -F "\t" '{ for (i=4;i<=NF;++i) {print $i} }' | \
while read header; do
echo "$header"
cat $pseudo_infile | awk -v header="$header" -F "\t" '
NR == 1 { for (i=1;i<=NF;++i) { if ($i==header) {c=i} } };
{print $c}
' > this_col.tmp
paste tmp this_col.tmp > tmp.tmp
mv tmp.tmp tmp
done
rm this_col.tmp
mv tmp $pseudo_infile
fi


#Cat the two tables together (without headers and row labels)
#And then sort so corresponding HOGs are on consecutive lines
#Then combine identical lines
#and print the output.
tail -n +2 $intact_infile | awk -F "\t" '{printf $1; for (i=4; i<=NF; ++i) {printf "\t"$i}; print ""}' > $outfile
#Add a -p to the HOG name for pseudogenes to make sure they come second after sort.
#Cat together
tail -n +2 $pseudo_infile | sed 's/\(^[^\t]*\)\t/\1-p\t/' >> $outfile
#sort
sort -k1,1 $outfile > tmp
mv tmp $outfile

#Go through, saving the previous line.
#If next line isn't the same HOG, then print the previous line with /0s for 0 pseudos.
#If it is the same, then print the intact (previous) line / current (pseudo) line values
cat $outfile | \
awk -F "\t" '
BEGIN {fresh=1};
{
if (fresh == 1) {fresh=0; split($0,prev,"\t"); next};
if ($1 == prev[1]"-p") {
        fresh=1; printf prev[1]; for (i=2; i<=NF; ++i) {printf "\t"prev[i]"/"$i }; print ""
}
else {
        printf prev[1]; for (i=2; i<=length(prev); ++i) {printf "\t"prev[i]"/0"}; print "";
        split($0,prev,"\t");
        next
}
}
END { if (fresh==0) {printf prev[1]; for (i=2; i<=length(prev); ++i) {printf "\t"prev[i]"/0"}; print ""} }
' \
> tmp
mv tmp $outfile

#Give it a header
head -n 1 $pseudo_infile > tmp
cat tmp $outfile > tmp.tmp
mv tmp.tmp $outfile
rm *tmp*

#############################
#### Split table into 2 #####
## with Taz and Nostocales ##
######## cols separate ######
#############################

cd $orthofinder_downstream_dir
intable=n_intact_and_pseudo.tsv

## Do the Taz cols, then the Nostocales
echo -e "\
Taz.n_intact_and_pseudo.tsv\t~/Taz.short_prefix.list
free_living.n_intact_and_pseudo.tsv\t~/subtree_freeLiving.short_prefix.list\
" | \
while read outtable prefix_list; do

#Print the hog column
cat $intable | awk -F "\t" '{print $1}' > $outtable

#Loop through the prefixes, pasting each column to the existing outtable
cat $prefix_list | \
while read prefix; do
cat $intable | \
awk -v prefix=$prefix -F "\t" '
NR==1 { for (i=2;i<=NF;++i) { if ($i==prefix) {c=i} }  };
{print $c}
' > tmp
paste $outtable tmp > tmp.tmp
mv tmp.tmp $outtable 

## Finish loop through prefixes
done
rm *tmp*

## Finish loop through Taz vs. Nostoc
done

## clean up some files
rm body.tmp first_col.tmp pseudo_counts.tmp
