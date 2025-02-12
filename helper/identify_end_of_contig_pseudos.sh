## Some pseudogenes maybe be inaccurately predicted because sequences are broken due to contig ends.
## this helper script identifies such pseudogenes and adds them as a third category to the 
## n_intact_pseudo table to be n_intact_pseudo_w_end_of_contig.tsv
## It takes a couple of variables so that it can be called from the main_pipeline and from the improve pipeline
## But quite a bit is hardcoded if I were to try to use it for another purpose.

## Set variables
#pseudofinder_dir=~/data/main_pipeline/pseudofinder_out
#assign_pseudos_dir=~/data/main_pipeline/assign_pseudogenes
#orthofinder_downstream_dir=~/data/main_pipeline/orthofinder_out/for_downstream


## Get variables
print_usage() {
  printf "usage requires:\n\
pseudofinder_dir -p\n\
assign_pseudos_dir -a\n\
orthofinder_downstream_dir -o\
"
}

while getopts 'p:a:o:' flag; do
  case "${flag}" in
    p) pseudofinder_dir="${OPTARG}" ;;
    a) assign_pseudos_dir="${OPTARG}" ;;
    o) orthofinder_downstream_dir="${OPTARG}" ;;
    *) print_usage
       exit 1 ;;
  esac
done



##################################
### Identify predicted-pseudos ###
#### that run off the end of #####
########### the contig ###########
##################################

cd ${pseudofinder_dir}
outfile=${pseudofinder_dir}/end_of_contig_pseudos.list
logfile=${pseudofinder_dir}/log.tmp
> $outfile
> $logfile

## Loop through
cat ~/subtree.shortened.tsv | \
while read long_prefix short_prefix; do

cd ${pseudofinder_dir}/${long_prefix}

## Get the contig lengths
## Could get this from inputs to pseudofinder (i.e. spades assemblies)
## But simpler to do it within the pseudofinder folder in case inputs are from different places
## (which they are for T. az. vs. free-living)
#cat ${long_prefix}_contigs.fasta | \
#awk ' NR%2==1 {printf $0"\t"; next};  {print length($0)} ' | \
#sed 's/^>//' \
#> contig_lengths.tmp

## Get the indices for each end of each contig that is outside of the last ORF
cat ${long_prefix}_cds.fasta | grep "^>" | tr -dc '[:alnum:][:space:]_:' | awk '{print $2"\t"$3}' | tr ":" "\t" | sort -k1,1 -k2,2n > tmp
echo "check that second index always greater than first index for ${short_prefix}. Should be blank:" >> $logfile
cat tmp | awk '$2 > $3' >> $logfile
cat tmp | \
awk '
BEGIN {contig="";first_start="";last_stop=""};
$1!=contig { print contig"\t"first_start"\t"last_stop; contig=$1; first_start=$2 };
{last_stop=$3};
END { print contig"\t"first_start"\t"last_stop }
' | \
tail -n +2 \
> contig_ends.tsv

## Get the intergenic pseudogenes. If a pseudogene has a closed ORF,
## I am not worried about it being at the very beginning or end of a contig.
## Some intergenic regions are included in "fragmented" pseudogenes,
## BUt I think I am also not worried about these because there is a STOP somewhere inside the predicted gene
cat ${long_prefix}_pseudos.gff | \
grep -v "^#" | \
grep "Intergenic region" \
> tmp

## Make a 4-col tsv of the locus_tag, contig, start, stop
cat tmp | sed 's/old_locus_tag=//' | sed 's/.*locus_tag=//' | sed 's/;.*//' > first_col.tmp
## make the locus_tags the same as in shorter_headers
sed -i 's/_[^_]\+_pseudo/_p/' first_col.tmp
cat tmp | awk '{print $1"\t"$4"\t"$5}' > other_cols.tmp
paste first_col.tmp other_cols.tmp > starts_stops_contigs.tmp
rm first_col.tmp other_cols.tmp tmp
echo "Check that all stops are higher than starts for ${long_prefix} (should be empty):" >> $logfile
awk '$3>$4' starts_stops_contigs.tmp >> $logfile


## list pseudogenes that have 1 or the length of the contig as start or stop.
## Loop through contigs
cat contig_ends.tsv | \
while read contig first_start last_stop; do
cat starts_stops_contigs.tmp | \
awk -v contig=$contig -v first_start=$first_start -v last_stop=$last_stop '
{
if ($2==contig) {
        if ( $3<first_start || $4>last_stop ) {print $1}
}
}
' >> ${outfile}
done

rm contig_ends.tsv starts_stops_contigs.tmp
## End loop through prefixes
done

## Check the logfile
grep -vi "^check" $logfile
## empty!
rm $logfile

## Add a slash value to the n_intact_pseudo table
cd ${pseudofinder_dir}
cat end_of_contig_pseudos.list ${assign_pseudos_dir}/pseudos_best_HOGs.tsv | \
sort -k1,1 | \
awk '
BEGIN {prev=""};
$1==prev {print $1"\t"$4}
{prev=$1}
' | \
sed 's/_p_.*\t/\t/' | \
sort | \
uniq -c | \
sort -k3,3 | \
awk '$3 != "-"' \
> pseudo_counts.tmp

prefix_tsv=~/subtree.shortened.tsv
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
paste first_col.tmp body.tmp > n_end_of_contig_pseudo.tsv
rm *tmp

## Before adding a new slash to n_intact_and_pseudo, make sure columns are in the same order
old_table=${orthofinder_downstream_dir}/n_intact_and_pseudo.tsv
> rearrange_table.tmp
head -n 1 $old_table | tr "\t" "\n" | \
while read header; do
cat n_end_of_contig_pseudo.tsv | \
awk -v header=$header -F "\t" ' NR==1 { for (i=1;i<=NF;++i) { if ($i==header) {c=i} } }; {print $c} ' \
> this_col.tmp
paste rearrange_table.tmp this_col.tmp > tmp
mv tmp rearrange_table.tmp
done
sed 's/^\t//' rearrange_table.tmp > n_end_of_contig_pseudo.tsv
rm rearrange_table.tmp
## Check the order now the same
head -n 1 $old_table > tmp
head -n 1 n_end_of_contig_pseudo.tsv > tmp.tmp
cmp tmp tmp.tmp
rm tmp tmp.tmp
## Check successful.
#Cat the two tables together (without headers and row labels)
#And then sort so corresponding HOGs are on consecutive lines
#Then combine identical lines
#and print the output.
outfile=${orthofinder_downstream_dir}/n_intact_pseudo_w_end_of_contig.tsv
end_of_contig_file=${pseudofinder_dir}/n_end_of_contig_pseudo.tsv
tail -n +2 $old_table > $outfile
#Add a -p to the HOG name for pseudogenes to make sure they come second after sort.
#Cat together
tail -n +2 $end_of_contig_file | sed 's/\(^[^\t]*\)\t/\1-p\t/' >> $outfile
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
head -n 1 $end_of_contig_file > tmp
cat tmp $outfile > tmp.tmp
mv tmp.tmp $outfile
rm *tmp*

## Split into Taz and free-living cols
cd ${orthofinder_downstream_dir}
intable=n_intact_pseudo_w_end_of_contig.tsv
## Do the Taz cols, then the Nostocales
echo -e "\
Taz.${intable}\t~/Taz.short_prefix.list
free_living.${intable}\t~/subtree_freeLiving.short_prefix.list\
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

