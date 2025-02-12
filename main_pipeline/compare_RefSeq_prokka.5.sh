## Run the early pipeline, through orthofinder
## Using both the prokka ORF predictions and the RefSeq ORF predictions
## The goal of this is to determine how much the gene prediction step, which we cannot ground truth, affects the results.
## Comments are more minimal in this doc, so see the main prokka/pseudofinder/orthofinder pipelines


###################################
######## RefSeq vs. prokka ########
########### annotations ###########
###################################

## In the main pipeline, the RefSeq genomes had ORF prediction redone with prokka
## In order to be consistent with the ORF predictions for the MAGs
## To determine how big of an effect this has
## I will here run orthofinder without the MAGs
## with both copies (RefSeq and prokka ORFs) of each RefSeq genome
## To see how well the pairs are self-consistent.
## The RefSeq ORFs must first be run through pseudofinder
## which will use the RefSeq ORFs as the database instead of the prokka database.

cd ~/sandbox
mkdir ~/data/main_pipeline/RefSeq_vs_prokka
outdir=~/data/main_pipeline/RefSeq_vs_prokka/pseudofinder_out
mkdir $outdir
indir=~/data/refseq_Nostocales_genomes/data
conda deactivate
conda activate pseudofinder

## Pseudofinder
## Loop through the 48 Nostocales
cat ~/Nostocales_incl0708.shortened.tsv | \
while read query short_query; do
## Make output folder
mkdir ${outdir}/$query
cd ${outdir}/$query
## Make db
db=${outdir}/$query/allNostocExceptTazAnd.${query}.faa
> $db
## Loop through the Nostocales again, adding the protein.faa to the database
## Except for the query protein.faa itself and Nostocales0708.
grep -v "^NostocAzollae" ~/Nostocales_incl0708.shortened.tsv  | \
while read prefix short_prefix; do
if [[ $prefix != $query ]]; then
cat ${indir}/${prefix}/protein.faa | \
awk -v prefix=$prefix -F " " ' {if ($0 ~ "^>") {printf $1"_"prefix" "; for (i=2; i<NF; ++i) {printf $i" "}; print $NF} else { print $0} }' \
>> $db;
fi
## End prefix loop for db build
done
in_gbk=${indir}/${query}/genomic.gbff
## Use use_deviation and -ce flags as in main pipeline
python ~/tools/pseudofinder/pseudofinder.py annotate -t 20 -ce --use_deviation -l 0.65 --diamond -skpdb -hc 46 --genome $in_gbk --outprefix ${query} --database $db | tee pseudofinder.RefSeq_vs_prokka.log;
## End query loop for pseudofinder run wrapper
done
## Change the headers in ${prefix}_pseudos.fasta to be shorter
## and easier to read and to contain the same locus prefix as in the intact RefSeq genomes
## Save in new file ${prefix}_pseudos.shorter_headers.fasta
cd $outdir
cat ~/Nostocales_incl0708.shortened.tsv  | \
while read long_prefix short_prefix; do
cat ${long_prefix}/${long_prefix}_pseudos.fasta | sed "s/^.*_pseudo/>${short_prefix}_RefSeq_p/" > ${long_prefix}/${long_prefix}_pseudos.shorter_headers.fasta
done


## orthofinder
## Run orthofinder using gloeobacter as the outgroup.
## Prepare input files for orthofinder!
## And rename input fastas to the shortened prefixes
conda deactivate
conda activate orthofinder
outdir=~/data/main_pipeline/RefSeq_vs_prokka/orthofinder_out
mkdir $outdir
indir=${outdir}/input_fastas
outdir=${outdir}/initial_run
mkdir $indir
## Get the free-living Nostocales fastas and rename
cat ~/Nostocales_incl0708.shortened.tsv | \
while read long_prefix short_prefix; do
## RefSeq
cp ~/data/main_pipeline/RefSeq_vs_prokka/pseudofinder_out/${long_prefix}/${long_prefix}_intact.faa ${indir}/${short_prefix}_RefSeq.fasta
## prokka
cp ~/data/main_pipeline/pseudofinder_out/${long_prefix}/${long_prefix}_intact.faa ${indir}/${short_prefix}_prokka.fasta
done
## Get the outgroup fastas and rename (these did not pass through pseudofinder)
cat ~/outgroup.shortened.tsv | \
while read long_prefix short_prefix; do
## RefSeq
cp ~/data/outgroup_refseq/${long_prefix}/protein.faa ${indir}/${short_prefix}_RefSeq.fasta
## prokka
cp ~/data/main_pipeline/prokka_out/${long_prefix}/${short_prefix}.faa ${indir}/${short_prefix}_prokka.fasta
done
## Run orthofinder!
orthofinder -y -M msa -o $outdir -f $indir | tee RefSeq_vs_Prokka.orthofinder.log

## Reroot species tree with gloeobacter as outgroup
conda deactivate
conda activate gotree
initial_orthofinder_outdir=~/data/main_pipeline/RefSeq_vs_prokka/orthofinder_out/initial_run/Results_Aug18
cd ${initial_orthofinder_outdir}/Species_Tree/
## make a list of the 4 gloeobacter outgroup genomes (2 genomes each annotated 2 ways)
cat ~/outgroup.short_prefix.list ~/outgroup.short_prefix.list > outgroup.prefix.list.tmp
echo -e "_RefSeq\n_RefSeq\n_prokka\n_prokka" > tmp
paste -d'\0' outgroup.prefix.list.tmp tmp > tmp.tmp
mv tmp.tmp outgroup.prefix.list.tmp
rm tmp
gotree reroot outgroup --strict -l outgroup.prefix.list.tmp -i SpeciesTree_rooted.txt -o SpeciesTree_RErooted.txt
cp SpeciesTree_RErooted.txt ${initial_orthofinder_outdir}/../..

## Re-run orthofinder with re-rooted tree
cd ${initial_orthofinder_outdir}/../..
conda deactivate
conda activate orthofinder
orthofinder -ft $initial_orthofinder_outdir -s SpeciesTree_RErooted.txt 
## Move and rename output folders
cd ${initial_orthofinder_outdir}/../..
mv initial_run/Results_Aug21 final_run
mv initial_run/Results_Aug18/* initial_run
rmdir initial_run/Results_Aug18
rm SpeciesTree_RErooted.txt

## From here down, look at how similar or divergent the results are for the prokka MAGs from the Refseq MAGs
orthofinder_outdir=~/data/main_pipeline/RefSeq_vs_prokka/orthofinder_out/final_run
orthofinder_downstream_dir=~/data/main_pipeline/RefSeq_vs_prokka/orthofinder_out/for_downstream
mkdir $orthofinder_downstream_dir
cd $orthofinder_downstream_dir

## Make a table with counts of how many loci from each genome are present in each N13.HOG
## For this analysis, I could alternatively have taken N1, which would have only eliminated the Gloeobacter outgroups
## But I want to keep it consistent with the main analyses and RefSeq_vs_prokka analyses
## So N13 corresponds to N9 and N8 in those analyses, respectively
outfile=N13_noOutgroups.tsv
infile=${orthofinder_outdir}/Phylogenetic_Hierarchical_Orthogroups/N13.tsv
##Remove columns corresponding to genomes that aren't in N13
##Get the names of the genomes in N13
awk -F "\t" '{print $1"\t"$2"\t"$3}' $infile > $outfile
tail -n +2 $infile | awk -F "\t" '{for (i=4;i<=NF;++i) {if ($i!="") {print i}}}' | sort -un > N13.cols.tmp
cat N13.cols.tmp | \
while read col; do
awk -v col="$col" -F "\t" '{print $col}' $infile > nextCol.tmp
paste $outfile nextCol.tmp > tmp
mv tmp $outfile
done

infile=N13_noOutgroups.tsv
head -n 1 $infile > header.tmp
tail -n +2 $infile | \
awk -F "\t" '{printf $1"\t"$2"\t"$3; for (i=4; i<=NF; ++i) {split($i,A,", "); printf "\t"length(A)}; print ""}'  \
> body.tmp
cat header.tmp body.tmp > n_intact.tsv
rm *.tmp


########################
## Assign pseudogenes ##
########################

## Make a prefix list (1st col long, 2nd col short prefix)
## Note that the long vs. short prefix thing is a bit funny for this RefSeq vs prokka pipeline so be careful!
cd ${orthofinder_downstream_dir}/../../
## Get the short prefixes
head -n 1 ${orthofinder_downstream_dir}/N13_noOutgroups.tsv | awk -F "\t" ' { for (i=4;i<=NF;++i) {print $i} } ' > N13.short_prefix.list
## Make this a double column to just use the short prefix as the long prefix for input to assign_pseudogenes.....sh
paste N13.short_prefix.list N13.short_prefix.list > N13.short_short.tsv
## Get the long prefixes for each of the short prefixes (there is overlap because of RefSeq and prokka having same long prefixes)
cat N13.short_prefix.list | sed 's/_.*//' | sort -u > short_prefixes.tmp
grep -Ff short_prefixes.tmp ~/All.shortened.tsv > N13.long_short.tsv
rm short_prefixes.tmp N13.short_prefix.list
## Double up the long prefix for each of the two short prefixes (RefSeq, prokka)
#cat N13.long_short.tsv | awk '{ print $1"\t"$2"_RefSeq"; print $1"\t"$2"_prokka" }' > tmp
#mv tmp N13.long_short.tsv

## Temporarily add the required pseudofinder outputs from the main pipeline
## to this pipeline's pseudofinder folder
## in the format of /pseudofinder_out/${long_prefix}/files
## so they can be input to assign_pseudogenes_to_hogs.sh
## Note that the long vs. short prefix thing is a bit funny for this RefSeq vs prokka pipeline so be careful!
this_pseudofinder_outdir=~/data/main_pipeline/RefSeq_vs_prokka/pseudofinder_out
main_pseudofinder_outdir=~/data/main_pipeline/pseudofinder_out
new_pseudofinder_outdir=~/data/main_pipeline/RefSeq_vs_prokka/tmp_pseudofinder_out
mkdir $new_pseudofinder_outdir
cat N13.long_short.tsv | \
while read long_prefix short_prefix; do
mkdir ${new_pseudofinder_outdir}/${short_prefix}_prokka
cp ${main_pseudofinder_outdir}/${long_prefix}/${long_prefix}_pseudos.shorter_headers.fasta \
${new_pseudofinder_outdir}/${short_prefix}_prokka/${short_prefix}_prokka_pseudos.shorter_headers.fasta
cp ${main_pseudofinder_outdir}/${long_prefix}/${long_prefix}_intact.faa \
${new_pseudofinder_outdir}/${short_prefix}_prokka/${short_prefix}_prokka_intact.faa
mkdir ${new_pseudofinder_outdir}/${short_prefix}_RefSeq
cp ${this_pseudofinder_outdir}/${long_prefix}/${long_prefix}_pseudos.shorter_headers.fasta \
${new_pseudofinder_outdir}/${short_prefix}_RefSeq/${short_prefix}_RefSeq_pseudos.shorter_headers.fasta
cp ${this_pseudofinder_outdir}/${long_prefix}/${long_prefix}_intact.faa \
${new_pseudofinder_outdir}/${short_prefix}_RefSeq/${short_prefix}_RefSeq_intact.faa
done

## Unlike in the main pipeline and Li_v_Dijk pipeline, this pipeline includes refseq genomes
## which do not have the genome prefixes in the feature headers
## I change that hear so downstream code will work properly.
## They have already been changed for the pseudogenes above so only intact need to change
cd $new_pseudofinder_outdir
cat ../N13.short_short.tsv | awk '{print $1}' | grep "RefSeq" | \
while read prefix; do
cat ${prefix}/${prefix}_intact.faa | sed -E "s/^>[^[:space:]]+[^0-9]([0-9]+[[:space:]])/>${prefix}_\1/" > tmp
mv tmp ${prefix}/${prefix}_intact.faa
done
## Also add "prokka" to the prefixes in headers for the prokka-annotated genomes
cat ../N13.long_short.tsv | awk '{print $2}' | \
while read prefix; do
cat ${prefix}_prokka/${prefix}_prokka_intact.faa | sed "s/${prefix}/${prefix}_prokka/" > tmp
mv tmp ${prefix}_prokka/${prefix}_prokka_intact.faa
cat ${prefix}_prokka/${prefix}_prokka_pseudos.shorter_headers.fasta | sed "s/${prefix}/${prefix}_prokka/" > tmp
mv tmp ${prefix}_prokka/${prefix}_prokka_pseudos.shorter_headers.fasta
done

## Now need to also update N13_noOutgrpus.tsv
cd ~/sandbox
infile=${orthofinder_downstream_dir}/${N}_noOutgroups.tsv
cp $infile ${orthofinder_downstream_dir}/${N}_noOutgroups_original.tsv
## Keep first three columns as are:
cat $infile | awk -F "\t" '{print $1"\t"$2"\t"$3}' > tmp
##Loop through rest of columns, modifying as above and pasting to tmp
col=4
nf=$(head -n 1 $infile | awk -F "\t" '{print NF}')
while (( $col <= $nf )); do
cat $infile | awk -F "\t" -v col=$col '{print $col}' > old_col.tmp
prefix=$(head -n 1 old_col.tmp)
echo $prefix > new_col.tmp
old_prefix=$(tail -n +2 old_col.tmp | awk '$0 !="" {print $0; exit}' | sed 's/,.*//' | sed -E "s/[0-9]+$//")
tail -n +2 old_col.tmp | sed "s/${old_prefix}/${prefix}_/g"  >> new_col.tmp
paste tmp new_col.tmp > tmp.tmp
mv tmp.tmp tmp
col=$(($col+1))
done
mv tmp $infile
rm *



## Inputs:
outdir=~/data/main_pipeline/RefSeq_vs_prokka/assign_pseudogenes
N=N13
prefix_tsv=~/data/main_pipeline/RefSeq_vs_prokka/N13.short_short.tsv
pseudofinder_outdir=~/data/main_pipeline/RefSeq_vs_prokka/tmp_pseudofinder_out
orthofinder_downstream_dir=~/data/main_pipeline/RefSeq_vs_prokka/orthofinder_out/for_downstream

## Call assign_pseudogenes_to_hogs.sh !
bash -i ~/scripts/helper/assign_pseudogenes_to_hogs.sh \
-p $pseudofinder_outdir \
-f $orthofinder_downstream_dir \
-n $N \
-l $prefix_tsv \
-o $outdir

## remove redundant pseudofinder outputs
rm -r $new_pseudofinder_outdir

##############
## Analyses ##
##############

## A couple things to do before the loop below:
mkdir ~/data/main_pipeline/RefSeq_vs_prokka/orthofinder_out/for_downstream/aminoAcid_HOGs
mkdir ~/data/main_pipeline/RefSeq_vs_prokka/orthofinder_out/for_downstream/nucleotide_HOGs
## Loop through the genome pairs
RefSeq_prokka_analysis_dir=~/data/main_pipeline/RefSeq_vs_prokka/analyses
mkdir $RefSeq_prokka_analysis_dir
orthofinder_downstream_dir=~/data/main_pipeline/RefSeq_vs_prokka/orthofinder_out/for_downstream
## Create a list of the genomes in N13
cd $RefSeq_prokka_analysis_dir
infile=~/data/main_pipeline/RefSeq_vs_prokka/orthofinder_out/final_run/Phylogenetic_Hierarchical_Orthogroups/N13.tsv
> N13_prefixes.tmp
head -n 1 $infile | \
awk -F "\t" '{ for (i=4;i<=NF;++i) {print $i} }' | \
grep 'RefSeq' | \
while read prefix; do
cat $infile | \
awk -v prefix=$prefix -F "\t" '
NR==1 { for (i=4;i<=NF;++i) { if ($i==prefix) {c=i; next} } };
$c != "" { print prefix; exit }
' \
>> N13_prefixes.tmp
done
sed -i 's/_RefSeq$//' N13_prefixes.tmp
grep -Ff N13_prefixes.tmp ~/Nostocales_incl0708.shortened.tsv > N13.shortened.tsv
rm N13_prefixes.tmp

## Loop through the genome pairs
cat ${RefSeq_prokka_analysis_dir}/N13.shortened.tsv | \
while read long_prefix short_prefix; do
refseq=${short_prefix}_RefSeq
prokka=${short_prefix}_prokka
cd ${RefSeq_prokka_analysis_dir}
mkdir $short_prefix
cd $short_prefix
log=${RefSeq_prokka_analysis_dir}/$short_prefix/comparison.log.txt
> $log

## First section compares basic genome stats
pseudofinder_refseq_dir=~/data/main_pipeline/RefSeq_vs_prokka/pseudofinder_out/$long_prefix
pseudofinder_prokka_dir=~/data/main_pipeline/pseudofinder_out/$long_prefix
echo "~~ First, let's compare some basic stats between the two genomes. ~~" >> $log
echo "" >> $log
echo "~~ How many total genes are predicted after prokka and pseudofinder ~~" >> $log
r_total=$(( $(grep -c "^>" ${pseudofinder_refseq_dir}/${long_prefix}_intact.faa) + \
$(grep -c "^>" ${pseudofinder_refseq_dir}/${long_prefix}_pseudos.fasta) ))
p_total=$(( $(grep -c "^>" ${pseudofinder_prokka_dir}/${long_prefix}_intact.faa) + \
$(grep -c "^>" ${pseudofinder_prokka_dir}/${long_prefix}_pseudos.fasta) ))
echo "~~ RefSeq: ~~" >> $log
echo $r_total >> $log
echo "~~ prokka: ~~" >> $log
echo $p_total >> $log
echo "" >> $log
echo "~~ How many of those are predicted intact (remainder are predicted pseudogenes)? ~~" >> $log
r_intact=$(grep -c "^>" ${pseudofinder_refseq_dir}/${long_prefix}_intact.faa)
p_intact=$(grep -c "^>" ${pseudofinder_prokka_dir}/${long_prefix}_intact.faa)
echo "~~ RefSeq: ~~" >> $log
echo $r_intact >> $log
echo "~~ prokka: ~~" >> $log
echo $p_intact >> $log
echo "" >> $log
echo "~~ So what percent of the total genes are predicted intact? ~~" >> $log
r_percent_intact=$( awk -v n=$r_intact -v d=$r_total ' BEGIN { if (n==0) {print "0"} else {printf "%.1f\n", 100*n/d} } ' )
p_percent_intact=$( awk -v n=$p_intact -v d=$p_total ' BEGIN { if (n==0) {print "0"} else {printf "%.1f\n", 100*n/d} } ' )
echo "~~ RefSeq: ~~" >> $log
echo $r_percent_intact >> $log
echo "~~ prokka: ~~" >> $log
echo $p_percent_intact >> $log
echo "" >> $log

## Second section
## Here, I am checking how similar the gene sets between the two are
## as determined by orthofinder-derived N13.HOGs.
## Get the appropriate columns from the gene count N13.HOG x genome table
## And narrow it down to only the rows that aren't 0's
cat ${orthofinder_downstream_dir}/n_intact.tsv | \
awk -v refseq=$refseq -v prokka=$prokka -F "\t" '
NR==1 {for (i=1; i<=NF; ++i) {
if ($i=="HOG") {A[1]=i}
else if ($i==refseq) {A[2]=i}
else if ($i==prokka) {A[3]=i} } };
{print $A[1]"\t"$A[2]"\t"$A[3]}
' | \
awk -F "\t" '$2 != 0 || $3 != 0 {print $0}' \
> comparison_cols.tsv
## Get the differences
cat comparison_cols.tsv | \
awk -F "\t" '$2 != $3 {print $0}' > diffs.tsv

echo "~~ How similar are the orthofinder HOG profiles of these two genomes? ~~" >> $log
echo "" >> $log
echo "~~ How many HOGs are present in each genome? ~~" >> $log
r_n_hogs=$( tail -n +2 comparison_cols.tsv | awk ' $2 != 0 ' | wc -l )
p_n_hogs=$( tail -n +2 comparison_cols.tsv | awk ' $3 != 0 ' | wc -l )
echo "~~ RefSeq: ~~" >> $log
echo $r_n_hogs >> $log
echo "~~ prokka: ~~" >> $log
echo $p_n_hogs >> $log
echo "" >> $log
echo "~~ How many HOGs are present in both genomes? ~~" >> $log
n_shared_hogs=$( tail -n +2 comparison_cols.tsv | awk ' $2 != 0 && $3 != 0 ' | wc -l )
echo $n_shared_hogs >> $log
echo "" >> $log
echo "~~ So what percent of the HOGs from each genome are present in the other genome? ~~" >> $log
r_p_shared_hogs=$( awk -v n=$n_shared_hogs -v d=$r_n_hogs ' BEGIN { if (n==0) {print "0"} else {printf "%.1f\n", 100*n/d} } ' )
p_p_shared_hogs=$( awk -v n=$n_shared_hogs -v d=$p_n_hogs ' BEGIN { if (n==0) {print "0"} else {printf "%.1f\n", 100*n/d} } ' )
echo "~~ RefSeq: ~~" >> $log
echo $r_p_shared_hogs >> $log
echo "~~ prokka: ~~" >> $log
echo $p_p_shared_hogs >> $log
echo "" >> $log
echo "~~ Of those shared HOGs, what percent have a different number of copies in the two genomes? ~~" >> $log
n_diff_counts=$( tail -n +2 diffs.tsv | awk ' $2 != 0 && $3 != 0 ' | wc -l )
p_diff_counts=$( awk -v n=$n_diff_counts -v d=$n_shared_hogs ' BEGIN { if (n==0) {print "0"} else {printf "%.1f\n", 100*n/d} } ' )
echo $p_diff_counts >> $log
echo "" >> $log

## Third section: Curious the effects of differential pseudogene calls and high copy number hogs
## Look at the presence/absence N13.HOGs in the table with intact and pseudogene counts.
echo "~~ How much of an effect do differential pseudogene calls have on the differences above? ~~" >> $log
echo "" >> $log
cd ${RefSeq_prokka_analysis_dir}/$short_prefix
cat ${orthofinder_downstream_dir}/n_intact_and_pseudo.tsv | \
awk -v refseq=$refseq -v prokka=$prokka -F "\t" '
NR==1 {for (i=1; i<=NF; ++i) {
if ($i=="HOG") {A[1]=i}
else if ($i==refseq) {A[2]=i}
else if ($i==prokka) {A[3]=i} } };
{print $A[1]"\t"$A[2]"\t"$A[3]}
' > tmp
head -n 1 tmp > presAbsDiff_wPseudos.tsv
tail -n +2 diffs.tsv | awk ' $2 == 0 || $3 == 0 {print $1}' > presence_absence_diffs.tmp
grep -Ff presence_absence_diffs.tmp tmp >> presAbsDiff_wPseudos.tsv 
rm tmp presence_absence_diffs.tmp
#Split pseudos into separate column
head -n 1 presAbsDiff_wPseudos.tsv | awk -F "\t" '{print $1"\t"$2"\t"$2"_pseudo\t"$3"\t"$3"_pseudo"}' > tmp
tail -n +2 presAbsDiff_wPseudos.tsv | awk -F "[\t/]" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' >> tmp
mv tmp presAbsDiff_wPseudos.tsv
echo "~~ How many HOGs are unique to each genome (only counting predicted-intact)? ~~" >> $log
r_n_unique=$( tail -n +2 diffs.tsv | awk '$3==0' | wc -l )
p_n_unique=$( tail -n +2 diffs.tsv | awk '$2==0' | wc -l )
echo "~~ RefSeq: ~~" >> $log
echo $r_n_unique >> $log
echo "~~ prokka: ~~" >> $log
echo $p_n_unique >> $log
echo "" >> $log
echo "~~ Of those unique HOGs, what percent have more pseudogenes in the other genome, suggesting differential pseudogene calls likely account for this percent of presence/absence discrepancies for predicted intact? ~~" >> $log
r_n_morePseudos=$(cat presAbsDiff_wPseudos.tsv | awk -F "\t" '$5 > $3 && $4 == 0 {print $0}' | wc -l)
p_n_morePseudos=$(cat presAbsDiff_wPseudos.tsv | awk -F "\t" '$3 > $5 && $2 == 0 {print $0}' | wc -l)
r_p_morePseudos=$( awk -v n=$r_n_morePseudos -v d=$r_n_unique ' BEGIN { if (n==0) {print "0"} else {printf "%.1f\n", 100*n/d} } ' )
p_p_morePseudos=$( awk -v n=$p_n_morePseudos -v d=$p_n_unique ' BEGIN { if (n==0) {print "0"} else {printf "%.1f\n", 100*n/d} } ' )
echo "~~ Intact in RefSeq: ~~" >> $log
echo $r_p_morePseudos >> $log
echo "~~ Intact in prokka: ~~" >> $log
echo $p_p_morePseudos >> $log
echo "" >> $log

echo "~~ How much of an effect do high copy number HOGs have on the differences above? ~~" >> $log
echo "" >> $log
echo "~~ Of the HOGs that were unique (only counting intact) to one or the other genome, what percent had >4 predicted-intact loci from that genome? ~~" >> $log
tail -n +2 diffs.tsv | awk -F "\t" ' $2 == 0 || $3 == 0 ' > tmp
n_unique=$( wc -l < tmp )
n_unique_highCopy=$( cat tmp | awk -F "\t" ' $2 > 4 || $3 > 4 ' | wc -l )
p_unique_highCopy=$( awk -v n=$n_unique_highCopy -v d=$n_unique ' BEGIN { if (n==0) {print "0"} else {printf "%.1f\n", 100*n/d} } ' )
echo $p_unique_highCopy >> $log 
echo "" >> $log
echo "~~ How many HOGs are present and intact in both genomes, but differ in number of copies? (displayed earlier as a percent of all shared HOGs) ~~" >> $log
echo $n_diff_counts >> $log
echo "" >> $log
echo "~~ Of those, what percent had >4 loci in at least one genome? ~~" >> $log
n_diffCounts_highCopy=$( tail -n +2 diffs.tsv | awk ' ( $2 > 4 && $3 != 0 ) || ( $2 != 0 && $3 > 4 ) ' | wc -l )
p_diffCounts_highCopy=$( awk -v n=$n_diffCounts_highCopy -v d=$n_diff_counts ' BEGIN { if (n==0) {print "0"} else {printf "%.1f\n", 100*n/d} } ' )
echo $p_diffCounts_highCopy >> $log
echo "" >> $log
rm tmp

## Last section
## Compare the nucleotide and amino acid identities of seqs in HOGs where refseq and prokka each have 1 locus
echo "~~ Compare the nucleotide and amino acid seqs for HOGs with exactly 1 predicted-intact locus present from each genome ~~" >> $log
echo "" >> $log
tail -n +2 comparison_cols.tsv | \
awk -F "\t" ' $2==1 && $3==1 {print $1} ' \
> refseq_prokka_single_copy.list
echo "~~ How many HOGs have 1 copy in each genome? ~~" >> $log
wc -l < refseq_prokka_single_copy.list >> $log
echo "" >> $log


## loop through the hogs and grab the amino acid seqs
cd ${orthofinder_downstream_dir}/aminoAcid_HOGs
mkdir $short_prefix
cd $short_prefix
indir=~/data/main_pipeline/RefSeq_vs_prokka/orthofinder_out/input_fastas
tail -n +2 ${RefSeq_prokka_analysis_dir}/${short_prefix}/comparison_cols.tsv | \
awk -F "\t" '{print $1}' | \
while read hog; do
> ${hog}.faa
this_refseq_locus=$(\
cat ${orthofinder_downstream_dir}/N13_noOutgroups_original.tsv | \
awk -v hog=$hog -v refseq=$refseq -F "\t" '
NR==1 { for (i=4;i<=NF;++i) { if ( $i == refseq ) { c=i; next } } };
$1 == hog {print $c; exit}
'\
)
this_prokka_locus=$(\
cat ${orthofinder_downstream_dir}/N13_noOutgroups_original.tsv | \
awk -v hog=$hog -v prokka=$prokka -F "\t" '
NR==1 { for (i=4;i<=NF;++i) { if ( $i == prokka ) { c=i; next } } };
$1 == hog {print $c; exit}
'\
)
echo "$this_refseq_locus" | \
sed 's/, /\n/g' | \
while read locus; do
grep -A 1 $locus ${indir}/${refseq}.fasta >> ${hog}.faa
done
echo "$this_prokka_locus" | \
sed 's/, /\n/g' | \
while read locus; do
grep -A 1 $locus ${indir}/${prokka}.fasta >> ${hog}.faa
done
## FInish loop through hogs
done

## loop through hogs and grab the nucleotide seqs
cd ${orthofinder_downstream_dir}/nucleotide_HOGs
mkdir $short_prefix
cd $short_prefix
refseq_indir=~/data/main_pipeline/RefSeq_vs_prokka/pseudofinder_out/$long_prefix
prokka_indir=~/data/main_pipeline/pseudofinder_out/$long_prefix
tail -n +2 ${RefSeq_prokka_analysis_dir}/${short_prefix}/comparison_cols.tsv | \
awk -F "\t" '{print $1}' | \
while read hog; do
> ${hog}.ffn
this_refseq_locus=$(\
cat ${orthofinder_downstream_dir}/N13_noOutgroups_original.tsv | \
awk -v hog=$hog -v refseq=$refseq -F "\t" '
NR==1 { for (i=4;i<=NF;++i) { if ( $i == refseq ) { c=i; next } } };
$1 == hog {print $c; exit}
'\
)
this_prokka_locus=$(\
cat ~/data/main_pipeline/RefSeq_vs_prokka/orthofinder_out/for_downstream/N13_noOutgroups_original.tsv | \
awk -v hog=$hog -v prokka=$prokka -F "\t" '
NR==1 { for (i=4;i<=NF;++i) { if ( $i == prokka ) { c=i; next } } };
$1 == hog {print $c; exit}
'\
)
echo "$this_refseq_locus" | \
sed 's/, /\n/g' | \
while read locus; do
grep -A 1 $locus ${refseq_indir}/${long_prefix}_intact.ffn >> ${hog}.ffn
done
echo "$this_prokka_locus" | \
sed 's/, /\n/g' | \
while read locus; do
grep -A 1 $locus ${prokka_indir}/${long_prefix}_intact.ffn >> ${hog}.ffn
done
## FInish loop through hogs
done

## Determine the overlap between the single copy HOG seqs
## Can do this here but not in Li vs. Dijkhuizen because these are identical assemblies and those were not.
cd ${RefSeq_prokka_analysis_dir}/$short_prefix
outfile=single_copy_Seqs_not_identical_startStop.txt
echo -e "refseq_start\trefseq_stop\trefseq_strand\tprokka_start\tprokka_stop\tprokka_strand\trefseq_len\tprokka_len\toverlap_len\tp_refseq\tp_prokka" > $outfile
indir=~/data/main_pipeline/RefSeq_vs_prokka/orthofinder_out/for_downstream/nucleotide_HOGs/$short_prefix
cat refseq_prokka_single_copy.list | \
while read hog; do
grep "^>" ${indir}/${hog}.ffn | \
sed 's/^.*\[//' | \
tr ':(' '\t' | \
tr -d "])" | \
tr '\n' '\t' | \
sed 's/\t$/\n/' | \
awk -F "\t" '
## Check data how I expect
$2<$1 || $5<$4 {print $0"\tindices not in order"; exit};
## gene lengths
{ refseq_len=$2-$1+1; prokka_len=$5-$4+1 };
{
## No overlap or on different strands.
	if ($1>$5 || $4>$2 || $3!=$6 ) { overlap=0 }
## Yes overlap
	else {
## find start of overlap 
		if ($1>=$4) { overlap_start=$1 } else { overlap_start=$4 };
## find end of overlap
		if ($2<=$5) { overlap_end=$2 } else { overlap_end=$5 };
## Get output values
		overlap=overlap_end-overlap_start+1;
	}
}
## print!
{ print $0"\t"refseq_len"\t"prokka_len"\t"overlap"\t"overlap/refseq_len"\t"overlap/prokka_len }
' \
>> $outfile
done
## Add a first column with the hog names
echo "hog" > first_col.tmp
cat refseq_prokka_single_copy.list >> first_col.tmp
paste first_col.tmp $outfile > tmp
mv tmp $outfile
rm first_col.tmp

echo "~~ Check that the indices are always in the expected order (should be 0):~~ " >> $log
grep -c "indices not in order$" $outfile >> $log
echo "" >> $log
echo "~~ NOTE THAT IF THE INDICES OVERLAPPED BUT WERE ON DIFFERENT CHROMOSOME/PLASMID, THAT WOULD STILL BE COUNTED AS OVERLAP ~~" >> $log
echo "~~ BUT THIS SEEMS EXCEEDINGLY UNLIKELY TO US GIVEN HOW LITTLE SEQUENCE IS ON PLASMIDS ~~" >> $log
echo "~~ What percent have identical starts and stops? ~~" >> $log
n=$( tail -n +2 $outfile | awk -F "\t" '$NF==1 && $(NF-1)==1' | wc -l )
d=$( wc -l < refseq_prokka_single_copy.list )
o=$( awk -v n=$n -v d=$d ' BEGIN { if (n==0) {print "0"} else {printf "%.1f\n", 100*n/d} } ' )
echo "$o" >> $log
echo "" >> $log
echo "~~ What percent overlap at least 95% of both genes? ~~" >> $log
n=$( tail -n +2 $outfile | awk -F "\t" '$NF>=0.95 && $(NF-1)>=0.95' | wc -l )
o=$( awk -v n=$n -v d=$d ' BEGIN { if (n==0) {print "0"} else {printf "%.1f\n", 100*n/d} } ' )
echo "$o" >> $log
echo "" >> $log
echo "~~ What percent overlap at least 90% of both genes? ~~" >> $log
n=$( tail -n +2 $outfile | awk -F "\t" '$NF>=0.90 && $(NF-1)>=0.90' | wc -l )
o=$( awk -v n=$n -v d=$d ' BEGIN { if (n==0) {print "0"} else {printf "%.1f\n", 100*n/d} } ' )
echo "$o" >> $log
echo "" >> $log
echo "~~ How many (number not percent) have <95% overlap of at least one sequence and neither sequence is fully covered by the other? ~~" >> $log
n_poorOverlap=$( tail -n +2 $outfile | awk -F "\t" '($NF<0.95 || $(NF-1)<0.95) && $NF!=1 && $(NF-1)!=1' | wc -l )
echo $n_poorOverlap >> $log
echo "" >> $log
echo "~~ How many (number not percent) have 0 overlap?~~ " >> $log
n_noOverlap=$( tail -n +2 $outfile | awk -F "\t" '($NF==0 && $(NF-1)==0)' | wc -l )
echo $n_noOverlap >> $log
echo "" >> $log
echo "~~ Are the above the two values above the exact same set of genes? (1 if yes, 0 if no) ~~" >> $log
awk -v a=$n_poorOverlap -v b=$n_noOverlap ' BEGIN { if (a==b) {print "1"} else {print "0"} }' >> $log

## Finish loop through refseq-prokka genome pairs
done

## combine the output files into one.
outfile=comparison_table.tsv
cd ${RefSeq_prokka_analysis_dir}
> $outfile
cat ${RefSeq_prokka_analysis_dir}/N13.shortened.tsv | \
while read long_prefix short_prefix; do
paste $outfile ${short_prefix}/comparison.log.txt > tmp
mv tmp $outfile
done
sed -i 's/^\t\+//' $outfile
sed -i 's/^\(~~[^~]*~~\).*/\1/' $outfile

## And get the summarize the results in mean/median/min/max format for each quantity
mv $outfile comparison_table_allValues.tsv
## Get the number of values that should be in each data row.
## This is helpful for the awk arrays used below to make sure the length of the array stays consistent
## because I don't know if an array can be fully deleted once created.
n_vals=$( wc -l < ${RefSeq_prokka_analysis_dir}/N13.shortened.tsv )
echo "~~ All data lines are in mean, median, min, max format ~~" > $outfile
cat comparison_table_allValues.tsv | \
awk -v n_vals=$n_vals -F "\t" '
## If not a data line, just reprint
$0=="" || $0 ~ "^~~" {print $0; next};
## Data line! Make sure NF is as expected
NF != n_vals {print NF" Wrong number of values"; next};
## get mean, median, min, max
{
	min=$1; max=$1; running_sum=0
	for (i=1; i<=NF;++i) {
		if ($i < min) {min=$i};
		if ($i > max) {max=$i};
		A[i]=$i;
		running_sum+=$i
	}
	mean=running_sum/n_vals
	## asort sorts an array in place
	asort(A)
	if (n_vals%2==1) {
		median=A[(n_vals+1)/2]
	} else {
	median=( ( A[n_vals/2] + A[n_vals/2+1] ) / 2 )
	}
	##Print!
	print mean"\t"median"\t"min"\t"max;
}
' \
>> $outfile

## Finally! (Unless I add something else after this)
## Last few analyses

## Get the counts of how many genome pairs each hog differs in
cd ${RefSeq_prokka_analysis_dir}

## differences by count and then by pres/abs
echo -e "\
diffs.tsv,hog_countDiffs_counts.tsv\n\
presAbsDiff_wPseudos.tsv,hog_presAbsDiffs_counts.tsv\
" | \
sed 's/,/\t/g' | \
while read infile outfile; do

> tmp
cat N13.shortened.tsv | \
while read long_prefix short_prefix; do
tail -n +2 ${short_prefix}/$infile | awk -F "\t" '{print $1}' >> tmp
done
sort tmp | uniq -c | sort -rnk1,1 | awk '{print $2"\t"$1}' > $outfile
rm tmp

done

## Get the lengths of ORFs in different categories
## to test if this affects the consistency of the ORF-callers in any of the above metrics
cd ${RefSeq_prokka_analysis_dir}
outfile=hog_lengths.tsv

## Loop through genome pairs
cat N13.shortened.tsv | \
while read long_prefix short_prefix; do
cd ${RefSeq_prokka_analysis_dir}/${short_prefix}
echo -e "hog\nRefSeq_len\nprokka_len" > $outfile

## loop through hogs and get hog name, refseq len, prokka len
ls ~/data/main_pipeline/RefSeq_vs_prokka/orthofinder_out/for_downstream/nucleotide_HOGs/${short_prefix}/N13*.ffn | \
while read fasta; do

## hog name
hog=${fasta##*/}
hog=${hog%.ffn}
echo $hog >> $outfile
## refseq length
grep "^>" $fasta | \
grep -v $short_prefix > tmp
if (( $(wc -l < tmp) > 0 )); then
cat tmp | \
sed 's/^.*\[//' | \
sed 's/].*$//' | \
tr ':' '\t' | \
sed 's/\t$/\n/' | \
awk -F "\t" 'BEGIN {c=0}; {c+=($2-$1+1)}; END {print c/FNR}' \
>> $outfile
else
echo "-" >> $outfile
fi
## prokka length
grep "^>" $fasta | \
grep $short_prefix > tmp
if (( $(wc -l < tmp) > 0 )); then
cat tmp | \
sed 's/^.*\[//' | \
sed 's/].*$//' | \
tr ':' '\t' | \
sed 's/\t$/\n/' | \
awk -F "\t" 'BEGIN {c=0}; {c+=($2-$1+1)}; END {print c/FNR}' \
>> $outfile
else
echo "-" >> $outfile
fi

## Finish loop through hogs
done
rm tmp
## Make this a 3-col tsv
cat $outfile | awk 'NR%3==0 {print $0; next}; {printf $0"\t"}' > tmp
mv tmp $outfile

## Finish loop through genome pairs
done

## Loop could be part of above, but ran separately.
## Loop through prefixes
cat ${RefSeq_prokka_analysis_dir}/N13.shortened.tsv | \
while read long_prefix short_prefix; do
cd ${RefSeq_prokka_analysis_dir}/$short_prefix

## Get the lists of hogs for...
## the hogs with different counts
tail -n +2 diffs.tsv | awk '{print $1}' > diffs_hogs.list
## the hogs that differ in presence/absence
tail -n +2 presAbsDiff_wPseudos.tsv | awk '{print $1}' > presAbs_hogs.list
## the single-copy hogs that significantly differ in sequence
tail -n +2 single_copy_Seqs_not_identical_startStop.txt | \
awk '$NF<0.95 || $(NF-1)<0.95 {print $1}' > poorOverlap_hogs.list
## Loop through the lists and get the lengths for the in and out hog sets.
echo -e "\
diffs_hogs.list,countDiffs_length_ttest.tsv
presAbs_hogs.list,presAbs_length_ttest.tsv
poorOverlap_hogs.list,poorOverlap_length_ttest.tsv\
" | \
sed 's/,/\t/g' | \
while read in_list out_suffix; do

tail -n +2 hog_lengths.tsv | grep -Ff $in_list > diffs_lengths.tmp
tail -n +2 hog_lengths.tsv | grep -vFf $in_list > no_diffs_lengths.tmp
awk '{print $2}' diffs_lengths.tmp | grep -v "^-$" > col1.tmp
awk '{print $2}' no_diffs_lengths.tmp | grep -v "^-$" > col2.tmp
paste col1.tmp col2.tmp > refseq_${out_suffix}
awk '{print $3}' diffs_lengths.tmp | grep -v "^-$" > col1.tmp
awk '{print $3}' no_diffs_lengths.tmp | grep -v "^-$" > col2.tmp
paste col1.tmp col2.tmp > prokka_${out_suffix}
rm *tmp
rm $in_list

## Finish loop through gene sets
done

## Finish loop through prefixes
done


## Create paired ttest tables for total ORFs, intact, % intact, # HOGs
cd ${RefSeq_prokka_analysis_dir}
intable=comparison_table_allValues.tsv
grep -A 4 "^~~ How many total genes are predicted after prokka and pseudofinder ~~$" $intable | grep -v "^~~" > tmp
cat tmp | head -n 1 | tr "\t" "\n" > col1.tmp
cat tmp | tail -n 1 | tr "\t" "\n" > col2.tmp 
echo -e "RefSeq\tprokka" > n_orfs_ttest.tsv
paste col1.tmp col2.tmp >> n_orfs_ttest.tsv
grep -A 4 "^~~ How many of those are predicted intact (remainder are predicted pseudogenes)? ~~$" $intable  | grep -v "^~~" > tmp
cat tmp | head -n 1 | tr "\t" "\n" > col1.tmp
cat tmp | tail -n 1 | tr "\t" "\n" > col2.tmp 
echo -e "RefSeq\tprokka" > n_intact_ttest.tsv 
paste col1.tmp col2.tmp >> n_intact_ttest.tsv
grep -A 4 "^~~ So what percent of the total genes are predicted intact? ~~$" $intable | grep -v "^~~" > tmp
cat tmp | head -n 1 | tr "\t" "\n" > col1.tmp
cat tmp | tail -n 1 | tr "\t" "\n" > col2.tmp 
echo -e "RefSeq\tprokka" > p_intact_ttest.tsv
paste col1.tmp col2.tmp >> p_intact_ttest.tsv
grep -A 4 "^~~ How many HOGs are present in each genome? ~~$" $intable | grep -v "^~~" > tmp
cat tmp | head -n 1 | tr "\t" "\n" > col1.tmp
cat tmp | tail -n 1 | tr "\t" "\n" > col2.tmp 
echo -e "RefSeq\tprokka" > n_hogs_ttest.tsv
paste col1.tmp col2.tmp >> n_hogs_ttest.tsv
rm *tmp

## SCP to local computer to use Rstudio to run ttests and paired ttests!
## Used ttestMain.sh to call Ttests.R both of which are in ~/scripts/helper/ now

## How many pseudogenes does pseudofinder find compared with refseq
## This should be in the comparison_table.tsv loop,
## but i am doing it after the fact and don't want to have to edit that code.
RefSeq_prokka_analysis_dir=~/data/main_pipeline/RefSeq_vs_prokka/analyses
refseq_genome_dir=~/data/refseq_Nostocales_genomes/data
refseq_pseudofinder_dir=~/data/main_pipeline/RefSeq_vs_prokka/pseudofinder_out
prokka_pseudofinder_dir=~/data/main_pipeline/pseudofinder_out
outfile=pseudogene_counts_from_various_approaches.tsv
echo -e "prefix\trefseq-refseq\trefseq-pseudofinder\tprokka-pseudofinder" > $outfile
## Loop through prefixes
cat ${RefSeq_prokka_analysis_dir}/N13.shortened.tsv | \
while read long_prefix short_prefix; do
## Count the number of refseq-refseq-predicted pseudos
refseq_refseq_pseudos=$(
grep -c "RefSeq[[:space:]]\+pseudogene" ${refseq_genome_dir}/${long_prefix}/genomic.gff
)
## Count the number of refseq-pseudofinder-predicted pseudos
refseq_pseudofinder_pseudos=$(
grep "Pseudogenes (total):" ${refseq_pseudofinder_dir}/${long_prefix}/${long_prefix}_log.txt | \
sed 's/.*[[:space:]]//'
)
## Count the number of prokka-pseudofinder-predicted pseudos
prokka_pseudofinder_pseudos=$(
grep "Pseudogenes (total):" ${prokka_pseudofinder_dir}/${long_prefix}/${long_prefix}_log.txt | \
sed 's/.*[[:space:]]//'
)
## add everything to the outfile
echo -e "${short_prefix}\t${refseq_refseq_pseudos}\t${refseq_pseudofinder_pseudos}\t${prokka_pseudofinder_pseudos}" >> $outfile
done

## Look at the differences:
## For refseq-pseudofinder vs refseq-refseq then prokka-pseudofinder vs refseq-refseq
echo "refseq-pseudofinder,3
prokka-pseudofinder,4" | \
sed 's/,/\t/g' | \
while read comparison target_col; do
echo "${comparison} percent change in pseudogenes vs. refseq-refseq"
cat $outfile | \
tail -n +2 | \
awk -v col=$target_col '{print 100*($col-$2)/$2}' | \
sort -g \
> sorted.tmp
echo "min = "
head -n 1 sorted.tmp
echo "max = "
tail -n 1 sorted.tmp
echo "median = "
head -n 21 sorted.tmp | tail -n 1
echo "mean = "
cat sorted.tmp | awk 'BEGIN {c=0}; {c+=$0}; END {print c/FNR}'
rm sorted.tmp
done

## Compare the number of intact genes from Refseq-Refseq vs prokka-pseudofinder
## This should be in the comparison_table.tsv loop,
## but i am doing it after the fact and don't want to have to edit that code.
RefSeq_prokka_analysis_dir=~/data/main_pipeline/RefSeq_vs_prokka/analyses
refseq_genome_dir=~/data/refseq_Nostocales_genomes/data
refseq_pseudofinder_dir=~/data/main_pipeline/RefSeq_vs_prokka/pseudofinder_out
prokka_pseudofinder_dir=~/data/main_pipeline/pseudofinder_out
outfile=intact_counts_from_various_approaches.tsv
echo -e "prefix\trefseq-refseq\trefseq-pseudofinder\tprokka-pseudofinder" > $outfile
## Loop through prefixes
cat ${RefSeq_prokka_analysis_dir}/N13.shortened.tsv | \
while read long_prefix short_prefix; do
## Count the number of refseq-refseq-predicted intact
refseq_refseq=$(
grep -c "^>" ${refseq_genome_dir}/${long_prefix}/protein.faa
)
## Count the number of refseq-pseudofinder-predicted intact
refseq_pseudofinder=$(
grep "Intact genes:" ${refseq_pseudofinder_dir}/${long_prefix}/${long_prefix}_log.txt | \
head -n 1 | \
sed 's/.*[[:space:]]//'
)
## Count the number of prokka-pseudofinder-predicted intact
prokka_pseudofinder=$(
grep "Intact genes:" ${prokka_pseudofinder_dir}/${long_prefix}/${long_prefix}_log.txt | \
head -n 1 | \
sed 's/.*[[:space:]]//'
)
## add everything to the outfile
echo -e "${short_prefix}\t${refseq_refseq}\t${refseq_pseudofinder}\t${prokka_pseudofinder}" >> $outfile
done

## Look at the differences:
## For refseq-pseudofinder vs refseq-refseq then prokka-pseudofinder vs refseq-refseq
echo "refseq-pseudofinder,3
prokka-pseudofinder,4" | \
sed 's/,/\t/g' | \
while read comparison target_col; do
echo "${comparison} percent change in intact genes vs. refseq-refseq"
cat $outfile | \
tail -n +2 | \
awk -v col=$target_col '{print 100*($col-$2)/$2}' | \
sort -g \
> sorted.tmp
echo "min = "
head -n 1 sorted.tmp
echo "max = "
tail -n 1 sorted.tmp
echo "median = "
head -n 21 sorted.tmp | tail -n 1
echo "mean = "
cat sorted.tmp | awk 'BEGIN {c=0}; {c+=$0}; END {print c/FNR}'
rm sorted.tmp
done
