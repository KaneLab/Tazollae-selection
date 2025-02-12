## Run early pipeline again using the prokka ORF predictions with both sets of T. az. MAGS
## One set is from Dijkhuizen et al, 2018, New Phytologist / Dijkhuizen et al, 2021, frontiers in Plant Science
## and the other set is from Li et al, 2018, Nature Plants
## All of the strains used by Li et al were also used by Dijkhuizen et al, but the DNA extraction was not from the exact same specimen
## So by comparing the two, we hope to get a sense of the variability due to intra-strain biological variability and extraction/assembly pipelines
## Comments are more minimal in this doc, so see the main prokka/pseudofinder/orthofinder pipelines


############################
#### Li vs. Dijkhuizen #####
########## MAGs ############
############################

## Li MAGs were reassembled by me, using a protocol as close to that described by Dijkhuizen et al as I could get
## Those scripts are in the assembly pipeline folder
## Here, the Li re-assemblies will be run through prokka and pseudofinder
## and then added into the orthofinder outputs.
## To see what the biological within-strain variance is and the effects of DNA extraction, sequencing, and assembly

## Collect and rename the reassemblies 
cd /home/liam/data/fastas
cp -r /home/liam/data/assembly_workflow/metabat2_out/final_MAGs/second_assemblies/ ./Li_reassemblies
cd Li_reassemblies
for file in *.fa; do mv $file ${file%%_final_MAG.fa}.fna; done

## Prokka
cd /home/liam/sandbox
conda deactivate
conda activate prokka
mkdir /home/liam/data/main_pipeline/Li_vs_Dijkhuizen
outdir=/home/liam/data/main_pipeline/Li_vs_Dijkhuizen/prokka_out
mkdir $outdir
indir=/home/liam/data/fastas/Li_reassemblies
## Loop through the reassemblies
cat /home/liam/Li.prefix.list | \
while read prefix; do
echo $prefix;
## Modify the headers because some are too long for prokka
cat $indir/${prefix}.fna | sed -E "s/^>(\S+).*$/>${prefix}_\1/" > infile.tmp
prokka --cpus 20 --compliant --rfam --kingdom Bacteria --outdir ${outdir}/$prefix --locustag $prefix --prefix $prefix infile.tmp;
rm infile.tmp
## Finish loop through the reassemblies
done


## Pseudofinder
indir=/home/liam/data/main_pipeline/Li_vs_Dijkhuizen/prokka_out
outdir=/home/liam/data/main_pipeline/Li_vs_Dijkhuizen/pseudofinder_out
mkdir $outdir
conda deactivate
conda activate pseudofinder
## Generated db in main pseudofinder pipeline
db=/home/liam/data/main_pipeline/pseudofinder_out/allFreeLivingNostoc.faa.dmnd
## Loop through prokka outputs for reassmblies
cat /home/liam/Li.prefix.list | \
while read query; do
mkdir ${outdir}/${query}
cd ${outdir}/${query} 
in_gbk=${indir}/${query}/${query}.gbk
## -ce flag adds contig ends as potential intergenic regions.
## Using this flag to get the most possible hits.
## Use use_deviation flag. Compared with and without in main pseudofinder pipeline.
python ~/tools/pseudofinder/pseudofinder.py annotate -t 20 -ce --use_deviation -l 0.65 --diamond -skpdb -hc 46 --genome $in_gbk --outprefix ${query} --database $db;
## End query loop for pseudofinder run wrapper
done
## Change the headers in ${prefix}_pseudos.fasta to be shorter and easier to read
## Save in new file ${prefix}_pseudos.shorter_headers.fasta
cd $outdir
cat ~/Li.prefix.list | \
while read prefix; do
cat ${prefix}/${prefix}_pseudos.fasta | sed 's/_[^_]\+_pseudo/_p/' > ${prefix}/${prefix}_pseudos.shorter_headers.fasta
done


## orthofinder
conda deactivate
conda activate orthofinder
outdir=/home/liam/data/main_pipeline/Li_vs_Dijkhuizen/orthofinder_out
mkdir $outdir
indir=${outdir}/input_fastas
outdir=${outdir}/initial_run
## Get the inputs to the main pipeline orthofinder
cp -r /home/liam/data/main_pipeline/orthofinder_out/input_fastas $indir
## Get the Li pseudofinder outputs
cat ~/Li.prefix.list | \
while read prefix; do
cp /home/liam/data/main_pipeline/Li_vs_Dijkhuizen/pseudofinder_out/${prefix}/${prefix}_intact.faa ${indir}/${prefix}.fasta
done
## Run orthofinder!
orthofinder -y -M msa -o $outdir -f $indir

## Reroot species tree with gloeobacter as outgroup
conda deactivate
conda activate gotree
initial_orthofinder_outdir=/home/liam/data/main_pipeline/Li_vs_Dijkhuizen/orthofinder_out/initial_run/Results_Aug21
cd ${initial_orthofinder_outdir}/Species_Tree/
gotree reroot outgroup --strict -l ~/outgroup.short_prefix.list -i SpeciesTree_rooted.txt -o SpeciesTree_RErooted.txt
cp SpeciesTree_RErooted.txt ${initial_orthofinder_outdir}/../..

## Re-run orthofinder with re-rooted tree
cd ${initial_orthofinder_outdir}/../..
conda deactivate
conda activate orthofinder
orthofinder -ft $initial_orthofinder_outdir -s SpeciesTree_RErooted.txt 
## Using ls -thor, it doesn't look like anything changed in the original results folder
## Move and rename output folders
cd ${initial_orthofinder_outdir}/../..
mv initial_run/Results_Aug22 final_run
mv initial_run/Results_Aug21/* initial_run
rmdir initial_run/Results_Aug21
rm SpeciesTree_RErooted.txt

## From here down, look at how similar or divergent the results are for the Li MAGs from the Dijkhuizen MAGs
orthofinder_outdir=/home/liam/data/main_pipeline/Li_vs_Dijkhuizen/orthofinder_out/final_run
orthofinder_downstream_dir=/home/liam/data/main_pipeline/Li_vs_Dijkhuizen/orthofinder_out/for_downstream
mkdir $orthofinder_downstream_dir
cd $orthofinder_downstream_dir

#Make a table with counts of how many loci from each genome are present in each N8.HOG
outfile=N8_noOutgroups.tsv
infile=${orthofinder_outdir}/Phylogenetic_Hierarchical_Orthogroups/N8.tsv
##Remove columns corresponding to genomes that aren't in N8
##Get the names of the genomes in N8
awk -F "\t" '{print $1"\t"$2"\t"$3}' $infile > $outfile
tail -n +2 $infile | awk -F "\t" '{for (i=4;i<=NF;++i) {if ($i!="") {print i}}}' | sort -un > N8.cols.tmp
cat N8.cols.tmp | \
while read col; do
awk -v col="$col" -F "\t" '{print $col}' $infile > nextCol.tmp
paste $outfile nextCol.tmp > tmp
mv tmp $outfile
done

infile=N8_noOutgroups.tsv
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
cd ${orthofinder_downstream_dir}/../../
## Get the short prefixes
head -n 1 ${orthofinder_downstream_dir}/N8_noOutgroups.tsv | awk -F "\t" ' { for (i=4;i<=NF;++i) {print $i} } ' > N8.short_prefix.list
## Get the long prefixes for everything except Li
grep -Ff N8.short_prefix.list ~/All.shortened.tsv > N8.shortened.tsv
## Make a double column for each Li genome (short same as long) and append
paste ~/Li.prefix.list ~/Li.prefix.list >> N8.shortened.tsv
##Check I got everything
wc -l N8.short_prefix.list
grep -Ff N8.short_prefix.list N8.shortened.tsv | wc -l
## Equal! Success! Can remove the single col.
rm N8.short_prefix.list 

## Temporarily add the required pseudofinder outputs from the main pipeline
## to this pipeline's pseudofinder folder
## in the format of /pseudofinder_out/${long_prefix}/files
## so they can be input to assign_pseudogenes_to_hogs.sh
this_pseudofinder_outdir=/home/liam/data/main_pipeline/Li_vs_Dijkhuizen/pseudofinder_out
main_pseudofinder_outdir=/home/liam/data/main_pipeline/pseudofinder_out
grep -vFf ~/Li.prefix.list N8.shortened.tsv | \
while read long_prefix short_prefix; do
mkdir ${this_pseudofinder_outdir}/${long_prefix}
cp ${main_pseudofinder_outdir}/${long_prefix}/${long_prefix}_pseudos.shorter_headers.fasta \
${this_pseudofinder_outdir}/${long_prefix}
cp ${main_pseudofinder_outdir}/${long_prefix}/${long_prefix}_intact.faa \
${this_pseudofinder_outdir}/${long_prefix}
done

## Inputs:
outdir=/home/liam/data/main_pipeline/Li_vs_Dijkhuizen/assign_pseudogenes
N=N8
prefix_tsv=/home/liam/data/main_pipeline/Li_vs_Dijkhuizen/N8.shortened.tsv
pseudofinder_outdir=/home/liam/data/main_pipeline/Li_vs_Dijkhuizen/pseudofinder_out
orthofinder_downstream_dir=/home/liam/data/main_pipeline/Li_vs_Dijkhuizen/orthofinder_out/for_downstream

## Call assign_pseudogenes_to_hogs.sh !
bash -i /home/liam/scripts/helper/assign_pseudogenes_to_hogs.sh \
-p $pseudofinder_outdir \
-f $orthofinder_downstream_dir \
-n $N \
-l $prefix_tsv \
-o $outdir

## Remove redundant pseudofinder outputs
grep -vFf ~/Li.prefix.list N8.shortened.tsv | \
while read long_prefix short_prefix; do
rm -r ${this_pseudofinder_outdir}/${long_prefix}
done

##############
## Analyses ##
##############

## A couple things to do before the loop below:
mkdir /home/liam/data/main_pipeline/Li_vs_Dijkhuizen/orthofinder_out/for_downstream/aminoAcid_HOGs
mkdir /home/liam/data/main_pipeline/Li_vs_Dijkhuizen/orthofinder_out/for_downstream/nucleotide_HOGs
## Check that there is a prefix on every header of the amino acid fastas.
cd /home/liam/data/main_pipeline/Li_vs_Dijkhuizen/orthofinder_out/input_fastas
> tmp
for file in *.fasta; do
prefix=${file%.fasta};
echo $prefix >> tmp;
if (( $(grep -c "^>" $file) == $(grep "^>" $file | grep -c "$prefix") && $(grep -v "^>" $file | grep -c "$prefix") == 0 )); then
echo "looks good!" >> tmp
else
echo "looks bad!" >> tmp
fi 
done
echo "How many look bad (hopefully 0):"
grep -B 1 "looks bad!" tmp
echo "How many look good:"
grep -c "looks good!" tmp
echo "This should be 2x look good number:"
wc -l < tmp
rm tmp
## Looks good!
## Make a 2-col prefix list of Dijk-Li for corresponding genome pairs
## The sort does put the carolinianas with their appropriate partners
grep -v "FiliculoidesGal" ~/Dijkhuizen.prefix.list | sort > ~/Dijk_Li.prefix.list
sort ~/Li.prefix.list | sort > tmp
paste ~/Dijk_Li.prefix.list tmp > tmp.tmp
mv tmp.tmp ~/Dijk_Li.prefix.list
rm tmp


## Loop through the Dijkhuizen-Li genome pairs (6)
Li_Dijk_analysis_dir=/home/liam/data/main_pipeline/Li_vs_Dijkhuizen/analyses
mkdir $Li_Dijk_analysis_dir
orthofinder_downstream_dir=/home/liam/data/main_pipeline/Li_vs_Dijkhuizen/orthofinder_out/for_downstream
cat ~/Dijk_Li.prefix.list | \
while read dijk li; do
cd ${Li_Dijk_analysis_dir}
mkdir $dijk
cd $dijk
log=${Li_Dijk_analysis_dir}/$dijk/comparison.log.txt
> $log

## First section compares basic genome stats
pseudofinder_dijk_dir=/home/liam/data/main_pipeline/pseudofinder_out/$dijk
pseudofinder_li_dir=/home/liam/data/main_pipeline/Li_vs_Dijkhuizen/pseudofinder_out/$li
echo "~~ First, let's compare some basic stats between the two genomes. ~~" >> $log
echo "~~ All values in this section in the format Dijk,Li. ~~" >> $log
echo "" >> $log
echo "~~ What are the genome names ~~" >> $log
echo "${dijk},${li}" >> $log
echo "" >> $log
dijk_genome=/home/liam/data/fastas/Dijkhuizen_genomes/${dijk}.fna
li_genome=/home/liam/data/fastas/Li_reassemblies/${li}.fna
echo "~~ How many nts are the genomes? ~~" >> $log
cat $dijk_genome | grep --no-group-separator -v "^>" | tr -d "\n" | awk '{print $0}' > dijk_tmp
cat $li_genome | grep --no-group-separator -v "^>" | tr -d "\n" | awk '{print $0}' > li_tmp
a=$(wc -c < dijk_tmp)
b=$(wc -c < li_tmp)
echo "$a,$b" >> $log
echo "" >> $log
echo "~~ How many contigs are the genomes ~~" >> $log
a=$(grep -c "^>" $dijk_genome)
b=$(grep -c "^>" $li_genome)
echo "$a,$b" >> $log
echo "" >> $log
echo "~~ What is the GC percent of the genomes ~~" >> $log
l_noN=$(cat dijk_tmp | tr -d "nN" | wc -c)
gc=$(cat dijk_tmp | tr -cd "GCgc" | wc -c)
a=$(( 100*$gc/$l_noN ))
l_noN=$(cat li_tmp | tr -d "nN" | wc -c)
gc=$(cat li_tmp | tr -cd "GCgc" | wc -c)
b=$(( 100*$gc/$l_noN ))
echo "$a,$b" >> $log
echo "" >> $log
echo "~~ How many total genes are predicted after prokka and pseudofinder ~~" >> $log
a=$(( $(grep -c "^>" ${pseudofinder_dijk_dir}/${dijk}_intact.faa) + \
$(grep -c "^>" ${pseudofinder_dijk_dir}/${dijk}_pseudos.fasta) ))
b=$(( $(grep -c "^>" ${pseudofinder_li_dir}/${li}_intact.faa) + \
$(grep -c "^>" ${pseudofinder_li_dir}/${li}_pseudos.fasta) ))
echo "$a,$b" >> $log
echo "" >> $log
echo "~~ How many of those are predicted intact (remainder are predicted pseudogenes)? ~~" >> $log
a=$(grep -c "^>" ${pseudofinder_dijk_dir}/${dijk}_intact.faa)
b=$(grep -c "^>" ${pseudofinder_li_dir}/${li}_intact.faa)
echo "$a,$b" >> $log
echo "" >> $log

## Second section
## Here, I am checking how similar the gene sets between the two are
## as determined by orthofinder-derived N8.HOGs.
## Get the appropriate columns from the gene count N8.HOG x genome table
## And narrow it down to only the rows that aren't 0's
cat ${orthofinder_downstream_dir}/n_intact.tsv | \
awk -v dijk=$dijk -v li=$li -F "\t" '
NR==1 {for (i=1; i<=NF; ++i) {
if ($i=="HOG") {A[1]=i}
else if ($i==dijk) {A[2]=i}
else if ($i==li) {A[3]=i} } };
{print $A[1]"\t"$A[2]"\t"$A[3]}
' | \
awk -F "\t" '$2 != 0 || $3 != 0 {print $0}' \
> comparison_cols.tsv
## Get the differences
cat comparison_cols.tsv | \
awk -F "\t" '$2 != $3 {print $0}' > diffs.tsv

echo "~~ How similar are the orthofinder HOG profiles of these two genomes? ~~" >> $log
echo "" >> $log
echo "~~ How many HOGs have at least 1 locus in both of these genomes? ~~" >> $log
echo $( tail -n +2 comparison_cols.tsv | awk ' $2 != 0 && $3 != 0 ' | wc -l ) >> $log
echo "" >> $log
echo "~~ Of those, how many have a different number from the two genomes? ~~" >> $log
echo $( tail -n +2 diffs.tsv | awk ' $2 != 0 && $3 != 0 ' | wc -l ) >> $log
echo "" >> $log
echo "~~ How many HOGs are present in only one of the genomes (so not included in the 2 preceding counts) (Dijkhuizen,Li)? ~~" >> $log
a=$( tail -n +2 diffs.tsv | awk '$3==0' | wc -l )
b=$( tail -n +2 diffs.tsv | awk '$2==0' | wc -l )
echo "$a,$b" >> $log
echo "" >> $log
echo "~~ Double check the above sum in a way that makes sense ~~" >> $log
echo "" >> $log

## Third section: Curious the effects of differential pseudogene calls and high copy number hogs
## Look at the presence/absence N8.HOGs in the table with intact and pseudogene counts.
echo "~~ How much of an effect do differential pseudogene calls and high copy number HOGs have on the differences above? ~~" >> $log
echo "" >> $log
cd ${Li_Dijk_analysis_dir}/$dijk
cat ${orthofinder_downstream_dir}/n_intact_and_pseudo.tsv | \
awk -v dijk=$dijk -v li=$li -F "\t" '
NR==1 {for (i=1; i<=NF; ++i) {
if ($i=="HOG") {A[1]=i}
else if ($i==dijk) {A[2]=i}
else if ($i==li) {A[3]=i} } };
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
echo "~~ Of HOGs with predicted-intact loci from only one genome, how many have more pseudogenes in the other genome (intact in Dijk, intact in Li)? ~~" >> $log
a=$(cat presAbsDiff_wPseudos.tsv | awk -F "\t" '$5 > $3 && $4 == 0 {print $0}' | wc -l)
b=$(cat presAbsDiff_wPseudos.tsv | awk -F "\t" '$3 > $5 && $2 == 0 {print $0}' | wc -l)
echo "$a,$b" >> $log
echo "" >> $log 
echo "~~ So differential pseudogene calls likely account for this portion of presence/absence discrepancies for predicted intact (intact in Dijk, intact in Li) ~~" >> $log
c=$(tail -n +2 diffs.tsv | awk '$3 == 0' | wc -l)
d=$(tail -n +2 diffs.tsv | awk '$2 == 0' | wc -l)
echo "$((100*$a/$c)),$((100*$b/$d))" >> $log
echo "" >> $log
echo "~~ Of the HOGs that only had predicted-intact loci from one genome, how many had >4 loci from that genome (count (percent))? ~~" >> $log
tail -n +2 diffs.tsv | awk ' $2 == 0 || $3 == 0 ' > tmp
t=$( wc -l < tmp )
c=$( cat tmp | awk -F "\t" ' $2 > 4 || $3 > 4 ' | wc -l ) 
echo "$c ($((100*$c/$t)))" >> $log
echo "" >> $log
echo "~~ Of the HOGs that had predicted-intact loci from both genomes, but differed in number, how many had >4 loci in at least one genome (count (percent))? ~~" >> $log
tail -n +2 diffs.tsv | awk ' $2 != 0 && $3 != 0 ' > tmp
t=$( wc -l < tmp )
c=$( cat tmp | awk ' $2 > 4 || $3 > 4 ' | wc -l )
echo "$c ($((100*$c/$t)))" >> $log
echo "" >> $log

## Last section
## Compare the nucleotide and amino acid identities of seqs in HOGs where Dijk and Li each have 1 locus
echo "~~ Compare the nucleotide and amino acid seqs for HOGs with exactly 1 predicted-intact locus present from each genoem ~~" >> $log
echo "" >> $log
tail -n +2 comparison_cols.tsv | \
awk -F "\t" ' $2==1 && $3==1 {print $1} ' \
> dijk_li_single_copy.list
echo "~~ How many HOGs have 1 copy in each? ~~" >> $log
wc -l < dijk_li_single_copy.list >> $log

## loop through the hogs that have a single copy from each of li and dijk and grab the amino acid seqs
cd /home/liam/data/main_pipeline/Li_vs_Dijkhuizen/orthofinder_out/for_downstream/aminoAcid_HOGs
mkdir $dijk
cd $dijk
indir=/home/liam/data/main_pipeline/Li_vs_Dijkhuizen/orthofinder_out/input_fastas
cat ${Li_Dijk_analysis_dir}/${dijk}/dijk_li_single_copy.list | \
while read hog; do
> ${hog}.faa
grep "^${hog}" /home/liam/data/main_pipeline/Li_vs_Dijkhuizen/orthofinder_out/for_downstream/N8_noOutgroups.tsv | \
awk -F "\t" '{ for (i=4;i<=NF;++i) { print $i } }' | \
awk '$0!=""' \
> this_hog_locus.list
this_dijk_locus=$(grep $dijk this_hog_locus.list)
this_li_locus=$(grep $li this_hog_locus.list)
rm this_hog_locus.list
grep -A 1 $this_dijk_locus ${indir}/${dijk}.fasta >> ${hog}.faa
grep -A 1 $this_li_locus ${indir}/${li}.fasta >> ${hog}.faa
done

## loop through the hogs that have a single copy from each of li and dijk and grab the nucleotide seqs
cd /home/liam/data/main_pipeline/Li_vs_Dijkhuizen/orthofinder_out/for_downstream/nucleotide_HOGs
mkdir $dijk
cd $dijk
dijk_indir=/home/liam/data/main_pipeline/pseudofinder_out/$dijk
li_indir=/home/liam/data/main_pipeline/Li_vs_Dijkhuizen/pseudofinder_out/$li
cat ${Li_Dijk_analysis_dir}/${dijk}/dijk_li_single_copy.list | \
while read hog; do
> ${hog}.ffn
grep "^${hog}" /home/liam/data/main_pipeline/Li_vs_Dijkhuizen/orthofinder_out/for_downstream/N8_noOutgroups.tsv | \
awk -F "\t" '{ for (i=4;i<=NF;++i) { print $i } }' | \
awk '$0!=""' \
> this_hog_locus.list
this_dijk_locus=$(grep $dijk this_hog_locus.list)
this_li_locus=$(grep $li this_hog_locus.list)
rm this_hog_locus.list
grep -A 1 $this_dijk_locus ${dijk_indir}/${dijk}_intact.ffn >> ${hog}.ffn
grep -A 1 $this_li_locus ${li_indir}/${li}_intact.ffn >> ${hog}.ffn
done

## Now get a list of which sequences are not identical using cmp on nucleotide!
cd ${Li_Dijk_analysis_dir}/$dijk
outfile=single_copy_Seqs_not_identical_nuc.txt
> $outfile
indir=/home/liam/data/main_pipeline/Li_vs_Dijkhuizen/orthofinder_out/for_downstream/nucleotide_HOGs/$dijk
cat dijk_li_single_copy.list | \
while read hog; do
grep -A 1 $dijk ${indir}/${hog}.ffn | tail -n 1 > dijk_seq.tmp
grep -A 1 $li ${indir}/${hog}.ffn | tail -n 1 > li_seq.tmp
if ! cmp -s dijk_seq.tmp li_seq.tmp; then
echo $hog >> $outfile
fi
done
rm *.tmp
echo "~~ How many of the single copy HOGs differ in nucleotide sequence? ~~" >> $log
wc -l < $outfile >> $log
echo "" >> $log

## Now get a list of which sequences are not identical using cmp on protein!
cd ${Li_Dijk_analysis_dir}/$dijk
outfile=single_copy_Seqs_not_identical_aa.txt
> $outfile
indir=/home/liam/data/main_pipeline/Li_vs_Dijkhuizen/orthofinder_out/for_downstream/aminoAcid_HOGs/$dijk
cat dijk_li_single_copy.list | \
while read hog; do
grep -A 1 $dijk ${indir}/${hog}.faa | tail -n 1 > dijk_seq.tmp
grep -A 1 $li ${indir}/${hog}.faa | tail -n 1 > li_seq.tmp
if ! cmp -s dijk_seq.tmp li_seq.tmp; then
echo $hog >> $outfile
fi
done
rm *.tmp
echo "~~ How many of the single copy HOGs differ in amino sequence? ~~" >> $log
wc -l < $outfile >> $log

rm *tmp
## Finish loop through Dijk-Li genome pairs
done

## combine the output files into one.
cd /home/liam/data/main_pipeline/Li_vs_Dijkhuizen/analyses
> comparison_table.tsv
cat ~/Dijk_Li.prefix.list | \
while read dijk li; do
paste comparison_table.tsv ${dijk}/comparison.log.txt > tmp
mv tmp comparison_table.tsv
done
sed -i 's/^\t\+//' comparison_table.tsv
sed -i 's/^\(~~[^~]*~~\).*/\1/' comparison_table.tsv
