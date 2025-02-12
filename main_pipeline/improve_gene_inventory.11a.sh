## Some KO's are only present in the eggnog-mapper outputs for NosAzo0708
## and not for any of the MAGs.
## Some of these are genes that are also present in the chloroplast,
## and thus the reads may have been filtered out during assembly.
## For instance
## K02703 (psbA): 3 (intact) copies in NosAzo0708, 1 in Nilotica5001. Also, 1 pseudogene in each MAG and 2 in NosAzo0708 although those 2 are right next to one another. Several of the MAG pseudogenes come very close to the edges of their contigs which suggests that the assembly might be falsely creating these pseudogenes.
## K02706 (psbD): 2 (intact) copies in NosAzo0708
## K02723 (psbY): 1 (intact) copy in all but Carolinianas (0). NOT IN THE CHLOROPLAST
## K02711 (psbJ): 1 (intact) copy in NosAzo0708
## So, we want to go back to the short-reads to try to assemble these gene sequences and determine if they are actually present
## There are two parts of this script.
## First, which is the majority of the script:
## Reassemble the MAGs using all reads that align best in bwa mem to the reference cyanobiont
## including reads that ALSO align to the nuclear or chloroplast host genome
## Look for more intact genes in these re-assemblies
## Second: determine which pseudogenes from the main pipeline are at the ends of contigs
## Defined as being between the first/last ORF and the start/stop of the contig.
## These will have their own classification because it is hard to confidently predict if these are actually pseudogenes.

## Prepare folder for this analysis
mkdir ~/data/chloroplast_genes
mkdir ~/data/chloroplast_genes/Dijkhuizen_SRAs

## Downloaded the raw reads from ENA for Dijkhuizen T. az. assemblies
cd ~/data/chloroplast_genes/Dijkhuizen_SRAs
## Filiculoides: https://www.ebi.ac.uk/ena/browser/view/SAMEA104284504?show=reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR211/007/ERR2114807/ERR2114807_1.fastq.gz &&\
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR211/007/ERR2114807/ERR2114807_2.fastq.gz
## Rubra: https://www.ebi.ac.uk/ena/browser/view/SAMN08375388?show=reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/000/SRR6480160/SRR6480160_1.fastq.gz &&\
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/000/SRR6480160/SRR6480160_2.fastq.gz
## Microphylla: https://www.ebi.ac.uk/ena/browser/view/SAMN08375389?show=reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/001/SRR6480161/SRR6480161_1.fastq.gz &&\
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/001/SRR6480161/SRR6480161_2.fastq.gz
## Mexicana: https://www.ebi.ac.uk/ena/browser/view/SAMN08375390?show=reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/009/SRR6480159/SRR6480159_1.fastq.gz &&\
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/009/SRR6480159/SRR6480159_2.fastq.gz
## Nilotica: https://www.ebi.ac.uk/ena/browser/view/SAMN08375391?show=reads
## Nilotica has 2 sets of paired read files. Looking at the assembled genome (https://www.ebi.ac.uk/ena/browser/view/CAJSYP010000000?show=sample-attributes)
## the SRA is the same as for both of those paired read sets. 
## So I will use the top one on the page, arbitrarily.
## Dijkhuizen et al, 2018 and 2021 did not mention anything relevant. And the sample accession + tax ID are each identical for both.
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/006/SRR6480196/SRR6480196_1.fastq.gz &&\
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/006/SRR6480196/SRR6480196_2.fastq.gz
## Caroliniana3017: https://www.ebi.ac.uk/ena/browser/view/SAMN08375392?show=reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/001/SRR6480201/SRR6480201_1.fastq.gz &&\
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/001/SRR6480201/SRR6480201_2.fastq.gz
## Caroliniana3004: https://www.ebi.ac.uk/ena/browser/view/SAMN08375393?show=reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/001/SRR6480231/SRR6480231_1.fastq.gz &&\
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/001/SRR6480231/SRR6480231_2.fastq.gz
## Gunzip and rename
## I think gunzip is actually unnecessary and trimmomatic can just take the .gz files

## Get a fasta of all of the reference Azolla filiculoides chloroplast genes
## and all of the reference Nostoc Azollae 0708 genes that are annotated with the same gene names as those genes
## There is a v1.1 and a v1.0 on fernbase.org for this genome,
## but using cmp, the .fasta and .gff for the chloroplast are identical between the two versions
mkdir ~/data/chloroplast_genes/chloroplast_genes
cd ~/data/chloroplast_genes/chloroplast_genes
wget https://fernbase.org/ftp/Azolla_filiculoides/Azolla_asm_v1.1/chloroplast_genome/Azolla_filiculoides.cp_genome_v1_4.fasta &&\
wget https://fernbase.org/ftp/Azolla_filiculoides/Azolla_asm_v1.1/chloroplast_genome/Azolla_filiculoides.cp_genome_v1_4.gff
mv Azolla_filiculoides.cp_genome_v1_4.fasta chloroplast.fasta
mv Azolla_filiculoides.cp_genome_v1_4.gff chloroplast.gff
## Get CDS's, not genes. Genes can have introns in the chloroplasts, but not in the cyanos.
## Also, genes appear to include tRNAs, while CDS's do not.
## sort -u because one CDS is repeated identically
cat chloroplast.gff | awk '$3=="CDS" {print $4"\t"$5"\t"$7"\t"$9}' | sed 's/Name=//' | sort -k1,1n -u > cds.tsv
## For genes with more than locus, name them 

## Get the nucleotide sequences for the CDS's
> chloroplast.ffn
i=0
cat cds.tsv | \
while read start stop dir name; do
i=$(($i+1))
echo ">chloroplast_${i}_${name}" >> chloroplast.ffn
if [[ "$dir" == "+" ]]; then
	tail -n +2 chloroplast.fasta | \
	awk -v start=$start -v stop=$stop '{print substr($0,start,(stop-start+1))}' \
	>> chloroplast.ffn
else
	tail -n +2 chloroplast.fasta | \
	awk -v start=$start -v stop=$stop '{print substr($0,start,(stop-start+1))}' | \
	tr ACGTacgt TGCAtgca | rev \
	>> chloroplast.ffn
fi
done
## NOT ALL OF THE ORFs have lengths DIVISIBLE BY 3! FOR INSTANCE, TWO LOCI NEAR EACH OTHER LABELED "petD" that I assume are exons...

## Get the gene names of those CDSs. Every one has a name.
cat cds.tsv | sed 's/^.*Name=//' | awk '{print $1}' | sort -u > cds_gene_names.list
## 85 unique gene names out of 106 CDSs
##

## Trim reads
## First, figure out the adapter to use by running fastqc on the reads.
cd ~/data/chloroplast_genes/Dijkhuizen_SRAs
conda activate fastqc
cat ~/Dijkhuizen.prefix.list | while read prefix; do fastqc -o ./ ${prefix}_*.fastq.gz; done
rm *.zip
## scp to local and open in web browser
## "Adapter Content" is empty for all runs of fastqc
## run trimmomatic
conda deactivate
conda activate trimmomatic
adapterFile=~/miniconda3/envs/trimmomatic/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa
cat ~/Dijkhuizen.prefix.list | while read prefix; do
in_1=${prefix}_1.fastq.gz
in_2=${prefix}_2.fastq.gz
outputBase=${prefix}_trimmed.fq.gz
trimmomatic PE -threads 20 $in_1 $in_2 -baseout $outputBase ILLUMINACLIP:${adapterFile}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 && \
done

## Align reads to the reference cyano, nuclear, and chloroplast genomes.
mkdir ~/data/chloroplast_genes/bwa_mem_out
cd ~/data/chloroplast_genes/bwa_mem_out
conda deactivate
conda activate samtools
## Put the 3 reference genomes together in one fasta
cat ~/data/refseq_Nostocales_genomes/data/NostocAzollae0708/NostocAzollae0708.fna \
~/data/reference/plantRefGenome_inclChlp.fasta \
> reference.fasta
## simplify the headers
cat reference.fasta | awk '{print $1}' > tmp
mv tmp reference.fasta

indir=~/data/chloroplast_genes/Dijkhuizen_SRAs
reference=reference.fasta
bwa index $reference
## Loop through MAGs, aligning to the reference
cat ~/Dijkhuizen.prefix.list | \
while read prefix; do
## Only use paired reads, which is the large majority.
## Makes code simpler and might make it more confident to which reference genome the paired read has aligned.
bwa mem -t 20 $reference ${indir}/${prefix}_trimmed_1P.fq.gz ${indir}/${prefix}_trimmed_2P.fq.gz | gzip > ${prefix}_SRA_combined_reference.sam.gz
done

## Get read pairs where at least one partner aligns best to the cyanobiont genome
## and neither partner aligns best to the nuclear or chloroplast host genomes
## For most reads, there are two lines in the SAM file.
## One for each paired forward/backward.
## But for "chimeric" reads, there could be more than 2.
## The T. az. chromosome and plasmids are all labeled with NC_0142*
cat ~/Dijkhuizen.prefix.list | \
while read prefix; do
echo $prefix
zcat ${prefix}_SRA_combined_reference.sam | \
grep -v "^@" | \
## I don't think sort is necessary but I am confident it can't hurt
## The sort is to make sure that paired reads are consecutive
sort -k1,1 | \
awk '
BEGIN {prev=""; p=0};
{
if ( $1 != prev ) {
	if (p>0) { print prev };
	prev=$1;
	p=0;
}
if ($3 ~ "^NC_0142") {p=p+1}
else if ($3 ~ "^Azfi" || $3 ~ "^Azolla_cp") {p=p-1};
}
' \
> ${prefix}_reads_for_assembly.list
done
## Add an @ symbol and " " to the filtered read names to make sure hits are exact when searching against all the reads.
cat ~/Dijkhuizen.prefix.list | \
while read prefix; do
cat ${prefix}_reads_for_assembly.list | \
sed 's/^/@/' | sed 's/$/ /' \
> tmp
mv tmp ${prefix}_reads_for_assembly.list
done


## Get the filtered reads!
mkdir ~/data/chloroplast_genes/filtered_reads
cd ~/data/chloroplast_genes/filtered_reads
conda deactivate
sam_dir=~/data/chloroplast_genes/bwa_mem_out
reads_dir=~/data/chloroplast_genes/Dijkhuizen_SRAs
cat ~/Dijkhuizen.prefix.list | \
while read prefix; do
zgrep --no-group-separator -A3 -Ff ${sam_dir}/${prefix}_reads_for_assembly.list ${reads_dir}/${prefix}_trimmed_1P.fq.gz | gzip > ${prefix}_filtered_1P.fq.gz && \
zgrep --no-group-separator -A3 -Ff ${sam_dir}/${prefix}_reads_for_assembly.list ${reads_dir}/${prefix}_trimmed_2P.fq.gz | gzip > ${prefix}_filtered_2P.fq.gz
done

## ASSEMBLE!
indir=~/data/chloroplast_genes/filtered_reads
mkdir ~/data/chloroplast_genes/assemblies
cd ~/data/chloroplast_genes/assemblies
conda deactivate
conda activate spades
cat ~/Dijkhuizen.prefix.list | \
while read prefix; do
~/miniconda3/envs/spades/bin/spades.py --meta -o $prefix -1 ${indir}/${prefix}_filtered_1P.fq.gz -2 ${indir}/${prefix}_filtered_2P.fq.gz
done
conda deactivate

## First look at the assemblies:
cd ~/data/chloroplast_genes/assemblies
echo "genome" > prefix.tmp
echo "n_contigs" > n_contigs.tmp
echo "n_nucleotides" > total_length.tmp
cat ~/Dijkhuizen.prefix.list | \
while read prefix; do
echo $prefix >> prefix.tmp
grep -c "^>" ${prefix}/scaffolds.fasta >> n_contigs.tmp
grep "^>" ${prefix}/scaffolds.fasta | sed 's/.*length_//' | sed 's/_.*//' | awk 'BEGIN {c=0}; {c+=$0}; END {print c}' >> total_length.tmp
done
paste prefix.tmp n_contigs.tmp total_length.tmp > assembly_stats_all_contigs.tsv
rm prefix.tmp n_contigs.tmp total_length.tmp
## Dijkhuizen assemblies have a contig length cut-off of 2500
echo "genome" > prefix.tmp
echo "n_contigs" > n_contigs.tmp
echo "n_nucleotides" > total_length.tmp
cat ~/Dijkhuizen.prefix.list | \
while read prefix; do
echo $prefix >> prefix.tmp;
grep "^>" ${prefix}/scaffolds.fasta | awk -F "_" '$4>=2500 {print $4}' > lengths.tmp;
wc -l < lengths.tmp >> n_contigs.tmp;
cat lengths.tmp | awk 'BEGIN {c=0}; {c+=$0}; END {print c}' >> total_length.tmp;
done
paste prefix.tmp n_contigs.tmp total_length.tmp > assembly_stats_gte2.5k_contigs.tsv
rm prefix.tmp n_contigs.tmp total_length.tmp lengths.tmp

## Prokka
conda deactivate
conda activate prokka
outdir=~/data/chloroplast_genes/prokka_out
mkdir $outdir
cd $outdir
cat ~/Dijkhuizen.prefix.list | \
while read prefix; do
## Change the headers to contain the prefix and the spades node number
cat ~/data/chloroplast_genes/assemblies/${prefix}/scaffolds.fasta | \
sed 's/\(^>NODE_[0-9]\+\)_.*/\1/' | sed "s/^>NODE/>${prefix}/" > infile.tmp
prokka --cpus 20 --compliant --rfam --kingdom Bacteria --outdir ${outdir}/$prefix --locustag $prefix --prefix $prefix infile.tmp;
rm infile.tmp
done

## Pseudofinder
conda deactivate
conda activate pseudofinder
indir=~/data/chloroplast_genes/prokka_out
outdir=~/data/chloroplast_genes/pseudofinder_out
mkdir $outdir
cd $outdir
db=~/data/main_pipeline/pseudofinder_out/allFreeLivingNostoc.faa.dmnd
cat ~/Dijkhuizen.prefix.list | \
while read prefix; do
mkdir $prefix
cd $prefix
in_gbk=${indir}/${prefix}/${prefix}.gbk
python ~/tools/pseudofinder/pseudofinder.py annotate -t 20 -ce --use_deviation -l 0.65 --diamond -skpdb -hc 46 --genome $in_gbk --outprefix $prefix --database $db;
cd ..
done
##Change the headers in pseudofinder_out/${prefix}/${prefix}_pseudos.fasta
##To be shorter and easier to read
##Save in new file pseudofinder_out/${prefix}/${prefix}_pseudos.shorter_headers.fasta
cd $outdir
cat ~/Taz.shortened.tsv | \
grep -v "NosAzo0708" | \
while read long_prefix short_prefix; do
cat ${long_prefix}/${long_prefix}_pseudos.fasta | sed 's/_[^_]\+_pseudo/_p/' > ${long_prefix}/${long_prefix}_pseudos.shorter_headers.fasta
done

## Orthofinder
conda deactivate
conda activate orthofinder
outdir=~/data/chloroplast_genes/orthofinder_out
mkdir $outdir
cd $outdir
outdir=${outdir}/initial_run
indir=input_fastas
mkdir $indir
## Populate the indir
pseudofinder_outdir=~/data/main_pipeline/pseudofinder_out
prokka_outdir=~/data/main_pipeline/prokka_out
## Get the free-living Nostocales fastas and rename
cat ~/Nostocales_incl0708.shortened.tsv | \
while read long_prefix short_prefix; do
cp ${pseudofinder_outdir}/${long_prefix}/${long_prefix}_intact.faa ${indir}/${short_prefix}.fasta
done
## Get the outgroup fastas and rename (these are prokka outputs b/c they did not pass through pseudofinder)
cat ~/outgroup.shortened.tsv | \
while read long_prefix short_prefix; do
cp ${prokka_outdir}/${long_prefix}/${short_prefix}.faa ${indir}/${short_prefix}.fasta
done
## Get the MAGs which do not have to be renamed and come from a separate pipeline
pseudofinder_outdir=~/data/chloroplast_genes/pseudofinder_out
prokka_outdir=~/data/chloroplast_genes/prokka_out
cat ~/Dijkhuizen.prefix.list | \
while read prefix; do
cp ${pseudofinder_outdir}/${prefix}/${prefix}_intact.faa ${indir}/${prefix}.fasta
done
#Run orthofinder!
#-z and -T raxml-ng flags might slightly improve performance, but at cost of time
##May need to increase soft limit on files open simultneously, but orthofinder error message is very clear for this.
orthofinder -y -M msa -o $outdir -f $indir
#Reroot species tree with gloeobacter as outgroup
conda deactivate
conda activate gotree
initial_orthofinder_outdir=~/data/chloroplast_genes/orthofinder_out/initial_run/Results_Dec23
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
cd ~/data/chloroplast_genes/orthofinder_out
mv initial_run/Results_Dec24 final_run
mv initial_run/Results_Dec23/* initial_run
rmdir initial_run/Results_Dec23

## Orthofinder post-processing
orthofinder_outdir=~/data/chloroplast_genes/orthofinder_out/final_run
cd ${orthofinder_outdir}
mkdir N9_HOGs_locus_lists
cd N9_HOGs_locus_lists
infile=${orthofinder_outdir}/Phylogenetic_Hierarchical_Orthogroups/N9.tsv
head -n 1 $infile | awk -F "\t" '{for (i=4;i<=NF;++i) {print $i}}' > genome_labels.tmp
tail -n +2 $infile | awk -F "\t" '{print $1"\t"$2}' | sort -u > HOGs2OGs.tsv
cat HOGs2OGs.tsv | \
while read hog og; do
grep $hog $infile | awk -F "\t" '{for (i=4;i<=NF;++i) {print $i}}' > this_HOG.tmp
paste genome_labels.tmp this_HOG.tmp > ${hog}_loci.tsv
done
rm *.tmp
## Move some stuff around
mkdir ~/data/chloroplast_genes/orthofinder_out/for_downstream
cd ~/data/chloroplast_genes/orthofinder_out/for_downstream
mv ${orthofinder_outdir}/N9_HOGs_locus_lists ./locus_lists_HOGs
## Get a table of the loci in each N9.HOG
outfile=N9_noOutgroups.tsv
infile=${orthofinder_outdir}/Phylogenetic_Hierarchical_Orthogroups/N9.tsv
##Remove columns corresponding to genomes that aren't in N9
##Get the names of the genomes in N9
awk -F "\t" '{print $1"\t"$2"\t"$3}' $infile > $outfile
tail -n +2 $infile | awk -F "\t" '{for (i=4;i<=NF;++i) {if ($i!="") {print i}}}' | sort -un > N9.cols.tmp
cat N9.cols.tmp | \
while read col; do
awk -v col="$col" -F "\t" '{print $col}' $infile > nextCol.tmp
paste $outfile nextCol.tmp > tmp
mv tmp $outfile
done
rm *.tmp


##############################
###### Compare HOGs from #####
### this and main pipeline ###
##############################

## Determine which HOGs are identical, excluding Dijkhuizen, in this pipeline and in the main pipeline. HOG names may vary.
## First, check that columns of N9_noOutgroups.tsv are in the same order in both pipelines
cd ~/data/chloroplast_genes/orthofinder_out/for_downstream
main_pipeline_N9_noOutgroups=~/data/main_pipeline/orthofinder_out/for_downstream/N9_noOutgroups.tsv
head -n 1 N9_noOutgroups.tsv > tmp
head -n 1 ${main_pipeline_N9_noOutgroups} > tmp.tmp
cmp tmp tmp.tmp
## Looks good!
rm *tmp
echo "N9_noOutgroups.tsv,this_pipeline
${main_pipeline_N9_noOutgroups},main_pipeline" | \
sed 's/,/\t/g' | \
while read file pipeline; do

## Get the table without the header or the second and third cols
cat $file | awk -F "\t" '{printf $1; for (i=4;i<=NF;++i) {printf "\t"$i}; print "" }' > in.tmp

## Get this table without the Dijkhuizen cols
cat ~/Dijkhuizen.prefix.list | \
while read prefix; do
cat in.tmp | \
awk -v prefix=$prefix -F "\t" '
NR==1 { for (i=1;i<=NF;++i) { if ($i==prefix) {c=i} else {printf "\t"$i}; }; print ""; next };
{ for (i=1;i<=NF;++i) { if (i!=c) {printf "\t"$i} }; print "" }
' | \
sed 's/^\t//' \
> tmp
mv tmp in.tmp
## Finish while read prefix loop
done
mv in.tmp ${pipeline}_table.tmp
cat ${pipeline}_table.tmp | tail -n +2 | awk -F "\t" '{ for (i=2;i<=NF;++i) {printf "\t"$i}; print "" }' | sed 's/^\t//'  >> combined_table.tmp

## Finish while read file pipeline loop
done



cat combined_table.tmp | grep "[[:alnum:]]" | sort | uniq -c | awk '$1==2' | sed 's/^      2 //' > matches.tmp
#awk -F "\t" ' { for (i=2;i<=NF;++i) {printf $i"\t"}; print "" } ' > matches.tmp
## How many are identical? I think -1 for each because of headers
## numbers are second attempt (looser read filter) with --meta and then with --isolate flags
## first attempt second attempt (second included reads where only one paired partner mapped to T. az. ref and the other mapped to nothing)
wc -l matches.tmp
## 17274 (94%) 17356 (94%)
## Out of?
wc -l N9_noOutgroups.tsv
## 18465-1 (-1 for header) 18470-1
wc -l ${main_pipeline_N9_noOutgroups}
## 18463-1 (-1 for header)
## Get the names of the hogs for matches
grep -Ff matches.tmp this_pipeline_table.tmp | awk '{print $1}' > matched_hogs_this_pipeline.list
grep -Ff matches.tmp main_pipeline_table.tmp | awk '{print $1}' > matched_hogs_main_pipeline.list
## Look at those that changed
grep -vFf matches.tmp main_pipeline_table.tmp | awk '{print $1}' > unmatched_hogs_main_pipeline.list
grep -Ff unmatched_hogs_main_pipeline.list ~/data/main_pipeline/orthofinder_out/for_downstream/Taz.n_intact.tsv | awk '$2>0' | awk '{ for (i=3;i<=NF;++i) { if ($i==0) {print $0; next} } }' | wc -l
## 835 835, so that is how many HOGs could potentially gain an "intact" from this side pipeline,
## meaning there is at least one intact locus in reference NosAzo0708 and at least one MAG with 0 intact loci
## but where the HOG is not identical between the two pipelines. 835/18462=4.5%
## This does not mean that all of these HOGs actually would gain any loci even if they were matched, just that I don't know if they would have and this is the upper bound.
rm *.tmp

## Instead of the non-MAG orthogroups needing to be identical
## Just use any pair of HOGs from the two pipelines where the NosAzo0708 loci are identical
## For each HOG that contains a NosAzo0708 locus in this pipeline
## make a table with cols the Dijkhuizen MAGs, rows NosAzo0708 loci
## and a 0 (no intact locus) or 1 (1+ intact loci) for each value
## Then swap the NosAzo0708 loci for main pipeline HOGs
## Then add a /0 or /1 to the n_intact_pseudo tables in the main pipeline

## Get the hogs from each pipeline that contain at least 1 NosAzo0708 locus
> hogs_NosAzo0708.tmp
echo "\
~/data/main_pipeline/orthofinder_out/for_downstream,?
./,~\
" | \
sed 's/,/\t/g' | \
while read indir symbol; do
cat ${indir}/N9_noOutgroups.tsv | \
awk -F "\t" '
NR==1 { for (i=1;i<=NF;++i) { if ($i=="NosAzo0708") {c=i} }; next};
$c!="" {print $c"\t"$1}
' | \
sed 's/, /,/g' | \
## Adding the symbol will keep the order of the two pipelines consistent when I use sort
sed "s/$/\t${symbol}/" \
>> hogs_NosAzo0708.tmp
done

## Sort by the NosAzo0708 loci and keep the orthogroups where identical
## If one or both orthogroups contain multiple NosAzo0708 loci and they are not 100% the same, throw out this comparison because interpretation is more confusing.
## Arrange multi-locus hogs locus lists to always have be in ascending order
cat hogs_NosAzo0708.tmp | \
awk -F "\t" '{split($1,A,","); asort(A,B); printf B[1]; for (i=2;i<=length(B);++i) {printf ","B[i]}; print "\t"$2"\t"$3}' \
> tmp
mv tmp hogs_NosAzo0708.tmp
## Check this worked.
cat hogs_NosAzo0708.tmp | \
sed 's/NosAzo0708_//g' | \
awk -F "\t" '{split($1,A,","); for (i=2;i<=length(A);++i) { if (A[i]<A[i-1]) {print $0} } }'
## Now sort
sort -k1,1 -k3,3 hogs_NosAzo0708.tmp | \
## And keep those hog pairs that match perfectly for NosAzo0708 loci
awk -F "\t" '
BEGIN {prev=""};
$1==prev {print $2"\t"prev_hog"\t"$3};
{prev=$1; prev_hog=$2}
' \
> hog_pairs.tsv
## wc -l hog_pairs.tsv = 3320 3309
## wc -l hogs_NosAzo0708.tmp = 6734 6745
## 6734/2-3320=47 -> 47/(6734/2)=1.4% lost due to loci only appearing in one hog set or to inconsistent multi-locus hogs
## Double check always in the same order
cat hog_pairs.tsv | \
awk '{print $3}' | \
sort -u
## Only ~s. Perfect. ~s correspond to this pipeline which is now in the first column.
sed -i 's/\t~$//' hog_pairs.tsv


##Now count how many loci from each genome are present in each N9.HOG
infile=N9_noOutgroups.tsv
head -n 1 $infile > header.tmp
tail -n +2 $infile | \
awk -F "\t" '{printf $1"\t"$2"\t"$3; for (i=4; i<=NF; ++i) {split($i,A,", "); printf "\t"length(A)}; print ""}'  \
> body.tmp
cat header.tmp body.tmp > n_intact.tsv
rm *.tmp

## Split n_intact.tsv into free living and T. az.
infile=n_intact.tsv
echo -e "\
Taz,~/Taz.shortened.tsv
free_living,~/subtree_freeLiving.shortened.tsv\
" | \
sed 's/,/\t/g' | \
while read outfile_prefix prefix_list; do

outfile=${outfile_prefix}.n_intact.tsv
cat $infile | awk -F "\t" '{print $1}' > $outfile

cat $prefix_list | awk '{print $2}' | \
while read prefix; do

cat $infile | \
awk -v prefix=$prefix -F "\t" '
NR==1 { for (i=4;i<=NF;++i) { if ($i==prefix) {c=i} } };
{ print $c }
' \
> this_col.tmp
paste $outfile this_col.tmp > tmp
mv tmp $outfile
rm this_col.tmp

## Finish loop through prefixes
done
## Finish loop through Taz / free-living
done

##########################
### Assign pseudogenes ###
##########################

## set variables
outdir=~/data/chloroplast_genes/assign_pseudogenes
N=N9
prefix_tsv=~/subtree.shortened.tsv
orthofinder_downstream_dir=~/data/chloroplast_genes/orthofinder_out/for_downstream
pseudofinder_outdir=~/data/chloroplast_genes/pseudofinder_out

## Copy the main_pipeline pseudofinder results for the free-living + NosAzo0708
main_pseudofinder_outdir=~/data/main_pipeline/pseudofinder_out
cat ~/subtree_freeLiving.shortened.tsv | \
awk '{print $1}; END {print "NostocAzollae0708"}' | \
while read prefix; do
cp -r ${main_pseudofinder_outdir}/${prefix} ${pseudofinder_outdir}
done

cd ~/sandbox
## Call assign_pseudogenes_to_hogs.sh !
bash -i ~/scripts/helper/assign_pseudogenes_to_hogs.sh \
-p $pseudofinder_outdir \
-f $orthofinder_downstream_dir \
-n N9 \
-l $prefix_tsv \
-o $outdir


##################################
### Identify predicted-pseudos ###
#### that run off the end of #####
########### the contig ###########
##################################

## Main pipeline
pseudofinder_dir=~/data/main_pipeline/pseudofinder_out
assign_pseudos_dir=~/data/main_pipeline/assign_pseudogenes
orthofinder_downstream_dir=~/data/main_pipeline/orthofinder_out/for_downstream
bash -i ~/scripts/helper/identify_end_of_contig_pseudos.sh -p $pseudofinder_dir -a $assign_pseudos_dir -o $orthofinder_downstream_dir

## improve pipeline (this pipeline)
pseudofinder_dir=~/data/chloroplast_genes/pseudofinder_out
orthofinder_downstream_dir=~/data/chloroplast_genes/orthofinder_out/for_downstream
assign_pseudos_dir=~/data/chloroplast_genes/assign_pseudogenes
bash -i ~/scripts/helper/identify_end_of_contig_pseudos.sh -p $pseudofinder_dir -a $assign_pseudos_dir -o $orthofinder_downstream_dir


###############################
##### Combine tables from #####
### this and main pipelines ###
###############################

## subset Taz.intact_pseudo table to only Dijkhuizen (no NosAzo0708)
cd ~/data/chloroplast_genes/orthofinder_out/for_downstream

intable=Taz.n_intact_pseudo_w_end_of_contig.tsv
outtable=Dijk.tsv
cat $intable | awk -F "\t" '{print $1}' > $outtable
## Loop through MAGs, adding a col for each 
cat ~/Dijkhuizen.prefix.list | \
while read prefix; do
cat $intable | \
awk -v prefix=$prefix -F "\t" '
NR==1 { for (i=1;i<=NF;++i) { if ($i==prefix) {c=i; print $i; next} } }
{ print $c };
' \
> this_col.tmp
paste $outtable this_col.tmp > tmp
mv tmp $outtable
done
rm this_col.tmp

## Subset Dijk.tsv table to only the hogs matched above to the main pipeline above
intable=$outtable
outtable=matched_hogs_Dijk.tsv
cp $intable $outtable
head -n 1 $outtable | sed 's/^HOG/HOG\tNosAzo0708/' > header.tmp
cat hog_pairs.tsv >> $outtable
sort -r $outtable | \
## sort should always put the line that bears only the hog names before the line that has the locus counts
awk -F "\t" '
BEGIN {prev=""};
$1==prev {print col2"\t"$0}
{prev=$1; col2=$2};
' \
>> header.tmp
mv header.tmp $outtable
## The second col is a place holder so the cols are in the same order as the old_Taz_table
## It is the N9.HOG from this pipeline, but it will be removed below.

## Check that columns are in the same order between the two tables.
main_pipeline_Taz_table=~/data/main_pipeline/orthofinder_out/for_downstream/Taz.n_intact_pseudo_w_end_of_contig.tsv
improve_pipeline_Taz_table=$outtable
head -n 1 $main_pipeline_Taz_table > tmp
head -n 1 $improve_pipeline_Taz_table > tmp.tmp
cmp tmp tmp.tmp
rm tmp tmp.tmp
## Check successful.
outtable=Taz.n_intact_pseudo_contig_ends_plus_improved.tsv
head -n 1 $main_pipeline_Taz_table > $outtable
tail -n +2 $main_pipeline_Taz_table > tmp
tail -n +2 $improve_pipeline_Taz_table >> tmp 
## Add new /values to lines from the old table
sort -r tmp | \
awk -F "\t" '
$2 ~ "^N9" { split($0,A,"\t"); next };
{
if ($1==A[1]) {
        printf $1"\t"$2"/-/-/-";
        for (i=3;i<=NF;++i) { printf "\t"$i"/"A[i] };
        print "";
        }
else {
        printf $1"\t"$2"/-/-/-";
        for (i=3;i<=NF;++i) { printf "\t"$i"/-/-/-" };
        print "";
}
}
' | \
sort -k1,1 \
>> $outtable
rm tmp

## How many genes now have at least 1 predicted-intact locus that previously had none?
cat $outtable | awk -F "\t" 'BEGIN {c=0}; {for (i=3;i<=NF;++i) {split($i,A,"/"); if (A[1]==0 && A[4]!="-" && A[4]!="0") {c+=1} } }; END {print c}'
## 370 (metaspades, looseFilter) 373 (isospades, loose filter)
## Across how many HOGs?
cat $outtable | awk -F "\t" '{for (i=3;i<=NF;++i) {split($i,A,"/"); if (A[1]==0 && A[4]!="-" && A[4]!="0") {print $0; next} } }' | wc -l
## 163 (metaspades, looseFilter) 178 (isospades, loose filter)

## Look at psbD
#N9.HOG0002331|N9.HOG0002332
#NosAzo0708_01318,NosAzo0708_04619

## Clean up some files
echo -e "improve_pipeline\tmain_pipeline" > tmp
cat hog_pairs.tsv >> tmp
mv tmp hog_pairs.tsv
rm Dijk.tsv matched_hogs_Dijk.tsv
mkdir perfect_hog_matches
mv matched_hogs_this_pipeline.list matched_hogs_main_pipeline.list unmatched_hogs_main_pipeline.list perfect_hog_matches
## Remove the pseudofinder outputs that were copied from the main pipeline
pseudofinder_dir=~/data/chloroplast_genes/pseudofinder_out
cat ~/subtree.shortened.tsv | awk '{print $1}' | grep -vFf ~/Dijkhuizen.prefix.list | \
while read prefix; do rm -r ${pseudofinder_dir}/${prefix}; done
## Move the whole folder
mv ~/data/chloroplast_genes ~/data/main_pipeline/improve_gene_inventory


#################################
### Make T. az. binary tables ###
#################################

## Make a T. az. table with binary values of possibly intact vs. likely not intact
## And another table of present vs. absent
## "gene_loss" dir is a relic of some older code.

main_dir=~/data/main_pipeline/gene_loss
mkdir $main_dir
cd $main_dir

## combine the the 6-slash tables from isospades and metaspades improve_pipelines
## Need to use both 6-slash (instead of 1-6-slash and 1-3-slash, for example) because 6-slash contain only hogs matched to the main_pipeline
intable=Taz.n_intact_pseudo_end_of_contig_ALL_PIPELINES.tsv
metaspades_table=~/data/main_pipeline/improve_gene_inventory/metaspades_looseFilter/orthofinder_out/for_downstream/Taz.n_intact_pseudo_contig_ends_plus_improved.tsv
isospades_table=~/data/main_pipeline/improve_gene_inventory/isospades_looseFilter/orthofinder_out/for_downstream/Taz.n_intact_pseudo_contig_ends_plus_improved.tsv
head -n 1 $metaspades_table > tmp
head -n 1 $isospades_table > tmp.tmp
cmp tmp tmp.tmp
## cmp empty. Good!
awk '{print $1}' $metaspades_table > tmp
awk '{print $1}' $isospades_table > tmp.tmp
cmp tmp tmp.tmp
## cmp empty. Good!
rm tmp tmp.tmp
## Check successful.
## initiate combined table
head -n 1 $metaspades_table > $intable
## cat input tables together
## Add ~ to keep sort consistent
tail -n +2 $metaspades_table > tmp
tail -n +2 $isospades_table | sed 's/\t/~\t/' >> tmp
sort -k1,1 tmp | \
awk -F "\t" '
NR%2==1 { split($0,A,"\t"); next };
{
if ($1==A[1]"~") {
	split($0,B,"\t");
	printf A[1];
        for (i=2;i<=NF;++i) {
		printf "\t"A[i]
		split(B[i],C,"/")
		printf "/"C[4]"/"C[5]"/"C[6]
	};
        print "";
} else {
        print "ERROR: "A[1]" doesnt match "$1
}
}
' \
>> $intable
rm tmp
## Check no errors
grep "ERROR" $intable
## Check that rows unchanged
awk '{print $1}' $metaspades_table > tmp
awk '{print $1}' $intable > tmp.tmp
cmp tmp tmp.tmp
rm tmp tmp.tmp

## Presence table:
head -n 1 $intable > Taz.presence_absence.tsv
cat $intable | \
tail -n +2 | \
awk '
{printf $1};
{
for (i=2;i<=NF;++i) {
	p=0;
	split($i,A,"/");
	for (j=1;j<=length(A);++j) {
		if (A[j]!=0 && A[j]!="-") {p=1};
	}
	printf "\t"p;
}
print ""
}
' \
>> Taz.presence_absence.tsv

## Make a table for intact_vs_not
head -n 1 $intable > Taz.intact_vs_not.tsv
cat $intable | \
tail -n +2 | \
awk '
{printf $1};
{
## First, check that at there is at least one predicted-intact locus for this hog in the T. az.
at_least_one_intact=0
for (i=2;i<=NF;++i) {
        split($i,A,"/");
        ## replace - with 0 so code below is more readable.
        for (j=1;j<=length(A);++j) { if (A[j]=="-") {A[j]=0} };
        ## if intact locus in any of the three pipelines, intact
        if (A[1]!=0 || A[4]!=0 || A[7]!=0) {at_least_one_intact=1}
}
## if none predicted intact, then all are 0s.
if (at_least_one_intact==0) {
	for (i=2;i<=NF;++i) {
		printf "\t0"
	}
print ""
} else {
## Otherwise, do the more involved check
	for (i=2;i<=NF;++i) {
		p=0;
		split($i,A,"/");
		## replace - with 0 so code below is more readable.
		for (j=1;j<=length(A);++j) { if (A[j]=="-") {A[j]=0} };
		## if intact locus in any of the three pipelines, intact
		if (A[1]!=0 || A[4]!=0 || A[7]!=0) {p=1}
		## if pipeline with the most pseudos has end-of-contig pseudos
		else if ( A[2]>A[5] && A[2]>A[8] ) { if ( A[3]>0 ) {p=1} }
		else if ( A[5]>A[2] && A[5]>A[8] ) { if ( A[6]>0 ) {p=1} }
		else if ( A[8]>A[5] && A[8]>A[2] ) { if ( A[9]>0 ) {p=1} }
		else if ( A[2]==A[5] && A[5]==A[8] ) { if ( A[3]>0 && A[6]>0 && A[9]>0 ) {p=1} }
		else if ( A[2]==A[5] && A[5]>A[8] ) { if ( A[3]>0 && A[6]>0 ) {p=1} }
		else if ( A[2]==A[8] && A[2]>A[5] ) { if ( A[3]>0 && A[9]>0 ) {p=1}  }
		else if ( A[5]==A[8] && A[8]>A[2] ) { if ( A[6]>0 && A[9]>0 ) {p=1}  }
		printf "\t"p;
	}
	print ""
}
}
' \
>> Taz.intact_vs_not.tsv

## Check first cols look good.
awk '{print $1}' Taz.presence_absence.tsv > tmp
awk '{print $1}' Taz.intact_vs_not.tsv > tmp.tmp
cmp tmp tmp.tmp
awk '{print $1}' $intable > tmp.tmp
cmp tmp tmp.tmp
rm tmp tmp.tmp

#################################
### Make Nostoc binary tables ###
#################################

## 55 hogs have a total of 60 instances of end-of-contig-pseudogenes in the free-living Nostocales
## But because these are complete genomes, we are still confident just treating these as pseudogenes.
intable=~/data/main_pipeline/orthofinder_out/for_downstream/free_living.n_intact_and_pseudo.tsv
## make free-living presence_absence table
head -n 1 $intable > free_living.presence_absence.tsv
cat $intable | \
tail -n +2 | \
awk '
{printf $1};
{
for (i=2;i<=NF;++i) {
        p=0;
        split($i,A,"/");
        for (j=1;j<=length(A);++j) {
                if (A[j]!=0 && A[j]!="-") {p=1};
        }
        printf "\t"p;
}
print ""
}
' \
>> free_living.presence_absence.tsv
## Add it to the Taz.presence_absence table
## Check That the rows are all the same
cat free_living.presence_absence.tsv | awk '{print $1}' > tmp
cat Taz.presence_absence.tsv | awk '{print $1}' > tmp.tmp
cmp tmp tmp.tmp
## Success.
rm tmp tmp.tmp
cat free_living.presence_absence.tsv | awk ' { for (i=2;i<NF;++i) {printf $i"\t"}; print $NF} ' > tmp
paste Taz.presence_absence.tsv tmp > all.presence_absence.tsv
rm tmp

## make free-living intact vs  not
head -n 1 $intable > free_living.intact_vs_not.tsv
cat $intable | \
tail -n +2 | \
awk '
{printf $1};
{
for (i=2;i<=NF;++i) {
        split($i,A,"/");
	if (A[1]==0) {printf "\t"0}
	else {printf "\t"1}
}
print ""
}
' \
>> free_living.intact_vs_not.tsv
## Add it to the Taz.intact_vs_not table
## Check That the rows are all the same
cat free_living.intact_vs_not.tsv | awk '{print $1}' > tmp
cat Taz.intact_vs_not.tsv | awk '{print $1}' > tmp.tmp
cmp tmp tmp.tmp
## Success.
rm tmp tmp.tmp
cat free_living.intact_vs_not.tsv | awk ' { for (i=2;i<NF;++i) {printf $i"\t"}; print $NF} ' > tmp
paste Taz.intact_vs_not.tsv tmp > all.intact_vs_not.tsv
rm tmp


#######################
### Check effect of ###
#### this pipeline ####
#######################

## How many HOGs/loci were affected by this whole script???
## Check how much of a difference there is between
## The intact binary table after orthofinder
## And the one that used the multiple assembly approaches.

cd ~/data/main_pipeline/gene_loss

## $orthofinder_intact is not binary, but anything other than a 0 can be considered a 1. It was generated using only the main assemblies.
## Taz.intact_vs_not.tsv is binary after the pipeline to look for pseudogenes that might be intact but broken by the partial assembly.
orthofinder_intact=~/data/main_pipeline/orthofinder_out/for_downstream/Taz.n_intact.tsv

echo -e "HOG\tGenomes" > intactness_unknown.tsv
paste Taz.intact_vs_not.tsv $orthofinder_intact | \
awk -F "\t" '
NR==1 { for (i=1;i<=NF/2;++i) {A[i]=$i}; next };
{
p=0;
for (i=2;i<=NF/2;++i) {
	if ( $i==1 && $(i+NF/2)==0 ) {
		if (p==0) {printf $1; p=1}
		printf ","A[i]
	}
}
if (p==1) {print ""; p=0}
}
' | \
sed 's/,/\t/' \
>> intactness_unknown.tsv


## Get the annotations for these
## Use the fact that header is "HOG" for intactness_unknown.tsv and for annotation_file
annotation_file=~/data/main_pipeline/annotations/eggnog-mapper_out/hog_annotations/hog_annotations.tsv
> annotation_col.tmp
cat intactness_unknown.tsv | \
awk '{print $1}' | \
while read hog; do
cat $annotation_file | awk -v hog=$hog -F "\t" '$1==hog' >> annotation_col.tmp
done
cat intactness_unknown.tsv | awk -F "\t" '{print $2}' > genome_col.tmp
paste annotation_col.tmp genome_col.tmp > paste.tmp
mv paste.tmp intactness_unknown.tsv
rm annotation_col.tmp genome_col.tmp


