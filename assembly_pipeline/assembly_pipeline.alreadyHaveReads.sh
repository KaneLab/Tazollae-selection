##Pipeline to assemble T. az. MAGs from Li et al, 2018 short reads
##Subsamples reads so that metaSPAdes can handle without running out of RAM
##Trimmomatic
##filter out reads that align to the reference Azolla host genome (bwa mem) 
##metaSPAdes
##metabat2
##Some custom bash to look at the bins
##T. az. bins are then analyzed for completeness and contamination using checkM on kbase
##mmseqs2 will be used to taxonomically ID the other bins
##A second round of read filtering may be added after ID'ing bins if any belong to streptophyta.

print_usage() {
  printf "Usage: requires cluster_list -l, indir -i, outdir -o"
}

while getopts 's:p:a:r:' flag; do
  case "${flag}" in
    s) sra_num="${OPTARG}" ;;
    p) prefix="${OPTARG}" ;;
    a) assembler="${OPTARG}" ;;
    r) sampleRatio="${OPTARG}" ;;
    *) print_usage
       exit 1 ;;
  esac
done

mkdir ~/data/assembly_workflow
mkdir ~/data/assembly_workflow/metaSPAdes
mkdir ~/data/assembly_workflow/SPAdes
mkdir ~/data/assembly_workflow/megahit

#Make the prefix an input.
#sra_num=SRR6480231
#prefix=caroliniana1
#assembler=metaSPAdes
#sampleRatio=3

##Download the reads!
echo -e "\
######################################\n\
### Downloading SRA! #################\n\
######################################\
"
#conda activate sra-tools
#mkdir ~/data/SRA_files/${prefix}
#cd ~/data/SRA_files/${prefix}
#prefetch $sra_num && \
#vdb-validate $sra_num > vdb-validate.out.txt && \
#fasterq-dump $sra_num #Don't need --split-files
##fasterq-dump can also be run without prefetch, but this is slower and maybe more prone to crash
#rm -r $sra_num
#gzip *.fastq

##Subsample reads so assembly runs!
echo -e "\
######################################\n\
### subsample reads! #################\n\
######################################\
"
mod_denominator=$((4*$sampleRatio))
cd ~/data/SRA_files/$prefix
gunzip *.gz
ls $SRA_*.fastq | \
while read file; do
echo $file
outfile=${file%_*}_${sampleRatio}xSubSampled_${file#*_}
if (( $sampleRatio==1 )); then
#using the whole sequence set! No filtering necessary!
cp $file $outfile
else
#filter!
cat $file | awk -v md=$mod_denominator -F "\n" '(NR - 1) % md < 4 {print $0}' > $outfile
fi
gzip $file
done

##Trim reads
##Run trimmomatic:
echo -e "\
######################################\n\
### trimmomatic! #####################\n\
######################################\
"
conda deactivate
conda activate trimmomatic
#mkdir ~/data/assembly_workflow/trimmomatic_out
mkdir ~/data/assembly_workflow/trimmomatic_out/${prefix}
cd ~/data/assembly_workflow/trimmomatic_out/${prefix}
adapterFile=~/miniconda3/envs/trimmomatic/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa
indir=~/data/SRA_files/${prefix}
in_1=${indir}/*_1.fastq
in_2=${indir}/*_2.fastq
outputBase=${prefix}_trimmed.fq.gz
trimmomatic PE -threads 20 $in_1 $in_2 -baseout $outputBase ILLUMINACLIP:${adapterFile}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 && \
gunzip *.gz
##Filter out reads that align to the reference plastid and plant nuclear genomes.
#Make a reference to align against
#mkdir ~/data/reference
#cd ~/data/reference
#cat Azolla_filiculoides.cp_genome_v1_4.fasta Azolla_filiculoides.genome_v1.2.fasta > plantRefGenome_inclChlp.fasta
#gzip Azolla_filiculoides.cp_genome_v1_4.fasta
#gzip Azolla_filiculoides.genome_v1.2.fasta

##Align reads and output SAM file.
#mkdir ~/data/assembly_workflow/filter_PlantChlpReads
echo -e "\
######################################\n\
### Align and filter reads! ##########\n\
######################################\
"
mkdir ~/data/assembly_workflow/filter_PlantChlpReads/${prefix}
cd ~/data/assembly_workflow/filter_PlantChlpReads
####SCREEN!###
conda deactivate
conda activate samtools
reference=~/data/reference/plantRefGenome_inclChlp.fasta
outdir=~/data/assembly_workflow/filter_PlantChlpReads/${prefix}
indir=~/data/assembly_workflow/trimmomatic_out/${prefix}
#bwa index $reference
##I should re-write this so that the paired reads are done together.
for infile in ${indir}/*.fq; do
suffix=${infile##*_}
suffix=${suffix%.fq}
bwa mem $reference $infile > $outdir/trimmedReads2Ref_${suffix}.sam
done
#Filter out paired reads that are unmapped and partner is unmapped
#and unpaired reads that are unmapped.
## samtools view -f 4 is equivalent to awk ' ($2 % 8) >= 4 {print $0}'
## samtools view -f 4,8 is equivalent to awk ' ($2 % 16) >= 12 {print $0}'
## The first prints all unmapped reads.
## The second prints all reads where both read and partner are unmapped.
## If doing this with paired reads, samtools view -f 4,8 could be used instead of the awk call.
## See code notes from Jan 26.
##The @ symbol before and space after $1 in awk is to make grep -Ff work better.
cd $outdir
for infile in trimmedReads2Ref_*.sam; do
echo $infile
suffix=${infile#trimmedReads2Ref_}
suffix=${suffix%.sam}
echo $suffix
cat $infile | awk ' ($2 % 8) >= 4 {print "@"$1" "}' > Unmapped_${suffix}.list
done && \
##gzip the original sam files now that no longer needed.
gzip trimmedReads2Ref_*.sam

##If ran 1P and 2P separate, remove reads that did not hit both:
#Something in this pipe erases the space at the end of the line
#So I add it back in the awk call
cat Unmapped_1P.list Unmapped_2P.list | \
sort | uniq -c | awk '$1 == 2 {print $2" "}' \
> both_unmapped.list
##make fastqs of only the filtered hits.
#First do some file renaming of trimmomatic outputs
cd ~/data/assembly_workflow/filter_PlantChlpReads/${prefix}
#indir=~/data/assembly_workflow/trimmomatic_out
#mkdir ${indir}/${prefix}
#mv ${indir}/${prefix}*.fq ${indir}/${prefix}
indir=~/data/assembly_workflow/trimmomatic_out/${prefix}
for file in $indir/*.fq; do mv $file ${indir}/${file#*${prefix}_}; done
#-m 1 saves time in grep, although it could also cover up bugs if multiple matches would have been found
grep --no-group-separator -A3 -Ff both_unmapped.list ${indir}/trimmed_1P.fq > ${indir}/filtered_1P.fq && \
grep --no-group-separator -A3 -Ff both_unmapped.list ${indir}/trimmed_2P.fq > ${indir}/filtered_2P.fq && \
grep --no-group-separator -A3 -Ff Unmapped_1U.list ${indir}/trimmed_1U.fq > ${indir}/filtered_1U.fq && \
grep --no-group-separator -A3 -Ff Unmapped_2U.list ${indir}/trimmed_2U.fq > ${indir}/filtered_2U.fq && \
gzip ${indir}/trimmed*.fq


####Assemble!
echo -e "\
######################################\n\
### Assemble! ########################\n\
######################################\
"
indir=~/data/assembly_workflow/trimmomatic_out/${prefix}
mkdir ~/data/assembly_workflow/$assembler/assemblies
cd ~/data/assembly_workflow/$assembler/assemblies
conda deactivate
if [[ $assembler == SPAdes ]]; then
conda activate spades
~/miniconda3/envs/spades/bin/spades.py -o $prefix -1 ${indir}/filtered_1P.fq -2 ${indir}/filtered_2P.fq -s ${indir}/filtered_1U.fq -s ${indir}/filtered_2U.fq
elif [[ $assembler == metaSPAdes ]]; then
conda activate spades
~/miniconda3/envs/spades/bin/spades.py --meta -o $prefix -1 ${indir}/filtered_1P.fq -2 ${indir}/filtered_2P.fq -s ${indir}/filtered_1U.fq -s ${indir}/filtered_2U.fq
elif [[ $assembler == megahit ]]; then
conda activate megahit
megahit -1 ${indir}/filtered_1P.fq -2 ${indir}/filtered_2P.fq -r ${indir}/filtered_1U.fq,${indir}/filtered_2U.fq -o $prefix
#rename contigs file to be the same as the spades outputs to make downstream easier
mv $prefix/final.contigs.fa $prefix/scaffolds.fasta
else
echo "Check assembler variable"
fi

##Align reads to assembly before binning!
#Set variables
echo -e "\
######################################\n\
### Align before binning! ############\n\
######################################\
"
outdir=~/data/assembly_workflow/$assembler/metabat2_out
mkdir $outdir
outfile=${prefix}.pairedReads2assembly.sam
reads_dir=~/data/assembly_workflow/trimmomatic_out/${prefix}
reads_1=filtered_1P.fq
reads_2=filtered_2P.fq
subject=~/data/assembly_workflow/$assembler/assemblies/${prefix}/scaffolds.fasta
# run bwa mem
conda deactivate
conda activate samtools
bwa index $subject
bwa mem -t 20 $subject ${reads_dir}/${reads_1} ${reads_dir}/${reads_2} > ${outdir}/${outfile}
#convert sam to bam file, sort BAM and create an index
cd ~/data/assembly_workflow/$assembler/metabat2_out
samtools view -b -o ${prefix}.pairedReads2assembly.bam -S ${prefix}.pairedReads2assembly.sam && \
samtools sort ${prefix}.pairedReads2assembly.bam -o ${prefix}.pairedReads2assembly.sorted.bam && \
samtools index ${prefix}.pairedReads2assembly.sorted.bam

##Bin!!!
##Get depth file with metabat2
echo -e "\
######################################\n\
############### Bin! #################\n\
######################################\
"
conda deactivate
conda activate metabat2
cd ~/data/assembly_workflow/$assembler/metabat2_out
infile=${prefix}.pairedReads2assembly.sorted.bam
jgi_summarize_bam_contig_depths --outputDepth ${prefix}.metabat2.depth.txt $infile
mkdir ${prefix}
mv *.* ${prefix}
#run metabat2
assembly=~/data/assembly_workflow/$assembler/assemblies/${prefix}/scaffolds.fasta
cd $prefix
depth_file=${prefix}.metabat2.depth.txt
metabat2 -i $assembly -a ${prefix}.metabat2.depth.txt -o ${prefix}
#calculate GC content for metabat2 bins
echo -e "bin\ngenomeSize\nGC%" > GC_content.tsv
for file in *.fa; do
echo $file;
cat $file | grep -v "^>" > tmp.tmp
n_char=$((60*$(cat tmp.tmp | wc -l)));
echo -e "${file%.fa}\n$n_char" >> GC_content.tsv
n_GC=$(cat tmp.tmp | tr -cd 'GC' | wc -c)
awk -v n_GC=$n_GC -v n_char=$n_char 'BEGIN {print n_GC/n_char}' >> GC_content.tsv
rm tmp.tmp
done
cat GC_content.tsv | \
awk -F "\n" 'NR % 3 == 0 {print $0; next}; {printf $0"\t"}' | \
sort -t $'\t' -k3,3 -n > tmp.tmp
mv tmp.tmp GC_content.tsv
#Blast bins against reference genome to find plasmids and chromosomes.
conda deactivate
db_dir=~/data/reference
echo -e "Taz_chrom.ffn\nTaz_plas1.ffn\nTaz_plas2.ffn" > db.list.tmp
cat db.list.tmp | \
while read db; do
#makeblastdb -in $db_dir/$db -dbtype nucl 
for file in *.fa; do 
bin=${file##*/}
echo $bin;
outfile=${bin%.fa}.${db%.ffn}.blastout.tsv;
query=$file; blastn -query $query -db ${db_dir}/${db} -outfmt 6 -out $outfile;
done;
done
rm *.tmp
echo "n_contigs" > n_contigs.col.tmp
echo "n_contigs_chrom" > chrom.col.tmp
echo "n_contigs_plas1" > plas1.col.tmp
echo "n_contigs_plas2" > plas2.col.tmp

tail -n +2 GC_content.tsv | awk '{print $1}' | \
while read bin; do
echo $bin;
file=${bin}.fa
grep -c "^>" $file >> n_contigs.col.tmp
cat ${bin}.Taz_chrom.blastout.tsv | awk '{print $1}' | sort -u | wc -l >> chrom.col.tmp
cat ${bin}.Taz_plas1.blastout.tsv | awk '{print $1}' | sort -u | wc -l >> plas1.col.tmp
cat ${bin}.Taz_plas2.blastout.tsv | awk '{print $1}' | sort -u | wc -l >> plas2.col.tmp
done
paste GC_content.tsv n_contigs.col.tmp chrom.col.tmp plas1.col.tmp plas2.col.tmp > bin_analysis.tsv
rm GC_content.tsv n_contigs.col.tmp chrom.col.tmp plas1.col.tmp plas2.col.tmp
##Get rid of all the empty blast files
for file in *blastout.tsv; do
if (( $(wc -l < $file) == 0 )); then
rm $file
fi
done

#blast all of the unbinned contigs with len >=2000
#remove contigs with len < 2500
#This is the default for metabat2
cd ~/data/assembly_workflow/${assembler}/metabat2_out/${prefix} 
cp ~/data/assembly_workflow/$assembler/assemblies/${prefix}/scaffolds.fasta ./
if [[ $assembler == SPAdes || $assembler == metaSPAdes ]]; then
#This works because SPAdes organizes sthe contigs from longest to shortest
awk -F "_" '$0 ~ "^>" && $4 < 2500 {exit}; {print $0}' scaffolds.fasta > ${prefix}.unbinned.tmp
elif [[ $assembler = megahit ]]; then
cat scaffolds.fasta | awk -F "len=" '
BEGIN {p=0};
{ if ($0 ~ "^>") { if ($2 < 2500) {p=0} else if ($2 >= 2500) {p=1} } };  
p==1 {print $0}' \
> ${prefix}.unbinned.tmp
fi
rm scaffolds.fasta
#unwrap to just one line per read (so not the normal 2-line unwrap)
cat ${prefix}.unbinned.tmp | \
awk ' {if ($0 ~ "^>") {printf "\n"$0"\t"} else { printf $0} }; END {printf "\n"} ' | \
tail -n +2 \
> ${prefix}.unbinned.unwrapped.tmp
#list binned contigs
> binned_contigs.tmp
for file in ${prefix}.*.fa; do
echo $file;
grep "^>" $file | sed 's/$/ /' >> binned_contigs.tmp
done
##Remove binned contigs
grep -v -Ff binned_contigs.tmp ${prefix}.unbinned.unwrapped.tmp > tmp.tmp
##Put into normal 2-line unwrapped
sed 's/\t/\n/g' tmp.tmp > ${prefix}.unbinned.fa
outfile=unbinned_contigs.Taz_complete.blastout.tsv;
db=~/data/reference/Taz_complete.ffn
query=${prefix}.unbinned.fa
blastn -query $query -db ${db} -outfmt 6 -out $outfile;
#look at hits with len >= 1000
awk '$4 > 1000 {print $0}' $outfile > unbinned_contigs.Taz_complete.minlen1K.blastout.tsv
rm *.tmp



###############################################
## Filter out reads before a second assembly ##
###############################################

##Filter out 
