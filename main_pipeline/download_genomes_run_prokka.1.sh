## Initial steps of the pipeline
## Download genomes:
## 7 T. az. MAGs from Dijkhuizen et al, 2021;
## RefSeq: 1 T. az. (Nostoc azollae 0708), 47 free-living Nostocales, 2 free-living Gloeobacter outgroup 
## Run prokka on all genomes to get Prodigal ORF predictions in a pseudofinder-friendly format

##Starting from assembled genomes or MAGs with all .fna files in one folder:
##All scripts and all code below should be written with full pathnames.
##This allows me to call the scripts from anywhere, which I do from an empty "sandbox" folder.
##However, pathnames may certainly change as data gets moved around.
##A list of names of MAGs and genomes is contained in prefix.list which is used
##to loop function calls.

#####################################
## Download free-living Nostocales ##
#### and Gloeobacter (outgroup) #####
###### genomes from RefSeq ##########
#####################################

##Download the 48 Nostocales genomes from ncbi that are:
##refseq, complete, annotated, 2010-present
## --assembly-source refseq (other option genbank or no flag) restricts to only refseq assemblies.
##Some notes in project notes about this tool vs. online browser.
##Note that this is not done with the most up-to-date version which seems to have slightly different commands.

mkdir /home/liam/data/ncbi_genomes
cd /home/liam/data/ncbi_genomes
#Note there is an update version of this package and the commands might be slightly different
#https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/datasets_download_genome_taxon/
conda activate ncbi_datasets
datasets download genome taxon nostocales --assembly-level complete_genome --released-since 01/01/2010 --assembly-source refseq --annotated --include-gbff
unzip ncbi_dataset.zip
rm ncbi_dataset.zip 
##Get the names for each accession number and save in 2-col tsv accessions_names.tsv
##col 1 is accession, col 2 is name.
cd ncbi_dataset/data
ls | grep ^GCF > accessions.tmp
rm name.tmp
while read acc; do
head -n 1 ${acc}/${acc}*_genomic.fna | \
sed 's/\,.*$//' | \
awk -F " " '{for (i=2; i<=NF; ++i) {printf $i" "}; print "" }' | \
sed 's/chromosome.*$//' >> name.tmp;
done < accessions.tmp
paste -d "\t" accessions.tmp name.tmp > accessions_names.tsv
rm *.tmp
#-----------------------------
#Add a third column to accessions_names.tsv by hand that is camelCase,
#with no punctuation.
#-----------------------------
#Rearrange the folders
cd /home/liam/data/ncbi_genomes
mv ncbi_dataset/* ./
rmdir ncbi_dataset
##Rename the nostocales genomes with the camelCase name from above.
cd /home/liam/data/ncbi_genomes/data
for folder in $(ls -d GCF*); do 
new_name=$(awk -v acc=$folder -F "\t" '$1==acc {print $3}' accessions_names.tsv);
mv ${folder}/${folder}*.fna ${folder}/${new_name}.fna
mv $folder $new_name;
done

##Download the outgroup genomes from refseq
##Chose gloeobacter because it is the sister clade to the rest of the cyanos, I believe
cd /home/liam/data
mkdir outgroup_refseq
conda activate ncbi_datasets
datasets download genome taxon "Gloeobacter" --assembly-level complete --released-after 01/01/2010 --assembly-source refseq --annotated --include genome,protein,cds,gff3,gbff
unzip ncbi_dataset.zip
rm ncbi_dataset.zip
mv ncbi_dataset/* outgroup_refseq
rmdir ncbi_dataset
cd outgroup_refseq
mv data/* ./
rmdir data
cd /home/liam/data
rm README.md
#Rename files and folders with genome names determined from the .fna headers
#make a 2-col file with first col the GCF and second col the new name/prefix
cat name_conversion.tsv | while read a b; do mv ${a}/${a}*.fna ${a}/${b}.fna; mv $a $b; done
#Make it a prefix list!
awk -F "\t" '{print $2}' ${outgroup_dir}/name_conversion.tsv > ~/outgroup.prefix.list





################
### Prokka #####
################

conda deactivate
conda activate prokka
mkdir /home/liam/data/main_pipeline
cd /home/liam/sandbox
outdir=/home/liam/data/main_pipeline/prokka_out
mkdir $outdir

## Before running prokka, some of the names that I am using for RefSeq genomes (as in /home/liam/Nostocales_incl0708.prefix.list)
## are too long. Prokka needs them to be max 23 characters.
## I therefore made a 2-col list for /home/liam/Nostocales_incl0708.prefix.list and /home/liam/outgroup.prefix.list
## with the full camelCase name that I am using and then a shortened version of that name.
## in /home/liam/Nostocales_incl0708.shortened.tsv /home/liam/outgroup.shortened.tsv
## and confirm with sort | uniq -c that all the shortened names are unique.

## The RefSeq Nostocales (including NoAz0708), MAGs, and RefSeq outgroup genomes are in different directories
## So loop through the three indirs and lists of genomes
echo -e "\
/home/liam/data/refseq_Nostocales_genomes/data\t\
/home/liam/Nostocales_incl0708.prefix.list
/home/liam/data/fastas/Dijkhuizen_genomes\t\
/home/liam/Dijkhuizen.prefix.list
/home/liam/data/outgroup_refseq\t\
/home/liam/outgroup.prefix.list\
" | \
while read indir genome_list; do

## Loop through list of genomes
cat $genome_list | \
while read prefix; do
echo $prefix;

## Get the inputfile
if [[ -f $indir/${prefix}.fna ]]; then
file=$indir/${prefix}.fna
else
file=$indir/${prefix}/${prefix}.fna
fi

## Get the shortened prefix (if there is one)
if [[ -f ${genome_list%.prefix.list}.shortened.tsv  ]]; then
short_prefix=$( grep $prefix ${genome_list%.prefix.list}.shortened.tsv | awk -F "\t" '{print $2}' )
else
short_prefix=$prefix
fi

## Modify the headers because some are too long for prokka
cat $file | sed -E "s/^>(\S+).*$/>${prefix}_\1/" > infile.tmp

prokka --cpus 20 --compliant --rfam --kingdom Bacteria --outdir ${outdir}/$prefix --locustag $short_prefix --prefix $short_prefix infile.tmp;
rm infile.tmp
done

## End loop through indirs and genome_lists
done
