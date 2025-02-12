##Run mmseqs taxonomy on all of the bins from all of the first assemblies for the re-assembled Li MAGs
##Within each bin, all contigs are
mkdir ~/data/mmseqs2_out/all_assemblies
cd ~/data/mmseqs2_out/all_assemblies

##first, have to create database for mmseqs2 to use to search nr which was already downloaded.
##Set up taxonomyDB for MMseqs
##Install mmseqs2
conda create -n mmseqs2
conda activate mmseqs2
conda install -c bioconda mmseqs2
conda activate mmseqs2
cd ~/data/nr
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
mkdir taxonomy && tar -xxvf taxdump.tar.gz -C taxonomy
rm taxdump.tar.gz
blastdbcmd -db nr -entry all > nr.fna
blastdbcmd -db nr -entry all -outfmt "%a %T" > nr.fna.taxidmapping
mmseqs createdb nr.fna nr.fnaDB && \
mmseqs createtaxdb nr.fnaDB tmp --ncbi-tax-dump taxonomy/ --tax-mapping-file nr.fna.taxidmapping


##Combine all the bins from all the meta-assemblies into one file for mmseqs taxonomy. For each bin, combine all contigs into a single sequence, with different contigs separated by runs of 300 N's.
##Note that for some of the meta-assemblies, there are bins that look no bueno in terms of being way too large.
##Note that this sed call leaves the line before the new contig as potentially shorter than full line length. I don't think this should matter.

outfile=allAssemblies.allBins.fa
> $outfile

cat ~/Li.prefix.list | \
while read prefix; do

ls ~/data/assembly_workflow/metaSPAdes/metabat2_out/${prefix}/${prefix}.[0-9]*.fa | \
while read file; do
bin=${file%.fa}
bin=${bin##*/}
echo ">${bin}" >> $outfile
tail -n +2 $file | \
sed "s/>.*$/$(printf '%.0sN' {1..60})\n$(printf '%.0sN' {1..60})\n$(printf '%.0sN' {1..60})\n$(printf '%.0sN' {1..60})\n$(printf '%.0sN' {1..60})\n$(printf '%.0sN' {1..60})/" \
>> $outfile
done

done
##Done creating combined fasta input for mmseqs taxonomy module.

##Run mmseqs taxonomy
cd ~/data/mmseqs2_out/all_assemblies
query=allAssemblies.allBins.fa
queryDB=${query}.DB
seqTaxDB=~/data/mmseqs2_out/nr_for_mmseqs2/nr_for_mmseqs2
mmseqs createdb $query $queryDB
mmseqs taxonomy $queryDB $seqTaxDB taxonomyResult tmp
mmseqs createtsv $queryDB taxonomyResult taxonomyResult.tsv
sort -t $'\t' -k1,1 taxonomyResult.tsv > taxonomyResult.sorted.tsv

