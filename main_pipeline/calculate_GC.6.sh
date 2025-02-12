##Get GC content of all 48 subtree genomes.
## As well as # ORFs, # intact
cd /home/liam/data
outfile=GC_content.tsv
echo -e "genome\tGC%\tN%\tGC_n\tAT_n\tN_n\ttotal_noNs\tgenome_size_(bp)\tintact_genes\ttotal_length_intact_genes_(bp)\tpseudogenes" > $outfile
echo "-----" >> $outfile

cat /home/liam/subtree.shortened.tsv | \
while read long_prefix short_prefix; do
file=/home/liam/data/main_pipeline/prokka_out/${long_prefix}/${short_prefix}.fna
echo $short_prefix >> $outfile

cat $file | \
grep --no-group-separator -v "^>" | tr -d "\n" > tmp
characters=$( wc -c < tmp )
total=$(tr -cd "GCgcATat" < tmp | wc -c)
n=$(tr -cd "nN" < tmp | wc -c)
gc=$(tr -cd "GCgc" < tmp | wc -c)
at=$(tr -cd "ATat" < tmp | wc -c)

if (( $gc + $at == $total && $total + $n == $characters )); then
#Use awk for division instead of native bash because of decimal places.
awk -v gc=$gc -v total=$total -v n=$n 'BEGIN {print 100*gc/total "\t" 100*n/total}' >> $outfile
else
echo "error in GC calculation" >> $outfile
fi
echo -e "$gc\t$at\t$n\t$total\t$characters" >> $outfile

## Get the number of predicted intact genes
cat /home/liam/data/main_pipeline/orthofinder_out/input_fastas/${short_prefix}.fasta | grep -c "^>" >> $outfile
## Get the total length of predicted intact genes
cat /home/liam/data/main_pipeline/pseudofinder_out/${long_prefix}/${long_prefix}_intact.ffn | grep -v "^>" | tr -d "\n" | wc -c >> $outfile
## Get the number of predicted pseudogenes
cat /home/liam/data/main_pipeline/pseudofinder_out/${long_prefix}/${long_prefix}_pseudos.gff | grep -v "^#" | wc -l >> $outfile
##Add a demarcation for where the analysis for this genome ends
echo "-----" >> $outfile
done

cat $outfile | tr "\n" "\t" | sed 's/\t-----\t/\n/g' > tmp
mv tmp $outfile
sort -nk2,2 $outfile > tmp
mv tmp $outfile

## Split the table in 2, one with nucleotide counts, the other with GC%, genome size, ORFs, intact
cat $outfile | awk -F "\t" '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > genome_nucleotide_counts.tsv
cat $outfile | \
awk -F "\t" '
{
	printf $1"\t"$2"\t"$8"\t"$9"\t"$11"\t"$10"\t";
	if (NR==1) { printf "%_nucleotides_in_intact_genes\t"} else { printf 100*$10/$8"\t" }
	if (NR==1) { print "mean_length_intact_genes_(bp)" } else { print $10/$9 }
}
' > tmp
mv tmp $outfile

## Get a 2-col tsv of genome and genome # intact ORFs (subtree-only)
cat GC_content.tsv | \
grep -Ff ~/subtree.prefix.list | \
awk '{print $1"\t"$4}' \
> n_intact_ORFs.tmp

pruned_tree=/home/liam/data/main_pipeline/hyphy_out/concatenated_msas/single_copy_HOGs.pruned_tree.txt
R_outfile=n_intact_ORFs_tree.txt
table_path=n_intact_ORFs.tmp
Rscript --vanilla /home/liam/scripts/helper/ape_ace_continuous.R table_path=$table_path tree_path=$pruned_tree outfile=$R_outfile
rm n_intact_ORFs.tmp

## Get itol metadata formatted
cat ~/iTol_metadata_template.txt $R_outfile > genome_size_intact_ORFs_iTOL_metadata.txt

## Get the minimum, median, and maximum for each column
## For free-living and for T. az.
## Loop free-living, T. az.
echo "/home/liam/subtree_freeLiving.short_prefix.list,free_living
/home/liam/Taz.short_prefix.list,Taz" | \
sed 's/,/\t/g' | \
while read prefix_list group; do
echo "quantity,min,med,max" | \
sed 's/,/\n/g' \
> ${group}_genome_stats.tsv
## Loop through columns
n_col=$( cat GC_content.tsv | awk '{print NF}' | sort -u )
col=2
while (( $col <= $n_col )); do
n_genomes=$( wc -l < $prefix_list )
cat GC_content.tsv | head -n 1 | awk -v col=$col '{print $col}' > this_col_out.tmp
cat GC_content.tsv | \
grep -Ff $prefix_list | \
awk -v col=$col '{print $col}' | \
sort -g | \
awk -v n=$n_genomes -F "\t" '
NR==1 || NR==n {print $0}
NR==n/2 { a=$0 };
NR==n/2+1 { print ($0+a)/2 }
' \
>> this_col_out.tmp
## paste
paste ${group}_genome_stats.tsv this_col_out.tmp > paste.tmp
mv paste.tmp ${group}_genome_stats.tsv
## end col loop
col=$(($col+1))
done
rm this_col_out.tmp
## end prefix_list loop
done

## Paste the two tables together
head -n 1 Taz_genome_stats.tsv | \
sed 's/$/\tgroup/' \
> genome_stats.tsv
tail -n +2 free_living_genome_stats.tsv | \
sed 's/$/\tFree-Living (n=40)/' \
>> genome_stats.tsv
tail -n +2 Taz_genome_stats.tsv | \
sed 's/$/\tT. az. (n=8)/' \
>> genome_stats.tsv
rm Taz_genome_stats.tsv free_living_genome_stats.tsv
## use R (local script genome_stats_table.r) to  add superscripts to values that come from NosAzo0708

## Reformat data in GC_content.tsv for R to create a multi-pane chart. (This is the Supp Mat dot plots)
## col1 = category; col2 = group (clade I, clade II, T. az. MAG, T. az. reference); col3 = value; col4 = rank
outfile=genome_stats_against_rank.tsv
echo -e "category\tgroup\tcomplete\tvalue\trank" > $outfile
## Loop through columns
n_col=$( cat GC_content.tsv | awk '{print NF}' | sort -u )
col=2
while (( $col <= $n_col )); do
cat GC_content.tsv | \
awk -v col=$col -F "\t" '
NR==1 { category=$col; next }
{ print category"\t"$1"\t"$col }
' | \
## sort so always in ascending order
sort -k3,3g | \
## Add a fifth column for the rank of this value in this category
awk '{print $0"\t"NR}' \
>> $outfile
## end col loop
col=$(($col+1))
done

## Convert "_" to " "
sed -i 's/_/ /g' $outfile
## Convert genome prefixes to group names and reference or MAG
## T. az. reference
sed -i 's/NosAzo0708/T. az.\tcomplete/g' $outfile 
## T. az. MAG
cat ~/Taz.short_prefix.list | \
grep -v "NosAzo0708" | \
while read prefix; do
sed -i "s/${prefix}/T. az.\tcontig/g" $outfile
done
## clade I
cat ~/clade_I.short_prefix.list | \
grep -vFf ~/Taz.short_prefix.list | \
while read prefix; do
sed -i "s/${prefix}/clade I (not incl. T. az.)\tcomplete/g" $outfile
done
## clade II
cat ~/clade_II.short_prefix.list | \
while read prefix; do
sed -i "s/${prefix}/clade II\tcomplete/g" $outfile
done


## Reformat for graphs, but this time with genome size as the x axis instead of rank as the x axis.
## col1 = category; col2 = group (clade I, clade II, T. az. MAG, T. az. reference); col3 = value; col4 = genome size
outfile=genome_stats_against_genome_size.tsv
echo -e "category\tgroup\tcomplete\tvalue\tgenome size" > $outfile
## Loop through columns
## 3rd is the genome size which is the x axis, but can just include it in the table for conciseness of code
n_col=$( cat GC_content.tsv | awk '{print NF}' | sort -u )
col=2
while (( $col <= $n_col )); do
cat GC_content.tsv | \
awk -v col=$col -F "\t" '
NR==1 { category=$col; next }
{ print category"\t"$1"\t"$col"\t"$3 }
' | \
## sort so always in ascending order of genome size
sort -k4,4g \
>> $outfile
## end col loop
col=$(($col+1))
done

## Convert "_" to " "
sed -i 's/_/ /g' $outfile
## Convert genome prefixes to group names and reference or MAG
## T. az. reference
sed -i 's/NosAzo0708/T. az.\tcomplete/g' $outfile
## T. az. MAG
cat ~/Taz.short_prefix.list | \
grep -v "NosAzo0708" | \
while read prefix; do
sed -i "s/${prefix}/T. az.\tcontig/g" $outfile
done
## clade I
cat ~/clade_I.short_prefix.list | \
grep -vFf ~/Taz.short_prefix.list | \
while read prefix; do
sed -i "s/${prefix}/clade I (not incl. T. az.)\tcomplete/g" $outfile
done
## clade II
cat ~/clade_II.short_prefix.list | \
while read prefix; do
sed -i "s/${prefix}/clade II\tcomplete/g" $outfile
done

## Want a subset of genome_stats_against_genome_size.tsv and genome_stats_against_rank.tsv
## For manuscript figure.
outfile=genome_stats_for_figure.tsv
head -n 1 genome_stats_against_rank.tsv > $outfile
awk -F "\t" '$1=="GC%"' genome_stats_against_rank.tsv >> $outfile
awk -F "\t" '$1=="genome size (bp)"' genome_stats_against_rank.tsv >> $outfile
awk -F "\t" '$1=="intact genes"' genome_stats_against_genome_size.tsv >> $outfile
awk -F "\t" '$1=="total length intact genes (bp)"' genome_stats_against_genome_size.tsv >> $outfile
awk -F "\t" '$1=="pseudogenes"' genome_stats_against_genome_size.tsv >> $outfile

###########################
#### Below is for the #####
### Supp table and the ####
### Manuscript box plot ###
###########################

## Fix the header of GCV_content.tsv for publication
## And add a "Genes per Mbp" column
## And convert some bp columns to Mbp
cd /home/liam/data
cat GC_content.tsv | \
awk '
NR==1 {print "Genome\tGC%\tGenome Size (Mbp)\tIntact Genes\tPseudogenes\tTotal Length Intact Genes (Mbp)\tCoding Density (%)\tGenes per Mbp\tMean Length Intact Genes (bp)"; next};
{ print $1"\t"$2"\t"$3/10^6"\t"$4"\t"$5"\t"$6/10^6"\t"$7"\t"$4/($3/10^6)"\t"$8 }
' \
> awk.tmp
mv awk.tmp GC_content.tsv

## Reformat data in GC_content.tsv for R to create a multi-pane chart. (This will be box plots)
## col1 = category; col2 = clade (Parent, Sister, T. azollae); col3 = NosAzo0708 (Y/N)
outfile=genome_stats_for_box_plots.tsv
echo -e "category\tclade\tvalue" > $outfile
## Loop through columns
n_col=$( cat GC_content.tsv | awk -F "\t" '{print NF}' | sort -u )
col=2
while (( $col <= $n_col )); do
cat GC_content.tsv | \
awk -v col=$col -F "\t" '
NR==1 { category=$col; next }
{ print category"\t"$1"\t"$col }
' \
>> $outfile
## end col loop
col=$(($col+1))
done

## Add col4 NosAzo0708 Y/N
cat $outfile | \
awk -F "\t" '
NR==1 {print $0"\tNosAzo0708"; next};
$2=="NosAzo0708" {print $0"\tY"; next}
{print $0"\tN"}
' \
> awk.tmp
mv awk.tmp $outfile
## Change prefix to clade
## T. az. MAG
cat ~/Taz.short_prefix.list | \
while read prefix; do
cat $outfile | \
awk -v prefix=$prefix -F "\t" ' $2==prefix { print $1"\tT. azollae\t"$3"\t"$4; next }; {print $0} ' > awk.tmp
mv awk.tmp $outfile
done
## Parent
cat ~/parent.short_prefix.list | \
while read prefix; do
cat $outfile | \
awk -v prefix=$prefix -F "\t" ' $2==prefix { print $1"\tParent\t"$3"\t"$4; next }; {print $0} ' > awk.tmp
mv awk.tmp $outfile
done
## Sister
cat ~/sister.short_prefix.list | \
while read prefix; do
cat $outfile | \
awk -v prefix=$prefix -F "\t" ' $2==prefix { print $1"\tSister\t"$3"\t"$4; next }; {print $0} ' > awk.tmp
mv awk.tmp $outfile
done
