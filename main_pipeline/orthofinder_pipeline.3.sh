##Run Orthofinder with outgroups!
##Rerppt the tree using gotree, and feed back to orthofinder.
##Then make .faa and .fna files of all the HOGs from the clade with all Nostocales + T. az., but no outgroups.
##Finally, identify the single copy N9.HOGs
##There are code blocks for ensuring that file manipulations perform as expected that require some manual work.

#############################
## Initial orthofinder run ##
#############################

## Run orthofinder using gloeobacter as the outgroup.
## Prepare input files for orthofinder!
## And rename input fastas to the shortened prefixes
conda deactivate
conda activate orthofinder
outdir=/home/liam/data/main_pipeline/orthofinder_out
mkdir $outdir
cd $outdir
indir=input_fastas
outdir=${outdir}/initial_run
mkdir $indir
pseudofinder_outdir=/home/liam/data/main_pipeline/pseudofinder_out
prokka_outdir=/home/liam/data/main_pipeline/prokka_out
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
## Get the MAGs (which don't have to be renamed)
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
initial_orthofinder_outdir=/home/liam/data/main_pipeline/orthofinder_out/initial_run/Results_Aug17
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
cd /home/liam/data/main_pipeline/orthofinder_out
mv initial_run/Results_Aug18 final_run
mv initial_run/Results_Aug17/* initial_run
rmdir initial_run/Results_Aug17 

## From here down, some initil code to prepare orthofinder outputs for downstream data analyses
## Get the amino acid and nucleotide fastas for the N9 HOGs
## and some other files that will be used downstream

## Originally intended to use N1, which would include all the Free-Living Nostocales and eliminate only the Gloeobacter outgroup.
## But this led to many instances of multiple loci from the same genome/MAG in the same N1.HOG
## Including nifH.
## This would make dN/dS - based analyses difficult.
## Using N9.HOGs splits the 3 copies of nifH from reference into 3 separate N9.HOGs. (double check this is true after re-run of pipeline)
## N9 removes only 7 Nostocales, in addition to the 2 Gloeobacter, whereas going any further woud eliminate the majority of Nostocales. 
## So all subsequent analyses will use the Phylogenetic_Hierarchical_Orthogrgroups/N9.tsv orthogroups.

##Check that loci are never repeated in a given OG, even if there are multiple HOGs.
orthofinder_outdir=/home/liam/data/main_pipeline/orthofinder_out/final_run
cd ${orthofinder_outdir}/Phylogenetic_Hierarchical_Orthogroups
for a in {0..9}; do
OG=OG000000${a}
echo $OG
tail -n +2 N9.tsv | awk -F "\t" -v OG=$OG '$2 == OG {for (i=4;i<=NF;++i) {print $i}}' | tr "," "\n" | sed 's/[[:space:]]//g' | awk '$0 !="" {print $0}' > tmp
#How many loci?
wc -l tmp
#How many unique loci?
sort -u tmp | wc -l
done
##All looks good! (for each OG, the two numbers match)
rm tmp

##Get the lists of loci in each HOG 
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
##Check that no loci from outside N9 are included in the locus lists.
cd /home/liam/data/main_pipeline/orthofinder_out/final_run/N9_HOGs_locus_lists
##Get list of genomes not in N9
grep -vFf ~/subtree.shortened.tsv ~/All.shortened.tsv > ~/outside_subtree.shortened.tsv
cat ~/outgroup.shortened.tsv ~/outside_subtree.shortened.tsv > tmp; mv tmp ~/outside_subtree.shortened.tsv 
##Check that no loci from outside N9 are included in the indir locus lists
> tmp
cat ~/outside_subtree.shortened.tsv | \
while read long_prefix short_prefix; do
echo $short_prefix
for file in *_loci.tsv; do grep "^${short_prefix}" $file >> tmp; done
done
sort -u tmp > tmp.tmp
##Looking at tmp.tmp, there are no loci! Just the genome names 
rm tmp tmp.tmp

##Make sure all the amino acid fasta inputs to orthofinder are unwrapped
indir=${orthofinder_outdir}/../input_fastas
for file in ${indir}/*.fasta; do
cat $file | \
awk ' {if ($0 ~ "^>") {printf "\n"$0"\n"} else { printf $0} }; END {printf "\n"} ' | \
tail -n +2 \
> tmp
mv tmp $file
done

##Get the amino acid sequences for each HOG
mkdir ${orthofinder_outdir}/aminoAcid_HOGs
cd ${orthofinder_outdir}/aminoAcid_HOGs
cat ~/subtree.shortened.tsv | \
while read long_prefix short_prefix; do
cp $indir/${short_prefix}.fasta ./
done

##Check that there is a prefix on every header of the fastas.
for file in *.fasta; do 
prefix=${file%.fasta}; 
echo $prefix >> tmp;
if (( $(grep -c "^>" $file) == $(grep "^>" $file | grep -c "$prefix") && $(grep -v "^>" $file | grep -c "$prefix") == 0 )); then
echo "looks good!" >> tmp
else
echo "looks bad!" >> tmp
fi 
done
grep -B 1 "looks bad!" tmp
grep -c "looks good!" tmp
wc -l tmp
rm tmp
## Looks good!

## Now make the amino acid fastas!
cat ${orthofinder_outdir}/N9_HOGs_locus_lists/HOGs2OGs.tsv | \
while read hog og; do
> ${hog}.faa
cat ~/subtree.shortened.tsv | \
while read long_prefix short_prefix; do
grep "^$short_prefix" ${orthofinder_outdir}/N9_HOGs_locus_lists/${hog}_loci.tsv | \
awk -F "\t" '{print $2}' | awk -F ", " '{for (i=1; i<=NF; ++i) {print $i}}' | \
while read locus; do
grep -A 1 ${locus} ${short_prefix}.fasta >> ${hog}.faa
done
done
done
rm *.fasta

## Get the nucleotide sequences for each HOG
mkdir ${orthofinder_outdir}/nucleotide_HOGs
cd ${orthofinder_outdir}/nucleotide_HOGs

## first, copy the nucleotide fastas and unwrap (they are already unwrapped from pseudofinder so not necessary.)
## and rename file with the short_prefix
cat ~/subtree.shortened.tsv | \
while read long_prefix short_prefix; do
file=${short_prefix}_intact.ffn
cp /home/liam/data/main_pipeline/pseudofinder_out/${long_prefix}/${long_prefix}_intact.ffn ./$file
## Unwrap
cat $file | \
awk ' {if ($0 ~ "^>") {printf "\n"$0"\n"} else { printf $0} }; END {printf "\n"} ' | \
tail -n +2 \
> tmp
mv tmp $file
done

##Check that every header contains the appropriate prefix
> tmp
for file in *_intact.ffn; do 
prefix=${file%_intact.ffn}
echo $prefix >> tmp;
if (( $(grep -c "^>" $file) == $(grep "^>" $file | grep -c "$prefix") && $(grep -v "^>" $file | grep -c "$prefix") == 0 )); then
echo "looks good!" >> tmp
else
echo "looks bad!" >> tmp
fi 
done 
grep -B 1 "looks bad!" tmp
grep -c "looks good!" tmp
wc -l tmp
rm tmp
##Looks good!


#Now make the HOG fastas!
cat ../N9_HOGs_locus_lists/HOGs2OGs.tsv | \
while read hog og; do
> ${hog}.ffn
cat ~/subtree.shortened.tsv | \
while read long_prefix short_prefix; do
grep "^$short_prefix" ../N9_HOGs_locus_lists/${hog}_loci.tsv | \
awk -F "\t" '{print $2}' | awk -F ", " '{for (i=1; i<=NF; ++i) {print $i}}' | \
while read locus; do
grep -A 1 ${locus} ${short_prefix}_intact.ffn >> ${hog}.ffn
done
done
done
rm *_intact.ffn


#Do some rearranging of orthofinder outputs that will be used downstream.
mkdir /home/liam/data/main_pipeline/orthofinder_out/for_downstream
cd /home/liam/data/main_pipeline/orthofinder_out/for_downstream
mv ${orthofinder_outdir}/aminoAcid_HOGs ./
mv ${orthofinder_outdir}/nucleotide_HOGs ./
mv ${orthofinder_outdir}/N9_HOGs_locus_lists ./locus_lists_HOGs

#Make a table with counts of how many loci from each genome are present in each N9.HOG
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

#Make a prefix list!
> ~/subtree.prefix.list
head -n 1 $infile > header.tmp
cat N9.cols.tmp | \
while read col; do
cat header.tmp | awk -F "\t" -v col="$col" '{print $col}' >> ~/subtree.prefix.list
done
rm *.tmp

##Now count how many loci from each genome are present in each N9.HOG
infile=N9_noOutgroups.tsv
head -n 1 $infile > header.tmp
tail -n +2 $infile | \
awk -F "\t" '{printf $1"\t"$2"\t"$3; for (i=4; i<=NF; ++i) {split($i,A,", "); printf "\t"length(A)}; print ""}'  \
> body.tmp
cat header.tmp body.tmp > n_intact.tsv
rm *.tmp

#Get a list of single copy N9.HOGs
tail -n +2 n_intact.tsv | \
awk -F "\t" '{ for (i=4; i<=NF; ++i) {if ($i != 1) {next} }; print $1}' \
> single_copy_genes.list
wc -l single_copy_genes.list
#There are 1042!
#Also, count the total number of N9.HOGs
tail -n +2 n_intact.tsv | wc -l
#18462

## Get a list of N9.HOGs that are present and predicted-intact in >=75% of the free-living (40*0.75=30) and T. az. (8*0.75) genomes
## To be used downstream.
cd /home/liam/data/main_pipeline/orthofinder_out/for_downstream

## Split n_intact.tsv into free living and T. az.
infile=n_intact.tsv
echo -e "\
Taz,/home/liam/Taz.shortened.tsv
free_living,/home/liam/N9_freeLiving.shortened.tsv\
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

## Get the list of hogs with at least 75% of both sets of genomes having exactly one locus
## And get the list of genomes in those sets.
echo "\
Taz
free_living\
" | \
while read file_prefix; do
file=${file_prefix}.n_intact.tsv

## for each hog, get the list of genomes with exactly one locus!
cat $file | \
awk -F "\t" '
NR==1 { for (i=1;i<=NF;++i) {A[i]= $i}; next };
{
	printf $1;
	for (i=2;i<=NF;++i) {if ($i==1) { printf "\t"A[i] } }
	print ""
}
' \
> ${file_prefix}.all_hogs_which_genomes_single_copy.tsv

## Now use the 0.75 cut-off to get the list of hogs
n=$(head -n 1 $file | awk '{print (NF-1)*0.75}')
cat ${file_prefix}.all_hogs_which_genomes_single_copy.tsv | \
awk -v n=$n -F "\t" 'NF >= n+1 {print $1}' \
> ${file_prefix}.hogs_passed_thresh.tmp

## Finish loop Taz vs free-living
done

## Get the hogs that passed both thresholds
cat Taz.hogs_passed_thresh.tmp free_living.hogs_passed_thresh.tmp | \
sort | \
uniq -c | \
awk '$1==2 {print $2}' \
> hogs_passed_75_thresh.list
rm Taz.hogs_passed_thresh.tmp free_living.hogs_passed_thresh.tmp

## Get all of the HOGs with at least 3 genomes from each of ingroup and outgroup
cat free_living.all_hogs_which_genomes_single_copy.tsv Taz.all_hogs_which_genomes_single_copy.tsv | \
awk -F "\t" 'NF>=4 {print $1}' | \
sort | \
uniq -c | \
awk '$1==2 {print $2}' \
> hogs_passed_3_thresh.list
