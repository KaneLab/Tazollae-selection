## This code is not meant to just be run
## but is more of a repository of useful code
## and a way to remember what I did.
## But it involves manual work, so the code alone is not useful.


################################
#### Look at RELAX results #####
### for HOGs with non-unique ###
###### gene_names or KOs #######
################################

cd ~/data/main_pipeline/annotations/eggnog-mapper_out/hog_annotations

## Get a list of HOGs intact in at least 1 T. az. genome
## This is just the Taz pan genome. I checked it matches the number in the Taz pan genome which doesn't already exist as a list.
cat ~/data/main_pipeline/super_table.tsv | tail -n +2 | awk -F "\t" ' $4!=0 {print $1} ' > Taz_hogs.tmp

## Get a list of gene_names that are only assigned to one HOG of the HOGs that are intact in the T. az. pan genome.
## And the same for KOs instead of gene_names
echo "gene_name
KO" | \
while read ontology; do
cat search_hog_by_${ontology}.tsv | \
awk -F "\t" '{split($2,A,","); for (i=1;i<=length(A);++i) {print $1"\t"A[i]}}' > gene_name_single_lines.tmp
> gene_hits.tmp
cat Taz_hogs.tmp | while read hog; do
echo $hog
cat gene_name_single_lines.tmp | \
awk -v hog=$hog -F "\t" '$2==hog {print $1}' \
>> gene_hits.tmp
done
cat gene_hits.tmp | sort | uniq -c | awk '$1!=1 {print $2}' > multicopy_Taz_${ontology}.list
rm gene_hits.tmp
rm gene_name_single_lines.tmp
## finish while read ontology loop
done
rm Taz_hogs.tmp
cat multicopy_Taz_KO.list multicopy_Taz_gene_name.list > ~/data/main_pipeline/multicopy_genes.tmp

## Look at the hyphy results for multicopy genes!
cd ~/data/main_pipeline
## Get KO, gene_name, and pathway annotations for all HOGs with hyphy relax results
cat super_table.tsv | tail -n +2 | awk -F "\t" '$11!="-" {print $1"\t"$2"\t"$3"\t"$11"\t"$14}' > tmp
grep -Ff multicopy_genes.tmp tmp | awk -F "\t" '{print $4}' | sort | uniq -c
#8 intensification
#140 not_significant
#385 relaxation
## in R: binom.test()

## How are the multicopy genes diestributed amongst the KEGG pathways (of those with relax results)?
grep -Ff multicopy_genes.tmp tmp | awk -F "\t" '{print $5}' | tr "," "\n" | sort | uniq -c | sort -k1,1n > smp
cat smp | while read num ko; do grep "^$ko" ~/KEGG_pathway_definitions.tsv >> smp.smp; done
echo "-" >> smp.smp
paste smp smp.smp | less -S
rm smp smp.smp





###########
###########

### For some of the pseudogenes, I want to confirm that they really do look like pseudogenes
### An example of something looking like it's not a pseudogene would be that all the loci for that HOG that are called pseudogenes are all "too short" but they are all the same length as each other.

## Get the HOGs thus:
less -S ~/data/main_pipeline/core_vs_shell_analyses/free_core_Taz_not_core_by_gene_name_and_KO.tsv
cd ~/data/main_pipeline/core_vs_shell_analyses/

## Get the pseudo gene locus tags for the pseudogene.
hog=N9.HOG0006391
cat ~/data/main_pipeline/assign_pseudogenes/pseudos_best_HOGs.tsv | \
awk -v hog=$hog '$4==hog {print $1}' \
> ~/data/main_pipeline/assign_pseudogenes/psbT.tmp

cd ~/data/main_pipeline/pseudofinder_out
## Look at the sequences of the pseudogenes.
cat ~/data/main_pipeline/assign_pseudogenes/psbT.tmp | \
while read short_locus; do
prefix=${short_locus%%_*};
echo $prefix;
long_prefix=$(grep $prefix ~/subtree.shortened.tsv | awk '{print $1}');
long_header=$(\
grep "^>$short_locus" ${long_prefix}/${long_prefix}_pseudos.shorter_headers.fasta | \
sed 's/^>//' | \
awk -F "[[:space:]]|_" '{print $4"_"$5"_pseudo_"$3}';\
)
grep -A1 $long_header ${long_prefix}/${long_prefix}_pseudos.fasta | \
tail -n 1
done

## Look at the .gff from pseudofinder that says why each was called a pseudo
cat ~/data/main_pipeline/assign_pseudogenes/psbT.tmp | while read short_locus; do prefix=${short_locus%%_*}; echo $prefix; long_prefix=$(grep $prefix ~/subtree.shortened.tsv | awk '{print $1}'); long_header=$(\
grep "^>$short_locus" ${long_prefix}/${long_prefix}_pseudos.shorter_headers.fasta | \
sed 's/^>//' | \
awk -F "[[:space:]]|_" '{print $4"_"$5"_pseudo_"$3}';\
); grep $long_header ${long_prefix}/${long_prefix}_pseudos.gff; done


## Look at length in nts of each pseudogene
cat ~/data/main_pipeline/assign_pseudogenes/psbT.tmp | while read short_locus; do prefix=${short_locus%%_*}; echo $prefix; long_prefix=$(grep $prefix ~/subtree.shortened.tsv | awk '{print $1}'); long_header=$(\
grep "^>$short_locus" ${long_prefix}/${long_prefix}_pseudos.shorter_headers.fasta | \
sed 's/^>//' | \
awk -F "[[:space:]]|_" '{print $4"_"$5"_pseudo_"$3}';\
); grep -A1 $long_header ${long_prefix}/${long_prefix}_pseudos.fasta | tail -n 1| awk '{print length($0)}'; done


##########
##########


## look at HOGs present in NosAzo0708 but not in any other T. az.
## Looking at ~/data/GC_content.tsv, NosAzo0708 has 337-462 more intact ORFs than the other T. az.
## First, look after the improve_gene_inventory runs (so end of contig are included as "intact")
cd ~/data/main_pipeline/gene_loss
> tmp
cat Taz.intact_vs_not.tsv | \
awk -F "\t" ' $2==1 { for (i=3;i<=NF;++i) { if ($i==1) {next} }; print $1 }' | \
while read hog; do
grep $hog ~/data/main_pipeline/annotations/eggnog-mapper_out/hog_annotations/hog_annotations.tsv >> tmp
done
rm tmp
## looking at tmp, there are only 6 hogs, none of which have super useful annotations except an endonuclease and a histidine kinase
## Now, look using the original orthofinder outputs for HOGs present in NosAzo0708 and not in any other T. az.
cd ~/data/main_pipeline/orthofinder_out/for_downstream
> tmp
cat Taz.n_intact.tsv | \
tail -n +2 | \
awk -F "\t" ' $2>0 { for (i=3;i<=NF;++i) { if ($i!=0) {next} }; print $1 } ' | \
while read hog; do
grep $hog ~/data/main_pipeline/annotations/eggnog-mapper_out/hog_annotations/hog_annotations.tsv >> tmp
done
## Of these 25 HOGs: 5 transposase, 9 very poor annotation, 3 photosynthesis
rm tmp
## How about all HOGs for which NosAzo0708 has the most intact copies of the T. az.?
cat Taz.n_intact.tsv | \
tail -n +2 | \
awk -F "\t" ' { a=$2; for (i=3;i<=NF;++i) { if ($i>=a) {next} }; print $0 } ' | \
awk 'NR==1 {for (i=2;i<=NF;++i) {A[i-1]=0}}; {for (i=2;i<=NF;++i) {A[i-1]=A[i-1]+$i}}; END {for (i=1;i<=length(A);++i) {print A[i]}}' | \
sort -g
## These 52 HOGs account for 266-295 extra intact ORFs for NosAzo0708
> tmp
cat Taz.n_intact.tsv | \
tail -n +2 | \
awk -F "\t" ' { a=$2; for (i=3;i<=NF;++i) { if ($i>=a) {next} }; print $1 } ' | \
while read hog; do
grep $hog ~/data/main_pipeline/annotations/eggnog-mapper_out/hog_annotations/hog_annotations.tsv >> tmp
done
## How many completely un-annotated?
cat tmp | awk '{for (i=2;i<=NF;++i) {if ($i!="-") {next}}; print $0}' | wc -l #19 = 15 + 4 more that have very minimal annotation
## How many nuclease or transposase?
cat tmp | grep -i "transposase\|nuclease" | wc -l #21 = 19 + 2 more ("transposition" and "chromosome replication")
## Of remaining 12: 4 definitely photosynthesis involved.
rm tmp

## look at annotations of the free-core that are shell or absent from T. az.
## intact/Taz_absent_free_core.list and intact/Taz_shell_free_core.list
## For manual inspection:
cd ~/data/main_pipeline/core_vs_shell_analyses/intact
> missing_core_annotations.tsv
cat Taz_absent_free_core.list Taz_shell_free_core.list | \
while read hog; do
echo $hog >> missing_core_annotations.tsv;
cat ~/data/main_pipeline/annotations/eggnog-mapper_out/hog_annotations/${hog}_annotations.tsv | \
awk -F "\t" '{ for (i=7;i<=14;++i) {printf $i"\t"}; print $21}' \
>> missing_core_annotations.tsv;
done



## Looking for specific photosynthesis-related annotated hogs:
cd ~/data/main_pipeline/annotations/eggnog-mapper_out/hog_annotations
> ccm.tmp
grep -i "ccm\|carbon dioxide\|carboxysome\|co2\|ribulose bisphosphate\|rubisco" * | \
sed 's/:.*//' | \
sort -u | \
grep "tsv$" | \
while read file; do
echo $file >> ccm.tmp;
cat $file >> ccm.tmp;
done
## Manually inspect ccm.tmp for hogs that are annotated as ccms
## Make some lists of hogs (i.e. rubisco.list) in manual_hog_annotations/
# Nope (but some look interesting!): N9.HOG0004422 N9.HOG0005050 N9.HOG0009300 N9.HOG0012804 N9.HOG0014332 N9.HOG0014509 N9.HOG0010597
## Look at how many T. az. loci in each:
> ccm.tmp
echo "N9.HOG0004481 (L) N9.HOG0004483 (S) N9.HOG0005223 (S - also annotated as ccmM.) N9.HOG0006302 (S - confirmed with ncbi blast bc annotation poor)" | \
tr " " "\n" | \
grep "^N9" | \
while read hog; do
echo $hog >> ccm.tmp;
cat ${hog}_annotations.tsv | \
awk '{print $1}' | \
#grep -Ff ~/Taz.short_prefix.list | \
sed 's/_.*//' | \
sort | \
uniq -c \
>> ccm.tmp;
echo "" >> ccm.tmp;
done
## Better to look in the binary tables that take the end-of-contig-pseudos and re-assembly results into account.


################################
### Looking at RELAX results ###
################################

## Look at KEGG analysis for different HOG sets
## including: only core HOGs and only HOGs that have unique names or KOs in annotations
## Can use T. az. core, free-living core, and the intersecting core
## I think the free-living core is the most interesting.
## Note: these are not necessarily single-copy in each genome, so this doesn't mean
## that every genome was included in the RELAX call for this HOG.
## Loop core definitions:
## Note either core/absent combo is irrelevant because none of thse could be run through RELAX
## Create a list of all the free-living core genes and of all the T. az. core genes
core_list_dir=~/data/main_pipeline/core_vs_shell_analyses/HOG/intact
unique_annotation_list_dir=~/data/main_pipeline/annotations/eggnog-mapper_out/hog_annotations
cat ${core_list_dir}/Taz_*_free_core.list > ${core_list_dir}/free_core.list
cat ${core_list_dir}/Taz_core_free_*.list > ${core_list_dir}/Taz_core.list
## Make dir for this analysis
mkdir ~/data/main_pipeline/hyphy_out/relax_out/core_only
echo "${core_list_dir}/Taz_core_free_core.list
${core_list_dir}/free_core.list
${core_list_dir}/Taz_core.list
${unique_annotation_list_dir}/HOGs_w_unique_ko.list" | \
while read hog_list; do
hog_set=${hog_list##*/}; hog_set=${hog_set%.list}
cd ~/data/main_pipeline/hyphy_out/relax_out/core_only
mkdir $hog_set
cd $hog_set
## Want to recalculate FDR q-values as set is smaller.
head -n 1 ~/data/main_pipeline/hyphy_out/relax_out/results_no_FDR.tsv > results_no_FDR.tsv
tail -n +2 ~/data/main_pipeline/hyphy_out/relax_out/results_no_FDR.tsv | \
grep -Ff $hog_list \
>> results_no_FDR.tsv
## Of 1171 core hogs, 1 does not have a RELAX result: N9.HOG0002332 (psbD)
## It appears to have 2 predicted-intact loci in many genomes, which would therefore not qualify for RELAX
## RE-calculate p-values! (Stolen from hyphy_relax.9.sh)
## calculate FDR q-values!
## add a column to the header!
head -n 1 results_no_FDR.tsv | \
awk '{print $1"\t"$2"\t"$3"\tq-val"}' \
> results.tsv

## Calculate q-vals
m=$(( $(wc -l < results_no_FDR.tsv) - 1 ))
cat results_no_FDR.tsv | \
tail -n +2 | \
awk '$2!="-"' | \
sort -k3,3g \
> tmp
m=$(wc -l < tmp)
cat tmp | \
awk -v m=$m '{print $1"\t"$2"\t"$3"\t"$3*m/NR}' \
>> results.tsv
rm tmp

## add back the lines for hogs that did not converge
cat results_no_FDR.tsv | \
tail -n +2 | \
awk '$2=="-" {print $0}' \
>> results.tsv

## Give a trinary result using q value.
q=0.1 ## arbitrary. 0.1 is used by Jones et al, 2023, Nat. Eco. Evo.
echo "result" > trinary_results.tmp
tail -n +2 results.tsv | \
awk -v q=$q '
$4=="-" {print "-"; next};
$4>q {print "not_significant"; next};
$2<1 {print "relaxation"; next};
$2>1 {print "intensification"; next};
$2==1 {print "equal"; next};
' \
>> trinary_results.tmp
## Did any come back equal? Hopefully 0
awk '$4=="equal"' trinary_results.tmp | wc -l
## Paste
paste results.tsv trinary_results.tmp > paste.tmp
mv paste.tmp results.tsv
rm trinary_results.tmp

## COG and KEGG analysis!
bash -i ~/scripts/helper/hyphy_relax_categorical_analysis.sh

## Also curious how many of each trinary result from the full analysis are core genes
## And using just this subset
## They should be the same unless the FDR affected things.
echo "full_analysis_trinary_percents.tsv,~/data/main_pipeline/hyphy_out/relax_out/results.tsv
this_set_only_trinary_percents.tsv,results.tsv" | \
sed 's/,/\t/g' | \
while read outfile infile; do
> $outfile
tail -n +2 $infile | \
grep -Ff $hog_list \
> table.tmp
n_r=$(cat table.tmp | awk '$NF=="relaxation"' | wc -l)
n_i=$(cat table.tmp | awk '$NF=="intensification"' | wc -l)
n_n_s=$(cat table.tmp | awk '$NF=="not_significant"' | wc -l)
n_t=$(( $n_r + $n_i + $n_n_s ))
echo "relaxation,$n_r
not_significant,$n_n_s
intensification,$n_i" | \
sed 's/,/\t/g' | \
while read category n; do
awk -v category=$category -v n=$n -v t=$n_t 'BEGIN {print category"\t"n"\t"n/t}' >> $outfile
done
## Finish outfile infile loop
rm table.tmp
done

## Finish hog_set loop
done

