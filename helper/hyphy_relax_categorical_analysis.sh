## Look at the distribution of hyphy relax results across COG categories and KEGG_pathways

main_dir=./

## Get the portion of all results that are relax or intense
## Can be used by both COG and KEGG_pathway analyses
read percent_relax percent_intense <<< $(\
cat results.tsv | \
tail -n +2 | \
awk '$2!="-"' | \
awk '
BEGIN {c_r=0;c_i=0};
{
if ($5=="relaxation") {c_r+=1; next};
if ($5=="intensification") {c_i+=1};
};
END {print c_r/FNR"\t"c_i/FNR}
'\
)

## KEGG_pathay column has ko00000 and map00000
## I confirmed that these always perfectly match.
## For example: an annotation might be ko01234,ko22222,map01234,map22222
## So I can just use one or the other. I arbitrarily chose ko over map.
## Are map and ko annotations always exactly the same?
#cd ${hog_annotation_dir}
#> map_ko_disagree.list
#ls *_annotations.tsv | while read file; do
#echo $file; cat $file | \
#awk -F "\t" '{print $13}' | \
#tr "," "\n" > tmp;
#grep "map" tmp | sort -u | sed 's/^map//' > map.tmp;
#grep "ko" tmp | sort -u | sed 's/^ko//' > ko.tmp;
#if ! cmp ko.tmp map.tmp; then echo $file >> map_ko_disagree.list; fi;
#done
## yes!
#rm map_ko_disagree.list
#cd $main_dir
## The removing of the map values was done in the annotate_OGs.sh script

hog_annotation_file=/home/liam/data/main_pipeline/annotations/eggnog-mapper_out/hog_annotations/hog_annotations.tsv
## loop COG_category, KEGG_pathway, and KEGG_module
echo "COG_category
KEGG_pathway
KEGG_module" | \
sed 's/,/\t/g' | \
while read ontology; do


## Get the ontology categories for hogs that went through RELAX
echo "${ontology}" > this_col.tmp
cat results.tsv | tail -n +2 | awk '{print $1}' | \
while read hog; do
cat ${hog_annotation_file} | \
awk -v ontology=$ontology -v hog=$hog -F "\t" '
NR==1 { for (i=1;i<=NF;++i) { if ($i==ontology) {c=i; next} } };
$1==hog {print $c};
' \
>> this_col.tmp
done
paste results.tsv this_col.tmp > results_with_${ontology}.tsv
rm this_col.tmp

## To get KEGG pathway definitions, I copy-pasted from
## https://www.genome.jp/kegg/pathway.html
## into a plain text file (kegg_pathway_definitions.txt) which I then processed into a 2-col tsv
#cat ~/kegg_pathway_definitions.txt | \
#awk -F "" ' BEGIN {p=0}; $0~"[0-9][0-9][0-9][0-9][0-9]" {printf $0; p=1; next}; p==1 {print $0}; {p=0}' | \
#sed 's/ M    /\t/' | \
#sed 's/ M R    /\t/' | \
#sed 's/ M N    /\t/' | \
#sed 's/ N    /\t/' | \
#sed 's/    /\t/' | \
#sed 's/^/ko/' \
#> ~/KEGG_pathway_definitions.tsv
#rm ~/kegg_pathway_definitions.txt
## To get KEGG pathway definitions, I copy-pasted from
## https://www.genome.jp/brite/ko00002
## into a plain text file (kegg_module_definitions.txt) which I then processed into a 3-col tsv
## (3rd col is map which idk what it is but I kept it)
#cat ~/kegg_module_definitions.txt | \
#sed 's/[[:space:]]\+\[PATH:/\t/' | \
#sed 's/\([^]]\)$/\1\t-/' | \
#sed 's/]$//' | \
#sed 's/ map/,map/g' | \
#sed 's/[[:space:]]\+map/\tmap/' | \
#sed 's/[[:space:]]/\t/' | \
#> ~/KEGG_module_definitions.tsv
#rm ~/kegg_module_definitions.txt
## Nevermind the above!
## I got the module and pathway ko constituents and descriptions in a local R script "get_kegg_modules.R"
## Using package KEGGREST


## For each KEGG_pathway, count how many hogs with each RELAX result belong to that pathway.
## For hogs with multiple categories, they will be fully counted in each.
## All pathways have a 5 digit identifier, so ~ should only hit perfect matches
cat results_with_${ontology}.tsv | tail -n +2 | awk '{print $NF}' | tr "," "\n" | sort -u > category_col.tmp
> n_relax.tmp
> n_not_sig.tmp
> n_intense.tmp
> n_total_inclusive.tmp
cat category_col.tmp | \
while read category; do
cat results_with_${ontology}.tsv | \
awk -v category=$category '
BEGIN {c_r=0;c_n=0;c_i=0};
$6~category {
if ($5=="relaxation") {c_r+=1; next};
if ($5=="not_significant") {c_n+=1; next};
if ($5=="intensification") {c_i+=1};
};
END { print c_r >> "n_relax.tmp"; print c_n >> "n_not_sig.tmp"; print c_i >> "n_intense.tmp" }
'
cat $hog_annotation_file | \
awk -v category=$category -v ontology=$ontology -F "\t" '
BEGIN { running=0 };
NR==1 { for (i=1;i<=NF;++i) { if ($i==ontology) {c=i; next} } };
{ split($c,A,","); for (i=1;i<=length(A);++i) { if (A[i]==category) {running+=1; next} } };
END { print running };
' \
>> n_total_inclusive.tmp
done
paste category_col.tmp n_relax.tmp n_not_sig.tmp n_intense.tmp n_total_inclusive.tmp > category_result_counts.tmp
## Instead of the total number of HOGs assigned to this category, inclusive of both those that were and weren't passed to hyphy RELAX
## Just give the number that weren't passed.
cat category_result_counts.tmp | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"($5-$4-$3-$2)}' > awk.tmp
mv awk.tmp category_result_counts.tmp
rm category_col.tmp n_relax.tmp n_not_sig.tmp n_intense.tmp n_total_inclusive.tmp

## Get p-values for the results with regards to both intensification and relaxation

## Loop relaxation then intensification
echo "$percent_relax,relaxed
$percent_intense,intensified" | \
sed 's/,/\t/g' | \
while read percent ternary; do

## initiate output columns
## p_val
> p_val_col.tmp
## enriched or depleted (regardless of p_val)
> enriched_col.tmp

## Loop through categories
cat category_result_counts.tmp | \
awk '{print $2"\t"$3"\t"$4}' | \
while read n_relax n_not_sig n_intense; do

## Total HOGs in category
n=$( awk -v a=$n_relax -v b=$n_not_sig -v c=$n_intense 'BEGIN {print a+b+c}' )
if (( $n == 0 )); then echo "-" >> enriched_col.tmp; echo "-" >> p_val_col.tmp; continue; fi

## HOGs in relax or intense
if [[ $ternary == "relaxed" ]]; then x=$n_relax; elif [[ $ternary == "intensified" ]]; then x=$n_intense; else echo "ternary not an option!"; fi

## Get enriched/depleted
awk -v a=$percent -v x=$x -v n=$n '
BEGIN {
## Figure out if enriched or depleted
if ( (x/n)>a ) {
        enrichment="more"
} else if ( (x/n)<a ) {
        enrichment="less"
} else {
        enrichment="equal"
}
print enrichment
}
' \
>> enriched_col.tmp

## call R binom.test() and get p.value
Rscript --vanilla /home/liam/scripts/helper/binom_test_p_val.R x=$x n=$n p=$percent outfile=p_val_col.tmp

## Finish loop through categories
done

paste category_result_counts.tmp  enriched_col.tmp p_val_col.tmp > paste.tmp
rm p_val_col.tmp enriched_col.tmp

## Add header
echo "${ontology},n_relaxed,n_not_significant,n_intensified,n_not_analyzed,${ternary}_portion_relative_to_genome,p" | \
sed 's/,/\t/g' \
> ${ternary}_${ontology}_enrichment.tsv
cat paste.tmp >> ${ternary}_${ontology}_enrichment.tsv
rm paste.tmp

## Add ontology definitions
echo "${ontology}_def" > definitions.tmp
cat ${ternary}_${ontology}_enrichment.tsv | \
tail -n +2 | \
awk '{print $1}' | \
while read category; do
cat ~/${ontology}_definitions.tsv | \
awk -v category=$category -F "\t" 'BEGIN {s=0}; $1==category {print $2; s=1}; END { if (s==0) {print "-"} } ' \
>> definitions.tmp
done
paste ${ternary}_${ontology}_enrichment.tsv definitions.tmp > paste.tmp
mv paste.tmp ${ternary}_${ontology}_enrichment.tsv
rm definitions.tmp

## and sort by significance
head -n 1 ${ternary}_${ontology}_enrichment.tsv > ${ternary}_${ontology}_enrichment_sorted.tsv
cat ${ternary}_${ontology}_enrichment.tsv | \
tail -n +2 | \
sort -k7,7g \
>> ${ternary}_${ontology}_enrichment_sorted.tsv
mv ${ternary}_${ontology}_enrichment_sorted.tsv ${ternary}_${ontology}_enrichment.tsv

## Finish relax/intense loop
done
rm category_result_counts.tmp

## End loop KEGG_pathway / KEGG_module
done


###########################
### How many HOGs with ####
### KEGG Pathway do not ###
### have a RELAX result ###
###########################

## Determine the percent of HOGs that are intact in at least one T. az. that have a RELAX result, and break it down by pathway.

cd /home/liam/data/main_pipeline/hyphy_out/relax_out

## Get a list of HOGs intact in at least 1 T. az.
cat ~/data/main_pipeline/gene_loss/Taz.intact_vs_not.tsv | tail -n +2 | awk '{for (i=2;i<=NF;++i) {if ($i==1) {print $1; next}}}' > Taz_intact.tmp
wc -l Taz_intact.tmp ## 4032
## How many HOGs have RELAX results?
cat results.tsv | tail -n +2 | awk '$NF!="-" {print $1}' > HOGs_with_RELAX_result.tmp
wc -l HOGs_with_RELAX_result.tmp ## 2812
## So, 2812/4032 = 69.7% of T. az. HOGs have RELAX results.

## Get a list of HOGs with a KEGG pathway annotation (and the pathway)
cat /home/liam/data/main_pipeline/annotations/eggnog-mapper_out/hog_annotations/hog_annotations.tsv | tail -n +2 | awk -F "\t" '$5!="-" {print $1"\t"$5}' > hogs_w_pathway.tmp

## Get list of HOGs intact in at least T. az. with KEGG pathway (and the pathway)
grep -Ff Taz_intact.tmp hogs_w_pathway.tmp > Taz_hogs_w_pathway.tmp

## How many of those HOGs have relax results?
cat Taz_hogs_w_pathway.tmp | awk '{print $1}' > Taz_hogs_w_pathway_just_hog.tmp
wc -l Taz_hogs_w_pathway_just_hog.tmp ## 1070
cat results.tsv | tail -n +2 | awk '$NF!="-"' | grep -Ff Taz_hogs_w_pathway_just_hog.tmp | wc -l ## 982
## So, 982 / 1070 = 91.8% of HOGs with a KEGG pathway that are intact in at least 1 T. az. have a KEGG result

## How many HOGs do not have RELAX results from each KEGG pathway?
grep -vFf HOGs_with_RELAX_result.tmp Taz_hogs_w_pathway.tmp | awk '{print $2}' | tr "," "\n" | sort | uniq -c | sort -k1,1n | awk '{print $2"\t"$1}' > n_HOGs_Taz_pathway_missing.tsv
cat Taz_hogs_w_pathway.tmp | awk '{print $2}' | tr "," "\n" | sort | uniq -c | sort -k1,1n | awk '{print $2"\t"$1}' > n_HOGs_Taz_pathway_total.tsv
cat n_HOGs_Taz_pathway_missing.tsv n_HOGs_Taz_pathway_total.tsv | awk '{print $1}' | sort | uniq -c | awk '$1==2 {print $2}' > pathways.tmp
grep -Ff pathways.tmp n_HOGs_Taz_pathway_missing.tsv | sort -k1,1 > paste.tmp
grep -Ff pathways.tmp n_HOGs_Taz_pathway_total.tsv | sort -k1,1 > paste2.tmp
paste paste.tmp paste2.tmp | awk '$3!=$1' ## empty! good!
echo -e "pathway\tp_not_analyzed\tnot_analyzed\ttotal" > unanalyzed_hogs.tsv
paste paste.tmp paste2.tmp | awk '{print $1"\t"$2/$4"\t"$2"\t"$4}' | sort -k2,2g >> unanalyzed_hogs.tsv
rm paste.tmp paste2.tmp
wc -l pathways.tmp ## 103
wc -l n_HOGs_Taz_pathway_total.tsv ## 238
## So, 135/238=56.7% of pathways have 100% of HOGs with RELAX results
## 173/238=72.7% of pathways have 90+% of HOGs with RELAX results
## 212/238=89.1% of pathways have 75+% of HOGs with RELAX results
## Of the remainder, all but 2 pathways have <=6 HOGs. The two exceptions are (ko00650,0.3,3,10 - Butanoate metabolism) (ko01503,0.416667,5,12 - Cationic antimicrobial peptide (CAMP) resistance)

## The four "significant" pathways:
## ko00195,0.101695,6,59 - Photosynthesis
## ko03010,ko03010,0.0357143,2,56 - Ribosome
## ko01110,0.0512821,16,312 - Biosynthesis of secondary metabolites
## ko00970,0,25,25 - Aminoacyl-tRNA biosynthesis

rm Taz_intact.tmp HOGs_with_RELAX_result.tmp hogs_w_pathway.tmp Taz_hogs_w_pathway.tmp Taz_hogs_w_pathway_just_hog.tmp n_HOGs_Taz_pathway_missing.tsv n_HOGs_Taz_pathway_total.tsv pathways.tmp
