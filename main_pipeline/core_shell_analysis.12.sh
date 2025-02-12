## Generate binary tables to allow annotation-based core/shell analysis
#bash -i inventory_by_annotation.sh

#########################
## Core/Shell analyses ##
#########################

main_dir=/home/liam/data/main_pipeline/core_vs_shell_analyses
mkdir $main_dir
cd $main_dir

## Loop through indirs
## binary tables by HOG, by KO, or by gene name
echo "/home/liam/data/main_pipeline/gene_loss,HOG
/home/liam/data/main_pipeline/inventory_by_annotation/KO,KO
/home/liam/data/main_pipeline/inventory_by_annotation/gene_name,gene_name" | \
sed 's/,/\t/g' | \
while read indir category; do
echo $category
outdir=${main_dir}/${category}
mkdir $outdir
cd $outdir
## sum the rows in the binary tables
## and output a table with each row a hog and each column the sum (n)
## or the sum over the number of genomes (p)
outfile=genome_counts_per_hog.tsv
## The Taz and freeliving binary tables don't contain all the HOGs if the given subset would be all 0s
## So start with all the HOGs.
## This is not true for KO and gene_name
if [[ $category == "HOG" ]]; then
cat ${indir}/all.intact_vs_not.tsv | awk '{print $1}' > $outfile
else
cat ${indir}/Taz.intact_vs_not.tsv | awk '{print $1}' > $outfile
fi
tail -n +2 $outfile > hogs.tmp
## Loop through the 4 binary tables
echo "\
${indir}/free_living.intact_vs_not.tsv
${indir}/Taz.intact_vs_not.tsv
${indir}/free_living.presence_absence.tsv
${indir}/Taz.presence_absence.tsv\
" | \
while read table; do
subset=${table%.tsv}
subset=${subset##*/}
#Get the count of 1s in this row as well as the percent.
echo -e "${category}\tn_${subset}\tp_${subset}" > col.tmp
tail -n +2 $table | \
awk -F "\t" '
{ c=0; for (i=2;i<=NF;++i) { c+=$i }; printf $1"\t"c"\t%.3f\n", c/(NF-1) }
' >> col.tmp
## Print 0\t0 for all the hogs that are missing from this binary table
cp hogs.tmp unmatched_hogs.tmp
cat $table | awk '{print $1}' | tail -n +2 >> unmatched_hogs.tmp
cat unmatched_hogs.tmp | sort | uniq -c | awk '$1==1 {print $2"\t"0"\t"0}' >> col.tmp
## sort to keep the order consistent
head -n 1 col.tmp > sort.tmp
cat col.tmp | tail -n +2 | sort -k1,1 >> sort.tmp
mv sort.tmp col.tmp

## add col to the outfile table
paste $outfile col.tmp > tmp
mv tmp $outfile
rm col.tmp unmatched_hogs.tmp

## finish loop through tables
done
rm hogs.tmp

## Check that all the cols are in the same order
echo "Any lines out of order? (should be 0):"
cat $outfile | awk ' $1!=$2 || $1!=$5 || $1!=$8 || $1!=$11 ' | wc -l
## Remove hog cols
cat $outfile | awk ' { printf $1; for (i=3;i<=NF;++i) { if (i%3!=2) {printf "\t"$i} }; print "" } ' > out.tmp
mv out.tmp $outfile



## List core, shell, and absent HOGs for T. az. and for Nostocales

## Initiate the output table
outfile=core_shell_absent.tsv
cat genome_counts_per_hog.tsv | \
head -n 1 | \
awk ' { printf $1; for (i=2;i<=NF;++i) { if (i%2==1) {printf "\t"$i} }; print ""; } ' | \
sed 's/\tp_/\t/g' \
> $outfile
## Loop through the column of interest
cat genome_counts_per_hog.tsv | \
tail -n +2 | \
awk '
{
printf $1;
for (i=2;i<=NF;++i) {
	printf "\t";
	if (i%2==1) {
		if ($i==0) {printf "absent"}
		else if ($i==1) {printf "core"}
		else {printf "shell"}
	}
};
print "";
}
' \
>> $outfile

## Get a list of hogs for each combination of shell/core/absent for free/Taz
infile=$outfile
echo "2,intact
4,present" | \
sed 's/,/\t/g' | \
while read free_col intact_present; do
mkdir $intact_present
echo "core,shell,absent" | \
sed 's/,/\n/g' | \
while read Taz_cat; do
echo "core,shell,absent" | \
sed 's/,/\n/g' | \
while read free_cat; do
cat $infile | \
tail -n +2 | \
awk -v free_col=$free_col -v Taz_cat=$Taz_cat -v free_cat=$free_cat '
$free_col==free_cat && $(free_col+1)==Taz_cat {print $1}
' \
> ${intact_present}/Taz_${Taz_cat}_free_${free_cat}.list
done
done
done


## Turn this into 2, 3x3 tables of core vs. absent vs. shell for Nostoc and T. az. for present or intact
echo "intact
present" | \
while read intact_present; do
outfile=${intact_present}_venn.tsv
echo -e "Taz\\Free\ncore\nshell\nabsent" > $outfile
echo "core,shell,absent" | \
sed 's/,/\n/g' | \
while read free_cat; do
## How many of each Taz_cat in this free_cat?
Taz_core=$(wc -l < ${intact_present}/Taz_core_free_${free_cat}.list)
Taz_shell=$(wc -l < ${intact_present}/Taz_shell_free_${free_cat}.list)
Taz_absent=$(wc -l < ${intact_present}/Taz_absent_free_${free_cat}.list)
## Turn into a column with percents
awk -v free_cat=$free_cat -v core=$Taz_core -v shell=$Taz_shell -v absent=$Taz_absent '
BEGIN {
total=core+shell+absent
print free_cat" ("total")";
printf core" (%.2f)\n", core/total;
printf shell" (%.2f)\n", shell/total;
printf absent" (%.2f)\n", absent/total;
}
' \
> new_col.tmp
paste $outfile new_col.tmp > out.tmp
mv out.tmp $outfile
done
done
rm new_col.tmp

## Finish indir loop
done

###############
### Figures ###
###############

## Make a table for an upset plot in R of intact HOGs for each T. az., free-livng, and combined (all)

cd ~/data/main_pipeline/core_vs_shell_analyses/HOG

## Make a binary table of core HOGs
#echo "Free-Living Core" > core_binary.tsv
#cat core_shell_absent.tsv | tail -n +2 | awk '{print $2}' | sed 's/shell/0/' | sed 's/core/1/' | sed 's/absent/0/' >> core_binary.tsv
## Paste onto T. az. intact binary table
#paste /home/liam/data/main_pipeline/gene_loss/Taz.intact_vs_not.tsv core_binary.tsv > venn_binary_Taz_intact_free_core.tsv

## Loop Taz, free-living, and combined (all)
echo "Taz
free_living
all" | \
while read group; do

## Reformat to be first col HOG and second col a tsv of which genomes/sets the HOG is in
outfile=upset_input_${group}_intact.tsv
echo -e "HOG\tGenomes" > ${outfile}
cat /home/liam/data/main_pipeline/gene_loss/${group}.intact_vs_not.tsv | \
awk -F "\t" '
NR==1 { for (i=1;i<=NF;++i) { A[i]=$i }; next }
{
printf $1"\t";
for (i=2;i<=NF;++i) {
	if ($i==1) {
		printf A[i]","
	}
};
print ""
}
' | \
sed 's/,$//' | \
awk -F "\t" '$2!=""' \
>> ${outfile}

## Do counts of intersections:
echo -e "Genomes\tCounts" > tmp
cat $outfile | \
tail -n +2 | \
awk -F "\t" '{print $2}' | \
sort | \
uniq -c | \
sort -k1,1nr | \
sed 's/^[[:space:]]\+//' | \
awk '{print $2"\t"$1}' \
>> tmp
mv tmp $outfile
## Upset plot using upset_input_Taz_intact.tsv
## Local script Taz_upset.R which uses UpSetR

## finish loop while read group; do
done


#############################
### Size of core and pan ####
### genome from each node ###
#############################

cd /home/liam/data/main_pipeline/core_vs_shell_analyses/HOG
binary_infile=/home/liam/data/main_pipeline/gene_loss/all.intact_vs_not.tsv
in_tree=/home/liam/data/main_pipeline/hyphy_out/concatenated_msas/single_copy_HOGs.pruned_tree.txt
core_outfile=core_size_each_node.tsv
pan_outfile=pan_size_each_node.tsv

## Get the number of HOGs in each genome
cat $binary_infile | \
awk '
NR==1 { for (i=1;i<=NF;++i) { A[i]=$i; C[i]=0 } }
{ for (i=2;i<=NF;++i) { C[i]+=$i } };
END { for (i=2;i<=length(A);++i) { print A[i]"\t"C[i] } }
' \
> $core_outfile
cp $core_outfile $pan_outfile

## Loop through internal nodes
conda activate gotree
gotree labels --internal=true --tips=false -i $in_tree | \
while read internal_node; do
echo $internal_node

> binary_subset.tsv

## Get all leaves in the clade with that as the root
gotree subtree -n "^${internal_node}$" -i $in_tree | \
gotree labels | \
while read leaf; do
## Subset binary table to leaves in this clade
cat $binary_infile | \
awk -v leaf=$leaf -F "\t" '
NR==1 { for (i=1;i<=NF;++i) { if ($i==leaf) {c=i; next} } };
{print $c}
' \
> this_leaf.tmp
paste binary_subset.tsv this_leaf.tmp > paste.tmp
mv paste.tmp binary_subset.tsv
rm this_leaf.tmp
## Finish loop while read leaf; do
done
sed -i 's/^\t//' binary_subset.tsv
## Count HOGs with positive for all leaves in this clade!
cat binary_subset.tsv | \
awk -v internal_node=$internal_node '
BEGIN {c=0};
{ for (i=1;i<=NF;++i) { if ($i==0) {next} }; c+=1 };
END {print internal_node"\t"c}
' \
>> $core_outfile
## Count HOGs with positive for at least 1 leaf in this clade!
cat binary_subset.tsv | \
awk -v internal_node=$internal_node '
BEGIN {c=0};
{ for (i=1;i<=NF;++i) { if ($i==1) {c+=1; next} } };
END {print internal_node"\t"c}
' \
>> $pan_outfile
## Finish loop while read internal_node; do
rm binary_subset.tsv
done
conda deactivate

## Make iTOL data sheets
cat ~/iTol_metadata_template.txt $core_outfile > iTOL_$core_outfile
cat ~/iTol_metadata_template.txt $pan_outfile > iTOL_$pan_outfile


########################
### Ontology of core ###
### and pan genomes ####
########################

cd /home/liam/data/main_pipeline/core_vs_shell_analyses/HOG/intact
## List the core and shell for free and Taz
cat Taz_core_free_*.list > Taz_core.list
cat Taz_shell_free_*.list > Taz_shell.list
cat Taz_*_free_core.list > free_core.list
cat Taz_*_free_shell.list > free_shell.list

## Loop ontologies
echo "KEGG_pathway,5
COG_category,6" | \
sed 's/,/\t/g' | \
while read ontology ontology_col; do

## Count how many HOGs belong to each ontology category for each core/shell
## For "-" (Unassigned), make it a percent of total HOGs in that core/shell
## For rest, make it a percent of assigned HOGs in that core/shell
annotation_file=/home/liam/data/main_pipeline/annotations/eggnog-mapper_out/hog_annotations/hog_annotations.tsv
echo "Taz_core,Taz_shell,free_core,free_shell" | \
sed 's/,/\n/g' | \
while read core_shell; do
## How many HOGs in this category?
n_HOGs=$(wc -l < ${core_shell}.list)
grep -Ff ${core_shell}.list $annotation_file | \
awk -v ontology_col=$ontology_col -F "\t" '{print $ontology_col}' | \
tr "," "\n" | \
sort | uniq -c > assigned_or_not.tmp
## Get the number of unassigned
n_unassigned=$( cat assigned_or_not.tmp | awk '$2=="-" {print $1; exit}' )
cat assigned_or_not.tmp | \
awk -v n_HOGs=$n_HOGs -v n_unassigned=$n_unassigned '
$2=="-" {print $2"\t"($1/n_HOGs); next}
{print $2"\t"($1/(n_HOGs-n_unassigned))}
' | \
sort -k2,2g \
> ${core_shell}_ontology_categories.tsv
rm assigned_or_not.tmp
## Finish loop while read core_shell
done

## Combine into a single table
## Get every pathway that was hit at least once
cat Taz_core_ontology_categories.tsv Taz_shell_ontology_categories.tsv free_core_ontology_categories.tsv free_shell_ontology_categories.tsv | \
awk '{print $1}' | \
sort -u \
> pathways.col
## Loop through pathways
rm Taz_core.col Taz_shell.col free_core.col free_shell.col definitions.col
cat pathways.col | \
while read pathway; do

## Reorder each core_shell_ontology_categories.tsv to be in the same order by pathway
## Loop through core_shell
echo "Taz_core,Taz_shell,free_core,free_shell" | \
sed 's/,/\n/g' | \
while read core_shell; do
a=$( cat ${core_shell}_ontology_categories.tsv | awk -v pathway=$pathway '$1==pathway {print $2; exit}' )
## if pathway did not appear in ${core_shell}_ontology_categories.tsv, then print "0"
## Otherwise, print the value
if [[ $a ]]; then echo $a >> ${core_shell}.col; else echo "0" >> ${core_shell}.col; fi
## Finish loop while read core_shell
done

## Get the ontology definitions
a=$( cat ~/${ontology}_definitions.tsv | awk -v pathway=$pathway -F "\t" '$1==pathway {print $2; exit}' )
if [[ $a ]]; then echo $a >> definitions.col; else echo "-" >> definitions.col; fi

## Finish loop while read pathway
done

## paste cols together and sort by percent of Taz_core
echo -e "${ontology}\tTaz_core\tfree_core\tTaz_shell\tfree_shell\tPathway_def" > core_shell_${ontology}.tsv
paste pathways.col Taz_core.col free_core.col Taz_shell.col free_shell.col definitions.col | \
sort -k2,2gr \
>> core_shell_${ontology}.tsv
## manually inspect compared to earlier files to make sure they are correct

## Finish loop while read ontology ontology_col
done

## Clean up
rm Taz_core_ontology_categories.tsv Taz_shell_ontology_categories.tsv free_core_ontology_categories.tsv free_shell_ontology_categories.tsv \
pathways.col Taz_core.col Taz_shell.col free_core.col free_shell.col definitions.col Taz_core.list Taz_shell.list free_core.list free_shell.list

## Make long form for a grouped bar plot
## Subset to only pathways that are at least 5% in one of the 4 sets
echo -e "Pathway\tPan Genome Subset\tPercentage of HOGs" > core_shell_KEGG_pathway.long.tsv
cat core_shell_KEGG_pathway.tsv | \
awk -F "\t" '
NR==1 { for (i=1;i<=NF;++i) {A[i]=$i}; next };
$2 >= 0.05 || $3 >= 0.05 || $4 >= 0.05 || $5 >= 0.05 {
print $6"\t"A[2]"\t"$2
print $6"\t"A[3]"\t"$3
print $6"\t"A[4]"\t"$4
print $6"\t"A[5]"\t"$5
}
' | \
sed 's/^-/Unassigned/' \
>> core_shell_KEGG_pathway.long.tsv



######################################
### Deeper look at core Nostocales ###
####### HOGs missing in T. az. #######
######################################

## Get a list of HOGs that are core to the free-living, but not to T. az.
cat Taz_shell_free_core.list Taz_absent_free_core.list > free_core_Taz_not_core.list
#66 HOGs

## Make a .dmnd database of the intact protein sequences for each T. az.
conda deactivate
conda activate orthofinder
indir=/home/liam/data/main_pipeline/orthofinder_out/input_fastas
cat ~/Taz.short_prefix.list | \
while read prefix; do
#diamond makedb --in ${indir}/${prefix}.fasta -d ${prefix}.dmnd
makeblastdb -dbtype "prot" -in ${indir}/${prefix}.fasta -out ${prefix}.blastdb
done

## Diamond blast these sequences against the T. az. to check for paralogs
## Create a header for output
indir=/home/liam/data/main_pipeline/orthofinder_out/for_downstream/aminoAcid_HOGs
#echo "HOG" > blast_hit_binary.tsv
#cat ~/Taz.short_prefix.list >> blast_hit_binary.tsv
#cat blast_hit_binary.tsv | tr "\n" "\t" | sed 's/\t$//' | awk '{print $0}' > tr.tmp
#mv tr.tmp blast_hit_binary.tsv

## Loop through T. az. dmnd dbs
cat ~/Taz.short_prefix.list | \
while read prefix; do
echo $prefix > ${prefix}.col

## Loop through hogs
cat free_core_Taz_not_core.list | \
while read hog; do

## diamond blast the HOG of interest against the T. az. intact protein sequences of interest
#diamond blastp --threads 20 --outfmt 6 "length nident qlen slen evalue" --query ${indir}/${hog}.faa --db ${prefix}.dmnd --out dmnd_out.tsv
## blast the HOG of interest against the T. az. intact protein sequences of interest
blastp -outfmt "6 length nident qlen slen evalue" -query ${indir}/${hog}.faa -db ${prefix}.blastdb -out blast_out.tsv
## Filter to evalue and pident and percent coverage. Use 0.65 percent coverage because that is the default in pseudofinder.
a=$( cat blast_out.tsv | awk ' $5 < 1e-10 && $1/$3 > 0.65 && $1/$3 < 1.65 && $1/$4 > 0.65 && $1/$4 < 1.65 ' | wc -l )
if (( $a > 0)); then echo "1" >> ${prefix}.col; else echo "0" >> ${prefix}.col; fi

## Finish loop while read hog
done

## Finish loop while read prefix
done

## Clean up and paste
rm blast_out.tsv *blastdb*
echo "HOG" > blast_hit_binary.tsv
cat free_core_Taz_not_core.list >> blast_hit_binary.tsv
paste blast_hit_binary.tsv NosAzo0708.col Caroliniana3004.col Caroliniana3017.col FiliculoidesGal.col Mexicana2001.col Microphylla4021.col Nilotica5001.col Rubra6502.col > paste.tmp
mv paste.tmp blast_hit_binary.tsv
rm NosAzo0708.col Caroliniana3004.col Caroliniana3017.col FiliculoidesGal.col Mexicana2001.col Microphylla4021.col Nilotica5001.col Rubra6502.col
conda deactivate

## Get the HOGs that do not have a hit for each T. az.
cat blast_hit_binary.tsv | awk 'NR==1 {print $0; ext}; {for (i=2;i<=NF;++i) { if ($i==0) {print $0; next} } }' > free_core_Taz_not_core_after_blast.tsv
cat free_core_Taz_not_core_after_blast.tsv | awk '{print $1}' | tail -n +2 > free_core_Taz_not_core_after_blast.list



############################
### Look at missing core ###
############################

## Look at the free_core that are not Taz_core

super_table=/home/liam/data/main_pipeline/super_table.tsv

## First, using HOG:
## How many such HOGs?
wc -l /home/liam/data/main_pipeline/core_vs_shell_analyses/HOG/intact/free_core_Taz_not_core.list #66

## How much overlap with the ontology names that are core to all free-living but not T. az.?
grep -Ff /home/liam/data/main_pipeline/core_vs_shell_analyses/HOG/intact/free_core_Taz_not_core.list $super_table | awk '{print $2}' | awk '$0!="-"' > names_of_HOGs.tmp
wc -l names_of_HOGs.tmp ##17 so 66-17=49 no name
grep -vFf /home/liam/data/main_pipeline/core_vs_shell_analyses/gene_name/intact/free_core_Taz_not_core.list names_of_HOGs.tmp
## mutT
grep "mutT" $super_table
## 2 HOGs. Both under relaxed selection. Both core for free-living. One core for T. az., one intact in 6 and pseudogenized in 2 (Mexicana2001_p_02210, Microphylla4021_p_02015)
rm names_of_HOGs.tmp

## How much overlap with the ko #s that are core to all free-living but not T. az.?
grep -Ff /home/liam/data/main_pipeline/core_vs_shell_analyses/HOG/intact/free_core_Taz_not_core.list $super_table | awk '{print $3}' | awk '$0!="-"' > names_of_HOGs.tmp
wc -l names_of_HOGs.tmp #28 so 66-28=38 no KO
cat names_of_HOGs.tmp | tr "," "\n" | \
grep -vFf /home/liam/data/main_pipeline/core_vs_shell_analyses/KO/intact/free_core_Taz_not_core.list
#ko:K08696 has 3 HOGs. ccmK1/2. Everybody's got at least 2, except FiliculoidesGal has 1.
#ko:K04077 has 5 HOGs. groL. Everybody's got one.
#ko:K00858 has 2 HOGs. nadK1 and nadK2. Both core-core except nadK1 pseudogenized in Nilotica5001
#ko:K01649 has 6 HOGs. leuA. Everybody has at least 1
#ko:K03574 has 5 HOGs. mutT. Everybody has at least 1
#ko:K03564 has 4 HOGs. 37 free-living hav all 4. Rest have 3. T. az. all have the same 2.
#ko:K01563 has six HOGs. copies present in all genomes
#ko:K07259 has two HOGs. Both core to free-living. One core to T. az., one absent from T. az. Penicillin-binding and involved in peptidoglycan biosynthesis.
#ko:K00342 and ko:K05575 are ndHD subunits. 4/5 are core to T. az. 1/5 is in 0 T. az.
rm names_of_HOGs.tmp

## Now, get the gene_names and KOs that are each free_core_Taz_not_core and that do not have a corresponding annotation (KO for gene_name or gene_name for KO) that are Taz_core.

## Get gene_names that look like they are free_core_Taz_not_core have KOs that are also free_core_Taz_not_core
## Get kos corresponding to gene_names
> super_table_by_name.tsv
cat /home/liam/data/main_pipeline/core_vs_shell_analyses/gene_name/intact/free_core_Taz_not_core.list | \
while read name; do
cat $super_table | awk -v name=$name -F "\t" '$2==name' >> super_table_by_name.tsv
done
cat super_table_by_name.tsv | \
awk -F "\t" '$3!="-" {print $3}' | \
tr "," "\n" | \
sort -u \
> kos_of_names.tmp

## Subset to the ko's that also qualify
cat /home/liam/data/main_pipeline/core_vs_shell_analyses/KO/intact/free_core_Taz_not_core.list kos_of_names.tmp | \
sort | uniq -c | \
awk '$1==2 {print $2}' \
> names_and_kos_qualify.tmp
## Get the super_table entry for these and for gene_names that qualified and do not have a KO
cat super_table_by_name.tsv | awk '$3=="-"' > free_core_Taz_not_core_by_gene_name_and_KO.tsv

cat names_and_kos_qualify.tmp | \
while read ko; do
cat $super_table | awk -v ko=$ko -F "\t" '$3~ko' >> free_core_Taz_not_core_by_gene_name_and_KO.tsv
done


## Now do the reverse (lead with KOs)

## Get gene_names corresponding to KOs
cat /home/liam/data/main_pipeline/core_vs_shell_analyses/KO/intact/free_core_Taz_not_core.list | \
while read ko; do
cat $super_table | awk -v ko=$ko -F "\t" '$3==ko' >> super_table_by_ko.tsv
done
cat super_table_by_ko.tsv | \
awk -F "\t" '$2!="-" {print $2}' | \
tr "," "\n" | \
sort -u \
> names_of_kos.tmp

## Subset to the names that also qualify
cat /home/liam/data/main_pipeline/core_vs_shell_analyses/gene_name/intact/free_core_Taz_not_core.list names_of_kos.tmp | \
sort | uniq -c | \
awk '$1==2 {print $2}' \
> kos_and_names_qualify.tmp
## Get the super_table entry for these and for kos that qualified and do not have a KO
cat super_table_by_ko.tsv | awk '$2=="-"' >> free_core_Taz_not_core_by_gene_name_and_KO.tsv
cat kos_and_names_qualify.tmp | \
while read name; do
cat $super_table | awk -v name=$name -F "\t" '$2==name' >> free_core_Taz_not_core_by_gene_name_and_KO.tsv
done

## Some lines will be repeated because found when leading with gene_name and when leading with KO. So, get rid of these, and then organize by gene_name
cat free_core_Taz_not_core_by_gene_name_and_KO.tsv | \
sort -u | \
sort -k2,2r -k3,3 \
> sort.tmp
mv sort.tmp free_core_Taz_not_core_by_gene_name_and_KO.tsv

## Clean up
rm super_table_by_name.tsv kos_of_names.tmp names_and_kos_qualify.tmp super_table_by_ko.tsv names_of_kos.tmp kos_and_names_qualify.tmp


## Now, look at the genes where by name they are core to free-living, but not to T. az.
grep -Ff /home/liam/data/main_pipeline/core_vs_shell_analyses/gene_name/intact/free_core_Taz_not_core.list $super_table | less -S
## ccmK2: One HOG. core except FiliculoidesGal. But FiliculoidesGal has ccmK1.
## groL2: One HOG. core except Nilotica5001. But Nilotica5001 has ccmK1.
## exoD: 1 intact in 5 T. az.: only pseudo in Nilotica, Mexicana, Microphylla
## nifJ: fully absent in all T. az.!
## nadK1: intact in 7 T. az. pseudo in Nilotica5001, but Nilotica1 has nadK2 (see ko:K00858)
## moaE: intact in 7 T. az. absent in Nilotica5001
## melB: 0 intact, all 8 pseudo. FiliculoidesGal and Caroliniana3004 both fragmented and fragments not the same lenght. Did not look at rest. No other HOGs have same KO
## comA: 0 intact, all 8 pseudo. There is another HOG with the same KO, but 0 T. az. in that one.
## rpoZ: 1 intact (Nilotica5001), 7 pseudo. No other HOGs have same KO
## phoA: 4 intact (), 4 pseudo (Nilotica, Filiculoides, NosAzo, Rubra). Not assigned a KO
## gap1: 0 intact, 2 pseudo (Nilotica, NosAzo). 
