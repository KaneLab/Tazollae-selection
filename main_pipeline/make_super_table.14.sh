## Combine several output tables into a single output table

cd /home/liam/data/main_pipeline

## Get the list of all HOGs
## And annotations

cat /home/liam/data/main_pipeline/annotations/eggnog-mapper_out/hog_annotations/hog_annotations.tsv > super_table.tsv

## Get hyphy RELAX results
echo -e "RELAX_result\tRELAX_q" > hyphy_cols.tmp
cat super_table.tsv | \
awk -F "\t" '{print $1}' | \
tail -n +2 | \
while read hog; do
cat /home/liam/data/main_pipeline/hyphy_out/relax_out/results.tsv | \
awk -v hog=$hog '
	BEGIN {f=0};
	$1==hog { print $5"\t"$4; f=1; exit }
	END { if (f==0) { print "-\t-" } }
' \
>> hyphy_cols.tmp
done

paste super_table.tsv hyphy_cols.tmp > paste.tmp
mv paste.tmp super_table.tsv
rm hyphy_cols.tmp

## Get the number of genomes intact, pseudogenized (and not also intact), absent for each HOG
## Loop free-living and Taz.
## "present" means that the improve_gene_inventory.11a.sh pipeline determined that the locus might be intact, but might be a pseudogene
indir=/home/liam/data/main_pipeline/gene_loss
orthofinder_dir=/home/liam/data/main_pipeline/orthofinder_out/for_downstream
echo "free_living
Taz" | \
while read group; do
echo -e "${group}_intact\t${group}_absent\t${group}_pseudo\t${group}_present" > ${group}_cols.tmp
cat super_table.tsv | \
awk -F "\t" '{print $1}' | \
tail -n +2 | \
while read hog; do
paste ${indir}/${group}.intact_vs_not.tsv ${indir}/${group}.presence_absence.tsv ${orthofinder_dir}/${group}.n_intact.tsv | \
awk -v hog=$hog '
BEGIN {intact=0; present=0; pseudo=0; absent=0};
$1==hog {
	for (i=2;i<=NF/3;++i) {
		if ($(i+(NF/3))==0) {absent+=1}
		else {
			if ($i==1) {
				if ($(i+(2*NF/3))!=0) {intact+=1}
				else {present+=1}
			}
			else {pseudo+=1}
		}
	};
	exit
};
END {print intact"\t"absent"\t"pseudo"\t"present}
' \
>> ${group}_cols.tmp
## Finish while read hog loop
done

paste super_table.tsv ${group}_cols.tmp > paste.tmp
mv paste.tmp super_table.tsv
rm ${group}_cols.tmp
## Finish while read group loop
done

## Check, there should be 0 "free_living_present"
cat super_table.tsv | awk -F "\t" '$(NF-4)!=0' #1 (header)
## Good! Don't need this column.


## Rearrange columns
cat super_table.tsv | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$14"\t"$15"\t"$16"\t"$17"\t"$10"\t"$11"\t"$12"\t"$8"\t"$9"\t"$4"\t"$5"\t"$6"\t"$7}' > awk.tmp
mv awk.tmp super_table.tsv
