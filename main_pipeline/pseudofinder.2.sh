## Run pseudofinder on prokka outputs, not including the outgroups which will are only needed to root the species tree for orthofinder.
## Run pseudofinder first without the --use_deviation flag, then with
## Running twice to make sure that the outputs look as expected.
## Because personal correspondence with the authors of pseudofinder revealed that this is a newer tool

indir=/home/liam/data/main_pipeline/prokka_out
outdir=/home/liam/data/main_pipeline/pseudofinder_out
mkdir $outdir
conda deactivate
conda activate pseudofinder


############
## RefSeq ##
############

##Run pseudofinder on 48 Nostocales genomes from NCBI, including T. az. 0708
##For T. az., use the other 47 genomes as the db.
##For all others, use the 46 genomes other than T. az. and themselves
##Run pseudofinder iteratively on the 48 nostocales genomes without the use_deviation flag
##Create database of the 47 Nostocales ncbi genomes other than T. az.
#Not de-replicating
#(would just unwrap to a single line, sort -u and split back into 2 lines)
#because I think I want to hit multiple times if the same gene is in multiple genomes.
#Add to the WP-id so all still unique in case that would mess up pseudofinder.
#Loop through the 47 free-living Nostocales
grep -v "^NostocAzollae" ~/Nostocales_incl0708.shortened.tsv  | \
while read query short_query; do
#Make output folder
mkdir ${outdir}/$query
cd ${outdir}/$query

#Make db
db=${outdir}/$query/allNostocExceptTazAnd.${query}.faa
> $db
##Loop through the free-living nostocales, adding the protein.faa to the database
##Except for the query protein.faa itself.
grep -v "^NostocAzollae" ~/Nostocales_incl0708.shortened.tsv  | \
while read prefix short_prefix; do
if [[ $prefix != $query ]]; then
cat ${indir}/${prefix}/${short_prefix}.faa | \
awk -v prefix=$prefix -F " " ' {if ($0 ~ "^>") {printf $1"_"prefix" "; for (i=2; i<NF; ++i) {printf $i" "}; print $NF} else { print $0} }' \
>> $db;
fi
done #prefix loop for db build


# Could call prep diamond database (if want to use skpdb flag in pseudofinder)
# But this is unnecessary because I am making a new database each time
# db=allNostocExceptTazAnd.${query}.faa
# diamond makedb --in $db -d ${db}.dmnd
# call pseudofinder
# db=${db}.dmnd
in_gbk=${indir}/${query}/${short_query}.gbk
## Don't use use_deviation flag
python ~/tools/pseudofinder/pseudofinder.py annotate -t 20 -ce --diamond -hc 46 --genome $in_gbk --outprefix ${query}_woDeviation --database $db;
## Use use_deviation flag
python ~/tools/pseudofinder/pseudofinder.py annotate -t 20 -ce --use_deviation -l 0.65 --diamond -skpdb -hc 46 --genome $in_gbk --outprefix ${query} --database $db;

## End query loop for pseudofinder run wrapper
done


############
## T. az. ##
############

## Generate database for all of the MAGs for pseudofinder
## Which is all of the free-living Nostocales (so 47 instead of 46 as above)
cd ${outdir}
#Make db
db=${outdir}/allFreeLivingNostoc.faa
> $db
##Loop through the free-living nostocales, adding the protein.faa to the database
##Except for the query protein.faa itself.
grep -v "^NostocAzollae" ~/Nostocales_incl0708.shortened.tsv  | \
while read prefix short_prefix; do
cat ${indir}/${prefix}/${short_prefix}.faa | \
awk -v prefix=$prefix -F " " ' {if ($0 ~ "^>") {printf $1"_"prefix" "; for (i=2; i<NF; ++i) {printf $i" "}; print $NF} else { print $0} }' \
>> $db;
done #prefix loop for db build
#Make a diamond database because this same database will be used by pseudofinder for all the T. az. genomes
diamond makedb --in $db -d ${db}.dmnd
db=${db}.dmnd

## Run the reference genome (separate because of the short prefix for prokka outputs)
query=NostocAzollae0708
mkdir ${outdir}/${query}
cd ${outdir}/${query}
in_gbk=${indir}/${query}/NosAzo0708.gbk
## Don't use use_deviation flag
python ~/tools/pseudofinder/pseudofinder.py annotate -t 20 -ce --diamond -hc 46 --genome $in_gbk --outprefix ${query}_woDeviation --database $db;
## Use use_deviation flag
python ~/tools/pseudofinder/pseudofinder.py annotate -t 20 -ce --use_deviation -l 0.65 --diamond -skpdb -hc 46 --genome $in_gbk --outprefix ${query} --database $db

## Loop through the MAGs
cat ~/Dijkhuizen.prefix.list | \
while read query; do
mkdir ${outdir}/${query}
cd ${outdir}/${query} 
in_gbk=${indir}/${query}/${query}.gbk
## -ce flag adds contig ends as potential intergenic regions.
## Using this flag to get the most possible hits.
## Don't use use_deviation flag
python ~/tools/pseudofinder/pseudofinder.py annotate -t 20 -ce --diamond -hc 46 --genome $in_gbk --outprefix ${query}_woDeviation --database $db;
## Use use_deviation flag
python ~/tools/pseudofinder/pseudofinder.py annotate -t 20 -ce --use_deviation -l 0.65 --diamond -skpdb -hc 46 --genome $in_gbk --outprefix ${query} --database $db;
## End query loop for pseudofinder run wrapper
done


##################################
## Check the use_deviation flag ##
##################################

##Check that the length-based pseudos with deviation are a subset of those without deviation
cd $outdir
> length_pseudogenes_wVSwoDeviation.txt
cat ~/All.prefix.list | \
while read prefix; do
grep "candidate" ${prefix}/${prefix}_woDeviation_pseudos.gff | grep "average length" | awk '{ for (i=1; i<7; ++i) printf $i"\t"; print $7 }' > lengthPseudos_woDeviation.tmp
grep "candidate" ${prefix}/${prefix}_pseudos.gff | grep "average length" | awk '{ for (i=1; i<7; ++i) printf $i"\t"; print $7 }' > lengthPseudos_wDeviation.tmp
overlap=$(cat lengthPseudos_woDeviation.tmp lengthPseudos_wDeviation.tmp | sort | uniq -c | sort -k1,1nr | awk 'BEGIN {c=0}; $1==2 {c+=1}; END {print c}')
initial_pseudos=$(grep "Initial pseudogenes:" ${prefix}/${prefix}_woDeviation_log.txt | sed 's/Initial pseudogenes:[[:space:]]*//')
echo $prefix >> length_pseudogenes_wVSwoDeviation.txt
echo "Initial pseudos: ${initial_pseudos}" >> length_pseudogenes_wVSwoDeviation.txt
echo "length pseudos w/o Deviation: $(wc -l < lengthPseudos_woDeviation.tmp)" >> length_pseudogenes_wVSwoDeviation.txt
echo "length pseudos w/ Deviation: $(wc -l < lengthPseudos_wDeviation.tmp)" >> length_pseudogenes_wVSwoDeviation.txt
echo "overlap: ${overlap}/$(wc -l < lengthPseudos_wDeviation.tmp)" >> length_pseudogenes_wVSwoDeviation.txt
done
rm *.tmp

#The deviation flag appears to have worked,
#so rename the various folders to reflect that wDeviation are to be used downstream
cat ~/All.prefix.list | \
while read prefix; do
rm ${prefix}/${prefix}_woDeviation*
done


###############
## Filenames ##
###############

##Change the headers in pseudofinder_out/${prefix}/${prefix}_pseudos.fasta
##To be shorter and easier to read
##Save in new file pseudofinder_out/${prefix}/${prefix}_pseudos.shorter_headers.fasta
cd /home/liam/data/main_pipeline/pseudofinder_out
cat ~/All.shortened.tsv | \
while read long_prefix short_prefix; do
cat ${long_prefix}/${long_prefix}_pseudos.fasta | sed 's/_[^_]\+_pseudo/_p/' > ${long_prefix}/${long_prefix}_pseudos.shorter_headers.fasta
done
