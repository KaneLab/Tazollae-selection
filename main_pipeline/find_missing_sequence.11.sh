bwa index /home/liam/data/fastas/reference.fna

cat MAGList | while read sample; do

	bwa mem -t 24 /home/liam/data/fastas/reference.fna /home/liam/data/fastas/Dijkhuizen_genomes/$sample | samtools view -bhu - | samtools sort -@ 24 - -o ${sample}.sorted.bam

	#Index the completed alignment map
	samtools index ${sample}.sorted.bam

	#samtools depth omits certain sites by default, so we can get all sites with the -aa (absolutely all) flag
	samtools depth -aa ${sample}.sorted.bam > ${sample}.depth

done


echo "Locus	Pos" > Header

counter=0
cat MAGList | while read sample; do

	if [[ $counter -eq 0 ]]; then
		#If first loop, then get the chrom and position columns
		awk '{print $1"	"$2"	"$3}' ${sample}.depth > TmpMaster

	else
		#If not the first loop, then only get the depth values and add them to the Master depths file
		awk '{print $3}' ${sample}.depth > DepthCol
		paste TmpMaster DepthCol > tmp; mv tmp TmpMaster

	fi

	awk -v name=$sample '{print $0"	"name}' Header > TmpHead; mv TmpHead Header
	counter=$(( $counter + 1 ))

done

cat Header TmpMaster > MAG_Alignments.depth
#rm Header; rm TmpMaster; rm DepthCol
