##Call assembly_pipeline.sh using variables from a .tsv

cat assemblyPipelineVariables.tsv | \
while read prefix sra assembler ratio; do
bash -i ~/scripts/assembly_pipeline.alreadyHaveReads.sh -s $sra -p $prefix -a $assembler -r $ratio
done
