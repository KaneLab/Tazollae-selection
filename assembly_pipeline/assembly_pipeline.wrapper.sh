##Call assembly_pipeline.sh using variables from a .tsv

cat assemblyPipelineVariables.tsv | \
while read prefix sra assembler ratio; do
bash -i ~/scripts/assembly_pipeline.sh -s $sra -p $prefix -a $assembler -r $ratio
done

##Call mmseqs taxonomy on the bins
##Not in the above loop because faster to call once on all of the bins than
##for each prefix separately.
##I ran this code by copy paste, not by calling the wrapper or even having a script
##So could have some bugs.
bash -i assembly_pipeline_1a.sh
