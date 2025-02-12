## Call assign_pseudogenes_to_hogs.sh for the main pipeline. 

## Make a N9 prefix list with both the long and short prefixes
## hashed out because doesn't need to re-run if this script reruns.
#grep -Ff ~/subtree.prefix.list ~/All.shortened.tsv > ~/subtree.shortened.tsv

## Inputs:
outdir=/home/liam/data/main_pipeline/assign_pseudogenes
N=N9
prefix_tsv=~/subtree.shortened.tsv
pseudofinder_outdir=/home/liam/data/main_pipeline/pseudofinder_out
orthofinder_downstream_dir=/home/liam/data/main_pipeline/orthofinder_out/for_downstream

## Call assign_pseudogenes_to_hogs.sh !
bash -i /home/liam/scripts/helper/assign_pseudogenes_to_hogs.sh \
-p $pseudofinder_outdir \
-f $orthofinder_downstream_dir \
-n N9 \
-l $prefix_tsv \
-o $outdir
