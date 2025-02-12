## Takes a phylogenetic tree and a vector of continuous values corresponding to
## a continuous trait for the leaves of that tree

## Code for passing bash variables to R from:
## pauljohn32 answer edited Jun 26, 2019 on
## https://stackoverflow.com/questions/56777529/how-to-pass-bash-variable-into-r-script

##Parse command line arguments
##Took the following code for parsing command line arguments from:
##https://stackoverflow.com/questions/56777529/how-to-pass-bash-variable-into-r-script
cli = commandArgs(trailingOnly = TRUE) #trailingOnly means only consider args after --args
args = strsplit(cli, "=", fixed = TRUE)
for (e in args) {
  argname <- e[1]
  if (! is.na(e[2])) {
    argval <- e[2]
    ## regular expression to delete initial \" and trailing \"
    argval <- gsub("(^\\\"|\\\"$)", "", argval)
  }
  else {
    # If arg specified without value, assume it is bool type and TRUE
    argval <- TRUE
  }
  assign(argname, argval)
}

library(ape)
#library(phytools)
## Read the tree
my_tree = read.tree(tree_path)
## Read the data
my_data = read.table( table_path, row.names = 1, sep = "\t" )
## initiate output table to append to
outtable=c()

j=1
print(j)
## Subset the data
in_data = my_data[,j]
names(in_data)=rownames(my_data)
a=ace( x = in_data , phy = my_tree )

## Get values at each internal node
## The internal node values in $ace are in an order that seems to correspond with their order in my_tree$node.label,
## but they are not explicitly labeled as such.
## I checked this by looking at the results for a tree and noting that the results made sense with this order correspondence.
## Check was done with the discrete wrapper.
outtable=rbind(
	cbind( names(in_data), in_data[] ),
	cbind( my_tree$node.label,a$ace )
)


write.table( outtable, file=outfile, append=FALSE, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE )
