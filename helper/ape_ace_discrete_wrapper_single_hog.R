## Takes a phylogenetic tree and a vector of binary values corresponding to
## a discrete trait for the leaves of that tree
## Returns the ML no-gene-gain probability that the each node in the tree
## including leaves (which is just the input) are positive for that trait.
## To be used in determining where genes are lost.

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
my_data = read.table( table_path, header = FALSE , row.names = 1, sep = "\t" )
in_data = my_data[,1]
names(in_data) = row.names(my_data)

## No gene gain allowed:
a=ace(x = in_data , phy = my_tree , type = "discrete" , model = matrix(c(0, 1, 0, 0), 2) , marginal=FALSE , use.expm = TRUE, use.eigen=FALSE )
## phytools rerootingMethod should be the same as ace without the use.flags flipped, I think.
	
## Get the probabilities for a positive trait at each node.
## The internal node values in $lik.anc are in an order that seems to correspond with their order in my_tree$node.label,
## but they are not explicitly labeled as such.
## I checked this by looking at the results for a tree and noting that the results made sense with this order correspondence.
new_row=unname(
	t(
	cbind(
		rbind( names(in_data), in_data[] ),
		rbind( my_tree$node.label,a$lik.anc[,colnames(a$lik.anc)=="1"] )
	)
	)
)


write.table( new_row, file=outfile, append=FALSE, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE )
