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
my_data = read.table( table_path, header = TRUE , row.names = 1, sep = "\t" )
## initiate output table to append to
outtable=c()

## Loop through columns of the input table
#for ( j in 1:10 ) {
for ( j in 1:ncol(my_data) ) {
	print(j)
	## Subset the data
	hog = colnames(my_data[j])
	in_data = my_data[,j]
	names(in_data)=rownames(my_data)
	if (approach=="no_gain") {
	## No gene gain allowed:
	a=ace(x = in_data , phy = my_tree , type = "discrete" , model = matrix(c(0, 1, 0, 0), 2) , marginal=FALSE , use.expm = TRUE, use.eigen=FALSE )
	} else if (approach=="ARD") {
	## Gene gain allowed. "ARD" makes gain/loss rates different, wile "ER" would make them equal.
	#a=ace(x = in_data , phy = my_tree , type = "discrete" , model = "ARD" , marginal=FALSE )
	## use.expm = TRUE, use.eigen=FALSE is recommended by the ace documentation for non-symmetric transition matrices, I think.
	## It is not 100% clear to me, though.
	a=ace( x = in_data , phy = my_tree , type = "discrete" , model = "ARD" , marginal=FALSE , use.expm = TRUE, use.eigen=FALSE )
	}
	## phytools rerootingMethod should be the same as ace without the use.flags flipped, I think.
	
	## Get the probabilities for a positive trait at each node.
	## The internal node values in $lik.anc are in an order that seems to correspond with their order in my_tree$node.label,
	## but they are not explicitly labeled as such.
	## I checked this by looking at the results for a tree and noting that the results made sense with this order correspondence.
	new_row=unname(
		t(
		rbind(
			cbind("HOG", hog),
			cbind( names(in_data), in_data[] ),
			cbind( my_tree$node.label,a$lik.anc[,colnames(a$lik.anc)=="1"] )
		)
		)
	)

outtable=rbind(outtable, new_row)
## End loop through columns of input
}


write.table( outtable, file=outfile, append=FALSE, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE )
