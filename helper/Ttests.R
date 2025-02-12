## Run through bash
## Run a on a 2-col tsv input
## and append to an output file.
## Specify with header= and paired= in command line
## if header on input and if paired ttest, respectively
## These input data are all saved *ttest.tsv in 
## liam@photobiont.colorado.edu:/home/liam/data/main_pipeline/RefSeq_vs_prokka/analyses
## and its subdirectories.

library(broom)
setwd("/Users/liamfriar/Desktop/Taz_stuff/analyses")

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

paired = as.logical(paired)
header = as.logical(header)


intable=read.table( infile , header=header, sep="\t")
ingroup=intable[!is.na(intable[,1]),1]
outgroup=intable[!is.na(intable[,2]),2]
a=tidy( t.test( ingroup, outgroup, paired=paired ) )
write (infile, file=outfile , append=TRUE)
write( a$statistic , file=outfile , append=TRUE)
write( a$p.value , file=outfile , append=TRUE)
     