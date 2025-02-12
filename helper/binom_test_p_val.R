## Calls R binom.text() with input x, n, p corresponding to what they are in the binom.test() definition
## appends p.value to outfile

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

n = as.integer(n)
x = as.integer(x)
p = as.numeric(p)


result=binom.test(x=x, n=n, p=p, alternative="two.sided")
write( result$p.value , file=outfile , append=TRUE)
