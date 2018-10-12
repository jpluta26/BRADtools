#!/usr/bin/env Rscript

# pluta 10/12/18
#
# example wrapper script for getHRDScore functions
# call this script from the command line or shell script:
# Rscript runHRD.R sub.id seq.fname ploidu.fname


source("/Users/jpluta/Desktop/CODE/getHRDScore.R")


args = commandArgs(trailingOnly = TRUE)
if(length(args) < 3)
{
  print("need to provide 3 arguments:")
  print("1: sub.id = character string, the subject ID")
  print("2: seq.dat   = character string, the filename of the seq data")
  print("3: ploidy.dat  = character string, the filename of the ploidy data")
  stop()
}

sub.id <- args[1]
seq.fname <- args[2]
ploidy.fname <- args[3]


ploidy.dat <- read.table(ploidy.fname, header = TRUE)
seq.dat <- read.table(seq.fname, header = TRUE)

hrd <- getHRD.Data( sub.id, seq.dat, ploidy.dat )

write.table(hrd, paste(sub.id, "hrd.txt", sep = "_"), quote = FALSE,
            col.names = TRUE, row.names = FALSE, append = FALSE)