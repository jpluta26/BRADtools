# pluta 12/3/19
# create a heat map comparing callers


#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(ggplot2)

args = commandArgs(trailingOnly = TRUE)
if( length(args) < 1)
{
  print("need to provide 2 arguments: INFILENAME OUTFILENAME")
  stop("INFILENAME = the caller data, eg somatic.high_confidence.hg19_multianno.report.txt")
  #stop("OUTFILENAME = the ROOT of the output image and text file, eg TEST for TEST.png and TEST-stats.txt")
}

INFILENAME <- args[1]
#OUTFILENAME <- args[2]

caller_stats <- function(INFILENAME){
  OUTFILENAME<-sub("txt","out",INFILENAME)
  dat <- read.table(INFILENAME, header = TRUE, sep = "\t", as.is = TRUE)
  #dat <- read.table("somatic.high_confidence.hg19_multianno.report.tsv", header = TRUE, sep = "\t", as.is = TRUE)

  caller.list <- c("MUTECT2", "STRELKA2", "VARDICTJAVA", "VARSCAN2")

  dat <- dat[ ,colnames(dat) %in% caller.list]

  # recode anything thats not PASS or . to FAIL
  for( i in caller.list )
  {
    dat[[i]][ !(dat[[i]] %in% c("PASS", ".")) ] <- "REJECT"
  }

  dat <- as.data.frame(apply(dat, 2, as.factor))


  # statistically compare the call rates of each pair of callers
  # use a fisher test- might be a smarter way to do this
  comp.mat <- matrix(0, nrow = 4, ncol = 4)

  for( i in 1:length(colnames(dat)) )
  {
    for( j in 1:length(colnames(dat)) )
    {
      t1 <- as.matrix(table(dat[[ colnames(dat)[i] ]]))
      t2 <- as.matrix(table(dat[[ colnames(dat)[j] ]]))
      comp.mat[i,j] <- chisq.test( cbind(t1, t2) )$p.value
    }
  }

  out <- c()
  for( i in 1:4 )
    {
      out <- cbind(out, table(dat[,i]))
    }


  rownames(comp.mat) <- caller.list 
  colnames(comp.mat) <- caller.list

  # format data for plotting
  # -convert to long format
  dat <- melt(dat, measure.vars = colnames(dat))
  n <- dim(dat)[1] / 4
  dat$x <- rep( seq(1:n), 4)
  dat$value <- as.factor(dat$value)
  levels(dat$value) <- c("NO CALL", "PASS", "REJECT")

  write.csv(table(dat$variable, dat$value),paste0(OUTFILENAME, "-count.csv"), row.names =T, col.names = T, append = F, quote = F)

  p1 <- ggplot(data = dat, aes(x = x, y = variable, fill = value)) + 
    geom_tile() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()
           ) +
    guides(fill = guide_legend(title= "Call Status")) +
    ylab("Caller") +
    scale_fill_manual(labels = c("NO CALL", "PASS", "REJECT"), values = c("blue", "green", "red"))

  png(paste0(OUTFILENAME, ".png"))
  print(p1)
  dev.off()

  write.table(format(comp.mat, scientific = TRUE, digits = 5), paste0(OUTFILENAME, "-stats.txt"), row.names =T, col.names = T, append = F, quote = F)
}

caller_stats(INFILENAME) 
