#!/usr/bin/env Rscript
# pluta 9/18/18

library(data.table)
library(qqman)

# script to create manhattan plot of either a single chromosome of the whole genome
# this script greatly increases the speed of plotting by subsampling a small fraction
# of the non-significant snps

# INPUT: INFILE (text file) with two columns, MarkerName (chr:pos) or rsid and 
# P-value
#
# OUTPUT: OUTFILE (string), the name of the pdf to output

args = commandArgs(trailingOnly=TRUE)
if( length(args) < 2 )
{
		print("need to provide 2 arguments: INFILE OUTFILE")
		print("INFILE should be a ggman format text file, with two columns: MarkerName (chrpos or rsid) and P-value")
        stop()
}

INFILE  = args[1]
OUTFILE = args[2]


# PARAMS
# the threshold for p-values; sub-sample any values with p > p.thresh
p.thresh = 0.01

# proportion of snps to keep
k.thresh = 0.2


MANOUTFILE <- paste(OUTFILE, "_man.png", sep="")

# ---

print("reading association testing data...")

dat <- fread(INFILE, header=TRUE)
dat <- as.data.frame(dat)

# check input validity
# important for consistency with locuszoom
if( colnames(dat)[1] != "MarkerName")
{
	stop("Colnames of INFILE should be 'MarkerName' and 'P-value' ")
}


dat$MarkerName <- as.character(dat$MarkerName)
dat$CHR <- unlist(lapply(strsplit(dat$MarkerName, ":"), function(x) x[1]))
dat$BP <- unlist(lapply(strsplit(dat$MarkerName, ":"), function(x) x[2]))
dat$BP <- as.numeric(dat$BP)

if( any(dat$CHR == "X"))
{
        dat$CHR[which(dat$CHR == "X")] <- 23
}

dat$CHR <- as.numeric(dat$CHR)
colnames(dat) <- c("rsid", "p", "chr", "bp")
chrlabs <- as.character(seq(1:length(unique(dat$CHR))))
dat$p <- as.numeric(dat$p)

# remove any rsids that snuck in
dat <- dat[!is.na(dat$chr),]
			
print("done!")




# ---

print("subsampling association data...")

# s is the sub sample
x <- dim(dat[which(dat$p > p.thresh),])[1]
s <- sample(dat$rsid[which(dat$p > p.thresh)], round(k.thresh * x), replace=F)

print("done!")


# ---

print("creating manhattan plot...")

png(MANOUTFILE, height=600, width=1200)
manhattan(dat[(dat$rsid %in% dat$rsid[match(s, dat$rsid)]) | (dat$p <= p.thresh),],
        chr="chr", bp="bp", p="p", snp="rsid", col=c("gray", "black"), ylim=c(0,12),
        logp=TRUE, chrlabs = chrlabs)
dev.off()

print("done!")



