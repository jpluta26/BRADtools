#!/usr/bin/env Rscript
# pluta 10/14/19

# step 2 of custom annotation with garfield. before running this, put summary stats in
# ~/garfield-data/pval/TECAC
# and put annotations in ~/garfield-data/bed

# the GARFIELD custom annotation script included in v2 doesnt work! 
# use this instead

args = commandArgs(trailingOnly = TRUE)
if( length(args) == 0 )
{
        stop()
}

# bed file containing annotation data
INPUTBEDFILE=args[1]

library(data.table)
library(IRanges)

bed <- fread(INPUTBEDFILE, skip = 1)

for( chr in c(1:22, "X"))
{
      # the variants that garfield has maf and TSSd for (from 1kg)
        tmp <- fread(paste0("~/garfield-data/maftssd/chr", chr))
        allVar <- cbind(tmp$V1, 0)
        dat <- bed[ bed$V1 == paste0("chr", chr),]

        # sometimes end points will come before end points- clear error
        # remove these
        if( any(dat$V2 > dat$V3) )
        {
                ind = dat$V2 < dat$V3
                dat <- dat[ind,]
        }

        # find snps overlapping annotations
        rangeA <- IRanges( allVar[,1], allVar[,1] )
        rangeB <- IRanges(dat$V2, dat$V3)

        overlap <- findOverlaps(rangeA, rangeB, type = "within")
        allVar[overlap@from,2] <- 1


        tmp = strsplit(INPUTBEDFILE, "/")
        BED = tmp[[1]][length(unlist(tmp))]
        write.table(allVar, paste0("~/garfield-data/annotation-custom/", BED, chr),
                col.names = FALSE,
                row.names = FALSE, sep = "\t", quote = FALSE)
}
