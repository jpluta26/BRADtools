rm(list = ls())

setwd("C:/Users/jpluta/Desktop/Melanoma/IpiNivo/")
library(ggplot2)
library(utils)

args = commandArgs(trailingOnly = TRUE)

if( length(args) < 1)
  {
    print("need to provide argument: MAF.THRESH")
    print("if the difference in MAF between batch1 and batch2 is > MAF.THRESH, tag snp for removal")
    stop()
  }

MAF.THRESH <- args[1]

# --------------------------- flipAllele ------------------------------------------- #
flipAllele <- function( allele )
  # function to orient strand to match reference (snplst)
{
  if( allele == "A")
  {
    return("T")
  } else
    if( allele == "G")
    {
      return("C")
    } else
      if( allele == "C" )
      {
        return("G")
      } else
        if( allele == "T")
        {
          return("A")
        } else
          return(allele)
}
# ---------------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------------- #
alignSNPs <- function( dat1, dat2 )
{
  # dat1 and dat2 are plink.frq files
  # only dat2 will be edited!!!
  ind = which(dat1$A1 != dat2$A1)
  pb <- txtProgressBar(min = 0, max = length(ind), style=3)
  
  for(i in ind)
  {
    # create a set of the alleles for snp j from both datasets
    set1 <- c(as.character(dat1$A1[i]), as.character(dat1$A2[i]))
    set2 <- c(as.character(dat2$A1[i]), as.character(dat2$A2[i]))
    
    # if both datasets have the same set of alleles, it means maj/min is flipped
    # flip the alleles and adjust MAF
    if( all(set1 %in% set2) )
    {
      temp <- dat2$A2[i]
      dat2$A2[i] <- dat2$A1[i]
      dat2$A1[i] <- temp
      dat2$MAF[i] <- 1 - dat2$MAF[i]
    } else
      # if they dont match, it means strand needs to be flipped
    {
      dat2$A1[i] <- flipAllele(dat2$A1[i])
      dat2$A2[i] <- flipAllele(dat2$A2[i])
      # then check alignment again and adjust
      if( dat2$A1[i] != dat1$A1[i] )
      {
        temp <- dat2$A2[i]
        dat2$A2[i] <- dat2$A1[i]
        dat2$A1[i] <- temp
        dat2$MAF[i] <- 1 - dat2$MAF[i]
      }
    }
    setTxtProgressBar(pb, i)
  }
  
  return(dat2)  
}
# ---------------------------------------------------------------------------------- #



# ------------------------------- MAIN --------------------------------------------------- #

# PLINK FREQ FILE:
# CHR SNP A1(min) A2(Maj) MAF NCHROBS

batch1 <- read.table("batch1.frq",
                     colClasses=c("character", "character", "character", "character", "numeric", "integer"), 
                    header=TRUE)
batch2 <- read.table("batch2.frq", 
                     colClasses=c("character", "character", "character", "character", "numeric", "integer"), 
                     header=TRUE)

# make sure unobserved snps are removed
batch1 <- batch1[which(batch1$MAF > 0) ,]
batch2 <- batch2[which(batch2$MAF > 0) ,]

# make sure indels and deletions havent been removed
batch1 <- batch1[ !(batch1$A1 %in% c("I", "D") | batch1$A2 %in% c("I", "D")), ]
batch2 <- batch2[ !(batch2$A1 %in% c("I", "D") | batch2$A2 %in% c("I", "D")), ]

# reduce to common set
keep <- intersect(batch1$SNP, batch2$SNP)

batch1 <- batch1[ batch1$SNP %in% keep,]
batch2 <- batch2[ batch2$SNP %in% keep,]



#
if( dim(batch1)[1] == 0)
{
  print("batch1 has 0 rows")
  stop()
}

if( dim(batch2)[1] == 0)
{
  print("batch2 has 0 rows")
  stop()
}


# make sure the reference allele is the same for both sets
batch2 = alignSNPs( batch1, batch2 )

# flag any snps with a MAF difference greater than 10% for removal
# this value is arbitrary- probably sufficient to capture random error
flagged.snps = batch2[which(abs(batch1$MAF - batch2$MAF) > MAF.THRESH),]

# plot
ind = which(abs(batch1$MAF - batch2$MAF) > MAF.THRESH)
dat <- data.frame(batch1=batch1$MAF, batch2=batch2$MAF)
p1 <- ggplot(data=dat,aes(x=batch1, y=batch2)) + geom_point(alpha = 0.25, color="blue") + 
      geom_point(data=dat[ind,], aes(x=batch1, y=batch2), color="red") +
  labs(x = "batch1 - MAF by SNP",
       y = "batch2 - MAF by SNP")  + xlim(0, 1) + ylim(0, 1) 

png("batch1.v.batch2-MAF-0.1_diff.png")
 print(p1)
dev.off()

# write snps to remove 
write.table(flagged.snps$SNP, "mismatchedsnps.txt", append=FALSE, quote=FALSE, col.names = FALSE, row.names = FALSE)
