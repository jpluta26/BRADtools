# pluta 3/19/19
rm(list = ls())
setwd("C:/Users/jpluta/Desktop/TECAC-fulldata/Validation/ImputationValidation/")

# check the concordance between imputed snps in the replication data, and their genotyped counter parts
# in the validatoin data
library(ggplot2)
library(epiR)
# the replication data
rep.dat <- read.table("repGeno.raw", header = TRUE, as.is = TRUE)
rep.dat$IID <- substr(rep.dat$IID, 1, 12)

# the validation data
val.dat <- read.table("impVal-QC.raw", header = TRUE, as.is = TRUE)

rep.dat <- rep.dat[ rep.dat$IID %in% val.dat$IID,]

ind <- match(rep.dat$IID, val.dat$IID)
rep.dat <- rep.dat[ind,]

out <- data.frame( rho.c = rep(0, dim(val.dat)[2] - 6), lower = 0, upper = 0, n = 0)


# for each snp...
for( j in 7:dim(val.dat)[2] )
{
  # compare the minor allele- flip if they dont match  
  repAllele <- strsplit(colnames(rep.dat)[j], "_")[[1]][2]
  valAllele <- strsplit(colnames(val.dat)[j], "_")[[1]][2]
  
  if( repAllele != valAllele)
  {
    # in both cases, MAF, is right around 0.5, so i think these assignments are wrong
    #X8.95661772_C
    #X18.703785_G   
    # dont flip these- they are using the wrong minor allele
    if( !(j %in% c(19,39)))
    {
      ind0 <- which(val.dat[,j] == 0)
      ind2 <- which(val.dat[,j] == 2)
      val.dat[ind0,j] <- 2
      val.dat[ind2,j] <- 0
    }
    
  }
  
  # drop subjects with missing observations
  ind = !(is.na(val.dat[,j]) | is.na(rep.dat[,j]))
  
  # concordance correlation coefficient
  x = epi.ccc(rep.dat[,j], val.dat[,j], rep.measure = TRUE, subjectid = rep.dat$IID[ind])
  out[ j - 6,1:3] <- x$rho.c
  out$n[j - 6] <- sum(ind)
}

rownames(out) <- colnames(rep.dat)[7: dim(rep.dat)[2]]


# could add bland and altman plot if need be




# format and write data
snp.dat <- data.frame( rsid = colnames(rep.dat)[7:dim(rep.dat)[2]], overlap = out)

snp.dat$rsid <- gsub("^X", "", snp.dat$rsid)
ind = substr(snp.dat$rsid, 1, 1) == "."
snp.dat$rsid[ind] <- paste("23", snp.dat$rsid[ind], sep = "")

allele <- unlist(lapply(strsplit(snp.dat$rsid, "_"), function(x) x[2]))

snp.dat$rsid <- gsub("\\.", ":", snp.dat$rsid)
snp.dat$rsid <- sub("_C", "", snp.dat$rsid)
snp.dat$rsid <- sub("_G", "", snp.dat$rsid)
snp.dat$rsid <- sub("_T", "", snp.dat$rsid)
snp.dat$rsid <- sub("_A", "", snp.dat$rsid)

pre <- substr(snp.dat$rsid, 1, 2)

if( any(pre == 23) )
{
  snp.dat$rsid[ pre == 23 ] <- sub("23:", "X:", snp.dat$rsid[ pre == 23 ])
}

write.table(snp.dat, "impVal-res.txt", col.names = T, row.names = F, quote = F)
