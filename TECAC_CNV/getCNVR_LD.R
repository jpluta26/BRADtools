

# pluta 8/31/2021

# script to compute LD between a CNVR and a given snp from the corresponding genotype data
# LD reported as pearson's correlation coefficient squared
library(data.table)
# hit on ATF7IP  is 12:14633223; CNVRID 1975

args = commandArgs(trailingOnly =  TRUE)
if( length(args) < 3 )
{
  print("USAGE:: ")
  print("CNVRMAT: The matrix of CNVR data, generated in runCNVRassoc.R (eg dup_cnvr.mat")
  print("CNVRID: The CNVR of interest in the form of CNVR_###")
  print("GENOFILE: file of genotypes from the genotype data- this is extracted as part of getCNRV_LD.sh")
  stop()
}

# have to combine dup and del

CNVRMAT = args[1]
CNVRID = args[2]
GENOFILE = args[3]

M  <- as.data.frame(fread(CNVRMAT, header = TRUE))
M$IID <- paste0(M$IID, "_", M$IID)
geno <- as.data.frame(fread(GENOFILE, header = TRUE ))

geno <- geno[ match(M$IID, geno$FID), ]
geno <- geno[ ,-grep("HET", colnames(geno))]


if(length(which(colnames(M) == CNVRID)) == 0)
{
  stop(paste0("CNVRID ", CNVRID, " not found."))
}

g <- geno[,2]
m <- M[ ,which(colnames(M) == CNVRID)  ] - 2

ind <- which(m != 0)
r <- cor(g[ind], m[ind], method = "pearson")

print(paste0("r2 = ", r^2))



