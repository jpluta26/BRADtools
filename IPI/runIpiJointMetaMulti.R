
# pluta 7/20/21
# v1.0


library(BEDMatrix) 
library(parallel)
library(pbapply)
library(matrixcalc)
library(Matrix)

source("~/IPI/jointMeta.R")
source("~/IPI/iraeAssoc.R")
source("~/gwas_tools/alignSnps.R")


# ================================================================================= #
# ======================================= functions =============================== #
# ================================================================================= #

# --------------------------------------------------------------------------------- #
prepData <- function(geno.dat, COVFILE, COVARS)
# function to set up genotypes with phenotype and covariate data
# input:
#   geno.dat (BEDMatrix object), genotype data from BED file
#   COVFILE (string), name of the file containing covariate data
#   COVARS (string), comma separated string of covariates, corresponding to columns in
#       COVFILE
#
# output:
#   list of data.frames, containing genotype and covariate data, truncated to the covariates
#   of interest
{
  # the list of variables we will always consider, plus the additional specified covariates
  if( COVARS == "NA")
  {
    covar.names <- c("PC1", "PC2", "PC3", "irae3", "Prior", "Surv_Months", "Vital_Status_2yrs")
  } else
  {
    covar.names <- c(strsplit(COVARS, ",")[[1]], "PC1", "PC2", "PC3", "irae3", "Prior", "Surv_Months", "Vital_Status_2yrs")
  }
  
  if( !file.exists(COVFILE))
  {
    stop(paste0(COVFILE, " not found. Exiting"))
  }
  
  # much more efficient  than trying to use the .raw file
  print("reading input...")
  
  cov.dat <- read.table(COVFILE, header = T, sep = ",")
  
  if(all(is.na(match(rownames(geno.dat), cov.dat$GWASID))))
  {
    stop("could not match geno.dat subject ids to cov.dat subjects ids")
  }
  
  if(!all(covar.names %in% colnames(cov.dat)))
  {
    print("the following covars were not found in the covariate data:")
    print(covar.names[!(covar.names %in% colnames(cov.dat))])
    print(paste0("colnames(cov.dat): ", colnames(cov.dat)))
    stop()
  }
  
  cov.dat <- cov.dat[match(rownames(geno.dat), cov.dat$GWASID),]
  cov.dat <- cov.dat[ ,colnames(cov.dat) %in% covar.names,]
  print("done")
  
  
  return(list(geno.dat,cov.dat))
  
}
# ---------------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------------- #
runJointMetaSNP <- function(snp, dat.list)
# run joint meta analysis for a given snp
# joint meta requires the full variance-covariance matrix between the two data sets;
# full model needs to be run
# 
# input: snp (string), snp of interest, should be a column of dat1 and dat2
# output: df (data.frame), summary statistics for the joint meta analysis of snp
{
  ind <- colnames(dat.list[[1]][[1]]) == snp
  
 
  k <- length(dat.list)
  
  fit <- list()
  for( i in 1:k )
  {
    fit[[i]] <- iraeAssoc.priorint( dat.list[[i]][[1]][,ind], dat.list[[i]][[2]], fit.only = TRUE)
  }

  # this assume all datasets have the same snps
  df <- jointMeta(fit, "geno.dat", "Prior", colnames(dat.list[[1]][[1]])[ind]) 
  
  # compute average maf
  n <- rep(0,k)
  maf <- rep(0,k)
  
  for( i in 1:k )
  {
    g <- dat.list[[i]][[1]][,ind]
    g <- g[!is.na(g)]
    n[i] <- length(g)
    maf[i] <- sum(g) / (length(g) * 2)
  }
  
  df$avg.maf <- ( t(n) %*% maf ) / ( sum(n) )
  
  return(df)
}
# ---------------------------------------------------------------------------------- #


# -------------------------------- substrRight ------------------------------------- #
substrRight <- function(x, n)
# helper function get  the last n character of x
{
  substr(x, nchar(x)-n+1, nchar(x))
}
# ---------------------------------------------------------------------------------- #

# ================================================================================== #
# ================================================================================== #
# ================================================================================== #



# ================================================================================== #
# ================================= MAIN  ========================================== #
# ================================================================================== #


args = commandArgs(trailingOnly=TRUE)

#  perform joint meta analysis of n datasets
if( length(args) < 2 )
{
  print("need to provide 2 arguments: PARAMFILE OUTNAME")
  print("PARAMFILE is a textfile of the format: ")
  print("BED file name")
  print("BIM file name")
  print("pheno file name")
  print("covariates")
  print("for each study")
  stop()
}


PARAMFILE <- args[1]
OUTNAME <- args[2]


setup.dat <- read.table(PARAMFILE, header = FALSE, stringsAsFactors = FALSE)
n.study <- dim(setup.dat)[1] / 4

bed.list <- list()
bim.list <- list()

for( i in 1:n.study )
{
  tmp <- setup.dat[(1 + (4 * (i - 1))):(4 + (4 * (i - 1))),]
  BEDFILE <- tmp[which(substrRight(tmp, 3) == "bed")]
  BIMFILE <- tmp[which(substrRight(tmp, 3) == "bim")]
  
  
  print(paste0("reading data set ", i))
  print(paste0("BEDFILE = ", BEDFILE))
  print(paste0("BIMFILE = ", BIMFILE))
  
  
  BIM <- read.table(BIMFILE, header = F)
  BED <- BEDMatrix(BEDFILE, simple_names = TRUE)
  
  # remove duplicated snps
  BED <- BED[,!duplicated(colnames(BED))]
  BIM <- BIM[!duplicated(BIM$V2),]
  
  if(dim(BED)[2] != dim(BIM)[1])
  {
    stop(paste0("dimension mistmatch between ", BED, " and ", BIM))
  }
  
  bim.list[[i]] <- BIM
  bed.list[[i]] <- BED
}

common.snps <- getCommonSnps( bim.list )
bed.list <- reduceBEDToCommonSet( bed.list, common.snps )
bim.list <- reduceBIMToCommonSet( bim.list, common.snps )


# bed.list[[1]] will be  used as the reference
ref  <- bim.list[[1]]

for( i in 2:length(bed.list))
{
  ind <- which(ref$V5 != bim.list[[i]]$V5 | ref$V6  != bim.list[[i]]$V6)
  
  # align snps in bim2 to bim1
  # dont need to align the bim files unless they are being written to output
  print(paste0("aligning bed", i))
  bed.list[[i]][ ,ind ] <- do.call(cbind, lapply(bim.list[[i]]$V2[ind], alignSnps, BIM1 = ref[ind,], BIM2 = bim.list[[i]][ind,], bed.list[[i]][,ind]))
  print("done")
}


# create the list of datasets
dat.list <- list()

# fetch the covariates for each study and attach to the genotype data
# sets up the data structure for runJointMetaSNP
for( i in 1:n.study )
{
  print(paste0("preparing study ", i, ": "))
  
  COVFILE <- setup.dat$V1[ 3 + 4 * (i - 1)]
  COVARS <- setup.dat$V1[ 4 + 4 * (i - 1)]
  
  print(paste0("COVFILE = ", COVFILE))
  print(paste0("COVARS = ", COVARS))
  dat.list[[i]] <- prepData(bed.list[[i]], COVFILE, COVARS)
  print("done")
}


print("setting up parallel processing..")
cl <- parallel::makeCluster(detectCores(), setup_strategy = "sequential")

# export necessary functions
clusterExport(cl, list("iraeAssoc.priorint", "jointMeta", "bdiag", "checkGenoFlip",
                       "flipGeno", "flipAllele", "speedglm", "alignSnps"))
print("done")


print("running meta analysis...")
out <- do.call(rbind.data.frame, pblapply(common.snps, runJointMetaSNP, dat.list, cl = cl))
print("done")

stopCluster(cl)

out$chr <- unlist(lapply(strsplit(as.character(out$snp), ":"), function(x) x[1]))

print(paste0("writing to file: ", OUTNAME))
write.table(out, OUTNAME, col.names = TRUE, row.names = FALSE, append = FALSE, quote = FALSE)
print("done")
# ================================================================================== #
