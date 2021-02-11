rm(list = ls())

# pluta 9/14/20
# general scriptto run all 4 kinds of ipi models (noprior - irae, priorint - irae, 
# noprior - coxph, priorint - coxph)
# writes output  with everything required for meta/pathway analysis

# hpc is particular about where packages are installed- doesnt always carry over to nodes
if(!require(pbapply)) 
{
  install.packages("pbapply", repos='http://cran.us.r-project.org')
}

if(!require(BEDMatrix))
{
  install.packages("BEDMatrix", repos='http://cran.us.r-project.org')
}

library(BEDMatrix) 
library(parallel)
library(pbapply)

args = commandArgs(trailingOnly = TRUE)

if( length(args) < 7)
{
  print("need to provide 7 arguments: BEDFILE COVARS COVFILE CHR MODEL MAF OUTNAME")
  print("BEDFILE: plink .bed file of genotype data, must have corresponding .bim and .fam files")
  print("COVARS: list of covariates to use in the model, in a comma seperated list")
  print("COVFILE: file of covariate data (ipi.nivo.pheno.txt)")
  print("CHR: numeric chromosome of interest")
  print("Model: choice of 1-5:")
  print("1: irae, no prior treatment")
  print("2: irae, prior interaction")
  print("3: coxph, no prior treatemt")
  print("4: coxph, prior interaction")
  print("5: irae, all subjects")
  print("OUTPREFIX: prefix attached to all output files")
  print("")
  print("If no covariates are used, set COVARS to NA")
  print("")
  print("")
  print(paste0("you provided ", length(args), " arguments:"))
  print(paste0("BEDFILE = ", args[1]))
  print(paste0("COVARS = ", args[2]))
  print(paste0("COVFILE = ", args[3]))
  print(paste0("CHR = ", args[4]))
  print(paste0("MDL = ", args[5]))
  print(paste0("MAF = ", args[6]))
  print(paste0("OUTPREFIX = ", args[7]))
  stop("missing necessary arguments")
}

BEDFILE = args[1]
COVARS = args[2]
COVFILE = args[3]
CHR = args[4]
MDL=args[5]
MAF=args[6]
OUTPREFIX=args[7]



# =============================================================================== #
# ================================= functions =================================== #
# =============================================================================== #


# -------------------------------- iraeAssoc ------------------------------------ #
# create a framework for generic processing code to be used across all studies; eg needs
# to accomdate different selected variables
iraeAssoc <- function( geno.dat, cov.dat )
  # can be used for prior, no prior, or all subjects - subset geno.dat when passing to function
  # input:
  #   geno.dat (data.frame) - snp matrix of genotype data
  #   cov.dat (data.frame) -  covariate data that sepcifies phenotype and all other covariates
  #         the actual covariates are selected a priori and are hard-coded into their
  #         respective functions
  #
  # output:
  #   a table of two columns, snp name and corresponding p-value of the genotype term from
  #   the logistic regression model
{
  
  # if MAF filter wasnt applied, this will throw an error
  if( sum(!is.na(unique(geno.dat))) == 1 )
  {
    return(NA)
  }
  
  if( all(is.na(geno.dat)) )
  {
    return(NA)
  }
  
  if( any(is.na(geno.dat)))
  {
    ind <- which(is.na(geno.dat))
    geno.dat <- geno.dat[-ind]
    cov.dat <- cov.dat[-ind,]
  }

  
  # convert categorical variables from characters to factors, so that they may be included
  # in the design matrix
  ind <- which(!(colnames(cov.dat) %in% c("PC1", "PC2", "PC3", "Prior", "irae3")))
  for( j in ind)
  {
    if( is.character(class(cov.dat[,j])))
    {
     cov.dat[,j] <- as.integer(as.factor(cov.dat[,j]))
    }
  }
  
  # pack all covariates into a single design matrix; exclude any covariates not wanted for
  # the model
  x <- as.matrix(cov.dat[,!(colnames(cov.dat) %in% c("Prior", "irae3", "Surv_Months", "Vital_Status_2yrs"))])
  y <- cov.dat$irae3
  
  print(paste0("run model with covariates; ", colnames(x)))
  # have to pass the whole df and specify the covariates of interest
  fit  <- glm(data = cov.dat, y ~ geno.dat + x , family = "binomial")
  
  b <- summary(fit)$coefficients[2,1]
  b.se <- summary(fit)$coefficients[2,2]
  p <- summary(fit)$coefficients[2,4]
  
  df <-data.frame(b = b, b.se = b.se, p = p)
  colnames(df) <- c("b", "b.se", "b.p")
  return(df)
}
# ------------------------------------------------------------------------------- #


# --------------------------- iraeAssoc.priorint ------------------------- #
iraeAssoc.priorint <- function( geno.dat, cov.dat)
  #   p will be significant if interaction or main effect of genotype are significant
{
  # when using multiple cores, the environment doesnt copy to each core
  # need to call libraries within function
  library(lmtest)
  if( sum(!is.na(unique(geno.dat))) == 1 )
    # all 1 genotype, cant test
  {
    return(NA)
  }
 
  if( all(is.na(geno.dat)) )
  {
    return(NA)
  }
  
  if( any(is.na(geno.dat)))
  {
    ind <- which(is.na(geno.dat))
    geno.dat <- geno.dat[-ind]
    cov.dat <- cov.dat[-ind,]
  }
  

  print(paste0("length(geno.dat) = ", length(geno.dat)))
  print(paste("dim(cov.dat[-ind]): ", dim(cov.dat)[1]))
  ind <- which(!(colnames(cov.dat) %in% c("PC1", "PC2", "PC3", "Prior", "irae3")))
  for( j in ind)
  {
    if( is.character(class(cov.dat[,j])))
    {
      cov.dat[,j] <- as.integer(as.factor(cov.dat[,j]))
    }
  }
  
  
  x <- as.matrix(cov.dat[,!(colnames(cov.dat) %in% c("Prior", "irae3", "Surv_Months", "Vital_Status_2yrs"))])
  y <- cov.dat$irae3
  
  print(paste0("run model with covariates; Prior ", colnames(x)))
  Prior <- cov.dat$Prior
  fit1 <- glm( y ~ x + Prior, family = "binomial")
  fit2 <- glm( y ~ geno.dat * Prior + x, family = "binomial")
  
  # interaction term will be last
  k <- dim(summary(fit2)$coefficients)[1]
  
  lrt.out <- lrtest(fit2,fit1)
  snp.b <- summary(fit2)$coefficients[2,1]
  snp.b.se <- summary(fit2)$coefficients[2,2]
  snp.p <- summary(fit2)$coefficients[2,4]
  int.b <- summary(fit2)$coefficients[k,1]
  int.b.se <- summary(fit2)$coefficients[k,2]
  int.p <- summary(fit2)$coefficients[k,4]
  df <-data.frame(snp.b = snp.b, 
                  snp.b.se = snp.b.se, 
                  snp.p = snp.p,
                  int.b = int.b,
                  int.b.se = int.b.se,
                  int.p = int.p,
                  lrt.p = lrt.out$`Pr(>Chisq)`[2])
  
  return(df)
}
# ------------------------------------------------------------------------------- #


# ---------------------------- coxPH ------------------------------------ #
coxPH <- function( geno.dat, cov.dat)
  # can be used for prior, no prior, or all subjects - subset geno.dat when passing to function
  # input:
  #   geno.dat (data.frame) - snp matrix of genotype data
  #   cov.dat (data.frame) -  covariate data that sepcifies phenotype and all other covariates
  #         the actual covariates are selected a priori and are hard-coded into their
  #         respective functions
  #
  # output:
  #   a table of two columns, snp name and corresponding p-value of the genotype term from the
  #   coxph model
{
  library(survival)
  # if MAF filter wasnt applied, this will throw an error
  if( sum(!is.na(unique(geno.dat))) == 1 )
  {
    return(NA)
  }
  
  if( all(is.na(geno.dat)) )
  {
    return(NA)
  }
  
  if( any(is.na(geno.dat)))
  {
    ind <- which(is.na(geno.dat))
    geno.dat <- geno.dat[-ind]
    cov.dat <- cov.dat[-ind,]
  }
  
  ind <- which(!(colnames(cov.dat) %in% c("PC1", "PC2", "PC3", "Prior", "irae3")))
  for( j in ind)
  {
    if( is.character(class(cov.dat[,j])))
    {
      cov.dat[,j] <- as.integer(as.factor(cov.dat[,j]))
    }
  }
  
  # cant have 0 survival time, but we know these patients died; assign an arbitrarily small number
  if( any(is.na(cov.dat$Surv_Months)))
  {
    cov.dat$Surv_Months[is.na(cov.dat$Surv_Months)] <- .1
  }
  
  if( any(cov.dat$Surv_Months == 0))
  {
    cov.dat$Surv_Months[ which(cov.dat$Surv_Months == 0)] <- .1
  }
  
  
  
  if( all(is.na(geno.dat)) )
  {
    return(NA)
  }
  
  if( any(is.na(geno.dat)))
  {
    ind <- which(is.na(geno.dat))
    geno.dat <- geno.dat[-ind]
    cov.dat <- cov.dat[-ind,]
  }
  
  # retain only descriptive covariates
  x <- as.matrix(cov.dat[,!(colnames(cov.dat) %in% c("Prior", "irae3", "Surv_Months", "Vital_Status_2yrs"))])
  y <- Surv(cov.dat$Surv_Months, cov.dat$Vital_Status_2yrs)
  print(paste0("run model with covariates; ", colnames(x)))
  fit <- coxph(y ~ geno.dat + x )
  b <- summary(fit)$coefficients[1,1]
  b.se <- summary(fit)$coefficients[1,3]
  p <- summary(fit)$coefficients[1,5]
  
  df <-data.frame(b = b, b.se = b.se, p = p)
  return(df)
}

# ------------------------------------------------------------------------------- #


# ---------------------------- coxPH.IpiNivo.priorint --------------------------- #
coxPH.priorint <- function( geno.dat, cov.dat )
  # input:
  #   geno.dat (data.frame) - snp matrix of genotype data
  #   cov.dat (data.frame) -  covariate data that sepcifies phenotype and all other covariates
  #         the actual covariates are selected a priori and are hard-coded into their
  #         respective functions
  #
  # output:
  #   a table of two columns, snp name and corresponding p-value of the likelihood ratio test,
  #   comparing a coxph model without genotype to a model with genotype * prior interaction term
  #   p will be significant if interaction or main effect of genotype are significant
{
  library(survival)
  library(lmtest)
  # if MAF filter wasnt applied, this will throw an error
  if( sum(!is.na(unique(geno.dat))) == 1 )
  {
    return(NA)
  }
  
  ind <- which(!(colnames(cov.dat) %in% c("PC1", "PC2", "PC3", "irae3")))
  for( j in ind)
  {
    if( is.character(class(cov.dat[,j])))
    {
      cov.dat[,j] <- as.integer(as.factor(cov.dat[,j]))
    }
  }
  
  # cant have 0 survival time, but we know these patients died; assign an arbitrarily small number
  if( any(is.na(cov.dat$Surv_Months)))
  {
    cov.dat$Surv_Months[is.na(cov.dat$Surv_Months)] <- .1
  }
  
  if( any(cov.dat$Surv_Months == 0))
  {
    cov.dat$Surv_Months[ which(cov.dat$Surv_Months == 0)] <- .1
  }
  
  if( all(is.na(geno.dat)) )
  {
    return(NA)
  }
  
  if( any(is.na(geno.dat)))
  {
    ind <- which(is.na(geno.dat))
    geno.dat <- geno.dat[-ind]
    cov.dat <- cov.dat[-ind,]
  }
  
 
  # retain only descriptive covariates
  x <- as.matrix(cov.dat[,!(colnames(cov.dat) %in% c("Prior", "irae3", "Surv_Months", "Vital_Status_2yrs"))])
  y <- Surv(cov.dat$Surv_Months, cov.dat$Vital_Status_2yrs)
  print(paste0("run model with covariates; Prior, ", colnames(x)))
 
  
  fit1 <- coxph(data = cov.dat, y ~ Prior + x )
  fit2 <- coxph(data = cov.dat, y ~ geno.dat * Prior + x)
  
  out <- lrtest(fit2,fit1)
  k <- dim(summary(fit2)$coefficients)[1]
  
  lrt.out <- lrtest(fit2,fit1)
  snp.b <- summary(fit2)$coefficients[2,1]
  snp.b.se <- summary(fit2)$coefficients[2,2]
  snp.p <- summary(fit2)$coefficients[2,4]
  int.b <- summary(fit2)$coefficients[k,1]
  int.b.se <- summary(fit2)$coefficients[k,2]
  int.p <- summary(fit2)$coefficients[k,4]
  df <-data.frame(snp.b = snp.b, 
                  snp.b.se = snp.b.se, 
                  snp.p = snp.p,
                  int.b = int.b,
                  int.b.se = int.b.se,
                  int.p = int.p,
                  lrt.p = lrt.out$`Pr(>Chisq)`[2])
  return(df)
}
# ------------------------------------------------------------------------------- #


# ------------------------------- attachMAF ------------------------------------- #
# function to attach minor allele frequencies to the data frame (used in meta-analysis)
# input: dat (data.frame), the results file
#        MAF (data.frame), the minor allele frequencies of each SNP
#
# output: dat (data.frame), results with MAF attached
attachMAF <- function( dat, MAF )
{
  if( substr(MAF$SNP[1], 1, 2) == "rs" )
  {
    ind <- match(dat$rsid, MAF$SNP)
  } else
  {
    ind <- match(dat$chrpos, MAF$SNP)
  }
  
  dat$MAF <- MAF$MAF[ind]
  return(dat)
}
# ------------------------------------------------------------------------------- #


# ------------------------------ runSNP.glm.parallel ---------------------------- #
# function to run the specified regression model for each SNP in the genotype data
# run on all available cores
# input:
#   geno.dat (data.frame) - snp matrix of genotype data
#   cov.dat (data.frame) -  covariate data that sepcifies phenotype and all other covariates
#         the actual covariates are selected a priori and are hard-coded into their
#         respective functions
#   BIM (data.frame) - the BIM file containing marker name information
#   model (string) - the function to call, one of: 
#     iraeAssoc
#     iraeAssoc.priorint
#     coxPH
#     coxPH.priorint
#   MAF (character) - the name of the file containing MAF information
#   outname (string) - name of the file where results are written
#
# output:
#   a table of two columns, snp name and corresponding p-value 
runSNP.glm.parallel <- function( geno.dat, cov.dat, BIM, model, MAF, outname )
{
  # assign function to call
  if( model == "iraeAssoc" )
  { func = iraeAssoc } else
    if( model == "iraeAssoc.priorint" )
    { func = iraeAssoc.priorint } else
      if( model == "coxPH" )
      { func = coxPH } else
        if( model == "coxPH.priorint" )
        { func = coxPH.priorint } else
        { 
          print(paste0("model ", model, " not an option."))
          stop()
        }
  
  # need to include all this info for pathway analysis
  # run glms in parallele
  out <- pbapply(geno.dat,   
                 2,              # apply to columns
                 func,           # called function
                 cov.dat,        # covariate data truncated to no prior subjects
                 cl = cl)                          # number of cores
  df <- matrix(0, nrow = length(out), ncol = length(out[[1]]))
  df <- as.data.frame(df)
  
  for( i in 1:length(out) )
  {
    df[i,] <- out[[i]]
  }
  
  colnames(df) <- colnames(out[[1]])
 
  # format and write data
  BIM <- BIM[ BIM$MarkerName %in% names(out),]
  BIM <- BIM[ match(names(out), BIM$MarkerName),]
  df$rsid <- BIM$MarkerName
  df$chr <- BIM$chr
  df$bp  <- BIM$bp
  df$chrpos <- paste(BIM$chr, BIM$bp, sep = ":")
  
  df$A1 <- BIM$A1
  df$A2 <- BIM$A2
  
  if(any(is.na(df$p)))
  {
    df <- df[ !is.na(df$p),]
  }
  
  MAF <- read.table(MAF, header = T)
  df <- attachMAF( df, MAF )
  write.table(df, outname, col.names = T, row.names = F, append = F, quote = F)
  
}
# ------------------------------------------------------------------------------- #


# =============================================================================== #
# =============================================================================== #
# =============================================================================== #






# =============================================================================== #
# ================================== MAIN ======================================= #
# =============================================================================== #


print(paste0("BEGIN: ", date()))

print("running with the following parameters:")
print(paste0("BEDFILE = ", BEDFILE))
print(paste0("COVARS = ", COVARS))
print(paste0("COVFILE = ", COVFILE))
print(paste0("CHR = ", CHR))
print(paste0("MAF = ", MAF))
print(paste0("OUTPREFIX = ", OUTPREFIX))
print("")

# the list of variables we will always consider, plus the additional specified covariates
if( COVARS == "NA")
{
  covar.names <- c("PC1", "PC2", "PC3", "irae3", "Prior", "Surv_Months", "Vital_Status_2yrs")
} else
{
  covar.names <- c(strsplit(COVARS, ",")[[1]], "PC1", "PC2", "PC3", "irae3", "Prior", "Surv_Months", "Vital_Status_2yrs")
}

if( !file.exists(BEDFILE))
{
  stop(paste0(BEDFILE, " not found. Exiting"))
}

if( !file.exists(COVFILE))
{
  stop(paste0(COVFILE, " not found. Exiting"))
}

# much more efficient  than trying to use the .raw file
print("reading input...")
geno.dat <- BEDMatrix(BEDFILE, simple_names = TRUE)
cov.dat <- read.table(COVFILE, header = T, sep = ",")

if(all(is.na(match(rownames(geno.dat), cov.dat$GWASID))))
{
  stop("could not match geno.dat subject ids to cov.dat subjects ids")
}

cov.dat <- cov.dat[match(rownames(geno.dat), cov.dat$GWASID),]
cov.dat <- cov.dat[ ,colnames(cov.dat) %in% covar.names,]
print("done")

print("setting up cores for parallel processing....")

# read BIM file for A1 and A2
BIM <- read.table( paste0( substr(BEDFILE, 1, nchar(BEDFILE) - 4), ".bim"), header = FALSE)
colnames(BIM) <- c("chr", "MarkerName", "GD", "bp", "A1", "A2")


# macOS workaround
cl <- parallel::makeCluster(detectCores(), setup_strategy = "sequential")
print("done")



# ------------- irae, no prior treatment ------------- #
if( MDL == 1)
{
  # make sure you specified the right covar names
  if(any(!(covar.names %in% colnames(cov.dat))))
  {
    ind <- which( !(covar.names %in% colnames(cov.dat)))
    print(paste0("ERROR: the following covar names were not found in covar data:", covar.names[ind]))
    stop()
  }
  
  # if irae, drop survival data
  cov.dat <- cov.dat[,!(colnames(cov.dat) %in% c("Surv_Months", "Vital_Status_2yrs"))]
  
  print("run glms for irae noprior model...")
  outname <- paste0(OUTPREFIX,"-chr", CHR, "-noprior.irae.txt")
  runSNP.glm.parallel( geno.dat[ cov.dat$Prior == 0,], 
                       cov.dat[ cov.dat$Prior == 0,], 
                       BIM,
                       "iraeAssoc",
                       MAF,
                       outname)
  print("done")
}
# ----------------------------------------------------- #


# ---------------- irae, prior treatment -------------- #
if( MDL == 2)
{
  print("run glms for irae prior-interaction model...")
  outname <- paste0(OUTPREFIX,"-chr", CHR, "-priorint.irae.txt")
  runSNP.glm.parallel( geno.dat, 
                       cov.dat, 
                       BIM,
                       "iraeAssoc.priorint",
                       MAF,
                       outname)
  print("done")
}
# ----------------------------------------------------- #


# ------------------ coxph no prior ------------------- #
if( MDL == 3)
{
  print("run glms for coxph no prior model...")
  outname <- paste0(OUTPREFIX,"-chr", CHR, "-noprior.coxph.txt")
  runSNP.glm.parallel( geno.dat[ cov.dat$Prior == 0,], 
                       cov.dat[ cov.dat$Prior == 0,], 
                       BIM,
                       "coxPH",
                       MAF,
                       outname)
  print("done")
}
# ----------------------------------------------------- #


# ----------------- coxph prior int ------------------- #
if( MDL == 4)
{
  print("run glms for coxph prior-interaction model...")
  outname <- paste0(OUTPREFIX,"-chr", CHR, "-priorint.coxph.txt")
  runSNP.glm.parallel( geno.dat, 
                       cov.dat, 
                       BIM,
                       "coxPH.priorint",
                       MAF,
                       outname)
}
# ----------------------------------------------------- #


# ----------------- irae ------------------- #
if( MDL == 5)
{
  # make sure you specified the right covar names
  if(any(!(covar.names %in% colnames(cov.dat))))
  {
    ind <- which( !(covar.names %in% colnames(cov.dat)))
    print(paste0("ERROR: the following covar names were not found in covar data:", covar.names[ind]))
    stop()
  }
  
  # if irae, drop survival data
  cov.dat <- cov.dat[,!(colnames(cov.dat) %in% c("Surv_Months", "Vital_Status_2yrs"))]
  
  print("run glms for irae noprior model...")
  outname <- paste0(OUTPREFIX,"-chr", CHR, "-all.irae.txt")
  runSNP.glm.parallel( geno.dat, 
                       cov.dat, 
                       BIM,
                       "iraeAssoc",
                       MAF,
                       outname)
  print("done")
}
# ----------------------------------------------------- #


print("stopping clusters...")
stopCluster(cl)
print("done")

print(paste0("analysis complete: ", date()))
# =============================================================================== #
# =============================================================================== #
# =============================================================================== #
