# pluta 4/9/21

# set of functions for analysis of Ipi/Nivo irAE data
# efficient modeling through the runSNP.glm.parallel function
# provides modeling of the data for all 4 model types:
# irae - no prior; irae - prior interaction; coxph - no prior; coxph - prior interaction
# intended to be invoked by a wrapper script


# ---------------------------  preprocessor ------------------------------------- #
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
library(lmtest)
library(speedglm)
# ------------------------------------------------------------------------------- #





# =============================================================================== #
# ================================= functions =================================== #
# =============================================================================== #


# -------------------------------- snpAssocSimple ------------------------------------ #
snpAssocSimple <- function( geno.dat, cov.dat )
# association of  genotype and irae, with no covariates
{
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
  
  fit <- glm(cov.dat$irae3 ~ geno.dat, family  = "binomial")
  
  b <- summary(fit)$coefficients[2,1]
  b.se <- summary(fit)$coefficients[2,2]
  p <- summary(fit)$coefficients[2,4]
  
  df <-data.frame(b = b, b.se = b.se, p = p)
  colnames(df) <- c("b", "b.se", "b.p")
  return(df)
}
# ----------------------------------------------------------------------------------- #


# -------------------------------- iraeAssoc ------------------------------------ #
iraeAssoc <- function( geno.dat, cov.dat )
  # input:
  #   geno.dat (data.frame) - snp matrix of genotype data
  #   cov.dat (data.frame) -  covariate data that contains phenotype and all other covariates
  #
  # output:
  #   a table of two columns, snp name and corresponding p-value of the genotype term from
  #   the logistic regression model
{
  
  # if there are no minor alleles, return NA
  # run MAF filter to address this
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

  
  # convert categorical variables from factors to integers, so that they may be included
  # in the design matrix
  ind <- which(!(colnames(cov.dat) %in% c("PC1", "PC2", "PC3", "Prior", "irae3")))
  for( j in ind)
  {
    if( is.character(class(cov.dat[,j])))
    {
     cov.dat[,j] <- as.integer(as.factor(cov.dat[,j]))
    }
  }
  
  # pack all covariates into a single design matrix; exclude any covariates not included in the model
  x <- as.matrix(cov.dat[,!(colnames(cov.dat) %in% c("Prior", "irae3", "Surv_Months", "Vital_Status_2yrs"))])
  y <- cov.dat$irae3
  
  # print(paste0("run model with covariates; ", colnames(x)))
  # have to pass the whole df and specify the covariates of interest
  fit  <- glm(y ~ geno.dat + x , family = "binomial")
  
  # stats for the genotype term
  b <- summary(fit)$coefficients[2,1]
  b.se <- summary(fit)$coefficients[2,2]
  p <- summary(fit)$coefficients[2,4]
  
  df <-data.frame(b = b, b.se = b.se, p = p)
  colnames(df) <- c("b", "b.se", "b.p")
  return(df)
}
# ------------------------------------------------------------------------------- #


# --------------------------- iraeAssoc.priorint -------------------------------- #
iraeAssoc.priorint <- function( geno.dat, cov.dat, fit.only = FALSE)
# input:
#   geno.dat (data.frame) - snp matrix of genotype data
#   cov.dat (data.frame) -  covariate data that contains phenotype and all other covariates
#   fit.only (boolean) - flag to return just the fit of the prior interaction model, rather
#       than summary statistics. used for meta analysis.  
  
# output:
#   df (data.frame) - summary statistics from the interaction model
#   p will be significant if interaction or main effect of genotype are significant
#
{
  # when using multiple cores, the environment doesnt copy to each core
  # need to call libraries within function
  
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
  
  # print(paste0("run model with covariates; Prior ", colnames(x)))
  Prior <- cov.dat$Prior
  
  fit1 <- speedglm( y ~ x + Prior, family = binomial(logit))
  fit2 <- speedglm( y ~ geno.dat * Prior + x, family = binomial(logit))
  
  if(fit.only ==  TRUE)
  {
    return(fit2)
  }
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
  # print(paste0("run model with covariates; ", colnames(x)))
  fit <- coxph(y ~ geno.dat + x )
  b <- summary(fit)$coefficients[1,1]
  b.se <- summary(fit)$coefficients[1,3]
  p <- summary(fit)$coefficients[1,5]
  
  df <-data.frame(b = b, b.se = b.se, p = p)
  return(df)
}

# ------------------------------------------------------------------------------- #


# ---------------------------- coxPH.priorint --------------------------- #
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
  # print(paste0("run model with covariates; Prior, ", colnames(x)))
 
  
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
#
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
#
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
  # run glms in parallel
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



