
#!/usr/local/bin/Rscript
rm(list = ls())

# pluta 6/26/23
# general script  to run all 4 kinds of ipi models (noprior - irae, priorint - irae, 
# noprior - coxph, priorint - coxph)
# writes output  with everything required for meta/pathway analysis

source("./iraeAssoc.R")

# ------------------------------------ user input  ------------------------------ #
suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c("-b", "--BEDFILE"), action="store", default = NA, type="character",
              help = "The plink .bed file containing the genotype data"),
  make_option(c("-c", "--COVARS"), action = "store", default = NA, type = "character",
              help = "The list of covariates, in quotes and comma separated. This must match columns in the phenotype file. \nExample: -c 'NDose.Nivo,ECOG,studyarm'"),
  make_option(c("-p", "--PHENOFILE"), action="store", default = NA, type="character",
              help = "The file containing phenotype data and covariates"),
  make_option(c("-C", "--CHR"), action="store", default = NA, type="integer",
              help = "The chromosome of interest."),
  make_option(c("-m", "--MODEL"), action="store", default = NA, type="integer",
              help = "1: irae, no prior\n
                      2: irae, prior interaction\n
                      3: coxph, no prior\n
                      4: coxph, prior interaction\n
                      5: all subjects, irae"),
  make_option(c("-M", "--MAF"), action="store", default = NA, type="character",
              help = "The .frq file from plink (not plink2)"),
  make_option(c("-f", "--MODELFORM"), action="store", default = NA, type="character",
              help = "the right-hand side of a glm formulation- this allows for the includsion of interaction terms beyond prior. Do not include genotype or prior, these are added automatically. \nExample: 'LDHBL + studyarm * Stage'."),
  make_option(c("-o", "--OUTPREFIX"), action="store", default = NA, type="character",
              help = "Name of output files.")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if( (opt$help))
{
  print_help(opt_parser)
  stop()
}

BEDFILE = opt$BEDFILE
COVARS = opt$COVARS
COVFILE = opt$PHENOFILE
CHR = opt$CHR
MDL=opt$MODEL
MAF=opt$MAF
MODELFORM=opt$MODELFORM
OUTPREFIX=opt$OUTPREFIX


# ------------------------------------------------------------------------------- #


# =============================================================================== #
# ================================== MAIN ======================================= #
# =============================================================================== #
print(paste0("BEGIN: ", date()))

print("running with the following parameters:")
print(paste0("BEDFILE = ", BEDFILE))
print(paste0("COVARS = ", COVARS))
print(paste0("COVFILE = ", COVFILE))
print(paste0("CHR = ", CHR))
print(paste0("MDL = ", MDL))
print(paste0("MAF = ", MAF))
print(paste0("MODELFORM = ", MODELFORM))
print(paste0("OUTPREFIX = ", OUTPREFIX))
print("")

# BEDFILE="nivo-chr6.qc2.bed"
# COVARS="studyarm,Stage,NDose.Nivo"
# COVFILE="nivo.pheno.txt"
# CHR=6
# MDL=2
# MAF="nivo-chr6.qc2.afreq"
# MODELFORM="studyarm * Stage + NDose.Nivo"

if( !file.exists(BEDFILE))
{
  stop(paste0(BEDFILE, " not found. Exiting"))
}

if( !file.exists(COVFILE))
{
  stop(paste0(COVFILE, " not found. Exiting"))
}

print("reading input...")
geno.dat <- BEDMatrix(BEDFILE, simple_names = TRUE)
geno.dat <- geno.dat[,1:100]
cov.dat <- read.table(COVFILE, header = T, sep = ",")
colnames(cov.dat)[2] <- "GWASID"
print("done")

print("matching genotype and covar ids")
geno.dat <- geno.dat[ rownames(geno.dat)  %in% cov.dat$GWASID,]
cov.dat <- cov.dat[match(rownames(geno.dat), cov.dat$GWASID),]

if( is.null(dim(geno.dat)))
{
  stop("geno.dat is empty, something went wrong.")
}

if(all(is.na(match(rownames(geno.dat), cov.dat$GWASID))))
{
  stop("could not match geno.dat subject ids to cov.dat subjects ids")
}

# the list of variables we will always consider, plus the additional specified covariates
if( COVARS == FALSE  )
{
  covar.names <- c("irae3", "Surv_Months", "Vital_Status_2yrs")
  
} else
{
  covar.names <- c(strsplit(COVARS, ",")[[1]], "irae3", "Prior", "Surv_Months", "Vital_Status_2yrs")
}

cov.dat <- cov.dat[ ,colnames(cov.dat) %in% covar.names,]
print("done")

print("setting up cores for parallel processing....")

# read BIM file for A1 and A2
BIM <- read.table( paste0( substr(BEDFILE, 1, nchar(BEDFILE) - 4), ".bim"), header = FALSE)
colnames(BIM) <- c("chr", "MarkerName", "GD", "bp", "A1", "A2")

# macOS workaround
cl <- parallel::makeCluster(detectCores(), setup_strategy = "sequential")
clusterExport(cl = cl, c("geno.dat", "cov.dat", "MODELFORM"))
clusterCall(cl = cl, fun = function() library("lmtest"))

print("done")

# simple case where  no coviarates are included
#  this overrides MDL variable
if( COVARS == FALSE  )
{
  print("running irae assoc simple...")
  out <- pbapply(geno.dat,   
                 2,              # apply to columns
                 snpAssocSimple,           # called function
                 cov.dat,        # covariate data truncated to no prior subjects
                 cl = cl)                          # number of cores
  df <- matrix(0, nrow = length(out), ncol = length(out[[1]]))
  df <- as.data.frame(df)
  
  for( i in 1:length(out) )
  {
    df[i,] <- out[[i]]
  }
  
  colnames(df) <- colnames(out[[1]])
  outname <- paste0(OUTPREFIX,"-chr", CHR, "-noprior.irae.txt")
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
  print('done"')
}



# ------------- irae, no prior treatment ------------- #
if( MDL == 1)
{
  # make sure you specified the right covar names
  if(COVARS != FALSE)
  {
    if(any(!(covar.names %in% colnames(cov.dat))))
    {
      ind <- which( !(covar.names %in% colnames(cov.dat)))
      print(paste0("ERROR: the following covar names were not found in covar data:", covar.names[ind]))
      stop()
    }
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
  if(COVARS != FALSE)
  {
    if(any(!(covar.names %in% colnames(cov.dat))))
    {
      ind <- which( !(covar.names %in% colnames(cov.dat)))
      print(paste0("ERROR: the following covar names were not found in covar data:", covar.names[ind]))
      stop()
    }
  }
  
  if( !exists("geno.dat"))
  {
    stop("prior to calling snpglmparallel, geno.dat doesnt exist")
  }
  print("run glms for irae prior-interaction model...")
  outname <- paste0(OUTPREFIX,"-chr", CHR, "-priorint.irae.txt")
  runSNP.glm.parallel( geno.dat, 
                       cov.dat, 
                       BIM,
                       "iraeAssoc.priorint",
                       MODELFORM,
                       MAF,
                       outname)
  print("done")
}
# ----------------------------------------------------- #


# ------------------ coxph no prior ------------------- #
if( MDL == 3)
{
  if(COVARS != FALSE)
  {
    if(any(!(covar.names %in% colnames(cov.dat))))
    {
      ind <- which( !(covar.names %in% colnames(cov.dat)))
      print(paste0("ERROR: the following covar names were not found in covar data:", covar.names[ind]))
      stop()
    }
  }
  
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
  if(COVARS != FALSE)
  {
    if(any(!(covar.names %in% colnames(cov.dat))))
    {
      ind <- which( !(covar.names %in% colnames(cov.dat)))
      print(paste0("ERROR: the following covar names were not found in covar data:", covar.names[ind]))
      stop()
    }
  }
  
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
  if(COVARS==FALSE)
  {
    if(any(!(covar.names %in% colnames(cov.dat))))
    {
      ind <- which( !(covar.names %in% colnames(cov.dat)))
      print(paste0("ERROR: the following covar names were not found in covar data:", covar.names[ind]))
      stop()
    }
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
