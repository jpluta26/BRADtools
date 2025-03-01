# pluta 6/30/21
# v1.0 9/8/21

# run in dir:
# /Users/johnpluta/Documents/nathansonlab/CNV/final

# ===================================================================================================== #
# ============================================= functions ============================================= #
# ===================================================================================================== #

# --------------------------------------- geneRanges -------------------------------------------------- #
geneRanges <- function( db , column = "ENTREZID")
{
    g <- genes(db, columns=column)
    col <- mcols(g)[[column]]
    genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
    mcols(genes)[[column]] <- as.character(unlist(col))
    genes
}
# ----------------------------------------------------------------------------------------------------- #

# --------------------------------------- splitByOverlap ---------------------------------------------- #
splitByOverlap <- function(query, subject, column="ENTREZID", ...)
{
    olaps <- findOverlaps(query, subject, ...)
    f1 <- factor(subjectHits(olaps),
                 levels=seq_len(subjectLength(olaps)))
    splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)
}
# ----------------------------------------------------------------------------------------------------- #


# ---------------------------------------------- getCNVRAssoc ------------------------------------------ #
getCNVRAssoc <- function(CNVR.name, M, pheno)
# function to get statistics for association between CNV state and phenotype, w/ adjustment for ancestry
# input: CNVR.name  (string), name of the CNVR
#        M (matrix), design matrix of CNVR values
#        pheno (data.frame), df containing phenotype and covariate data
# output: out, summary statistics from the association test
{
  if( sum(pheno$SSID == M$IID) != dim(M)[1])
  {
    print("pheno/M id mismatch")
    stop("exiting")
  }
  
  dat <- M[,colnames(M) == CNVR.name]
  out <- data.frame(CNVR = CNVR.name, B = NA, B.se = NA, p = NA)
  
  # if all CNVRs values are the same, no association test can be performed
  if( dim(table(dat))  ==  1)
  {
    return(out)
  } else
  {
    # cant use site as a covariate, observations are too sparse
    fit <- glm(data = pheno, PHENO ~ as.integer(dat) + EV1 + EV2 + EV3 , family = "binomial")
    out$B <- summary(fit)$coefficients[2,1]
    out$B.se <- summary(fit)$coefficients[2,2]
    out$p <- summary(fit)$coefficients[2,4]
    return(out)
  }
}
# ----------------------------------------------------------------------------------------------------- #



# ------------------------------------- createDesignMatrix -------------------------------------------- #
createDesignMatrix <- function(cnvr.name, cnvr, dat, subject.names)
# create the CNVR design matrix; intended to  be used with apply
# input:
#   cnvr.name (string), name of the CNVR of interest (eg 'CNVR_101')
#   cnvr (data.frame), data.frame of cnvr regions
#   dat (data.frame), data.frame of cnv regions and values
#   subject.names (string), the vector of  unique subject IDs
#
# output:
#   out (data.frame), matrix of CNV_values mapped to CNVR
{
  
  CNVR <- cnvr[ cnvr$CNVR_ID == cnvr.name, ]
  
  chr <- dat[ dat$Chr == CNVR[["Chr"]], ]
  
  # map  cnvs to cnvrs
  cnvr.range <- IRanges( as.integer(CNVR[["Start"]]), as.integer(CNVR[["End"]]))
  sub.range <- IRanges(chr$Start, chr$End)
  ind <- findOverlaps(sub.range, cnvr.range, type = "within")
  
  subInCNVR <- chr[ind@from,]
  ind <- match(subject.names, subInCNVR$Sample_ID )
  
  out <- rep(2, length(subject.names))
  out[!is.na(ind)] <- subInCNVR$CNV_Value[ind[!is.na(ind)]]
  
  return(out)
}
# ----------------------------------------------------------------------------------------------------- #

# ===================================================================================================== #
# ===================================================================================================== #
# ===================================================================================================== #





# ===================================================================================================== #
# ============================================== MAIN ================================================= #
# ===================================================================================================== #


# the CNV file with qc already performed; no further qc steps taken 


# this file is created from penncnv output after qc. HandyCNV expects 8 columns, which is what is produced
# from PennCNV1.04. PennCNV1.05 produces 7 columns. this file has a dummy column added so that it can be read,
# and 'cases'  replace with 'controls' so that penn_id_sep is uniform
# more penncnvINPUT.rawcnv | awk '{print $0, "conf=1"} > penncnv_qc.rawcnv
#  sed -i 's/cases/controls/g' penncnv_qc.rawcnv

# if you  get object 'V2' not found error - run dos2unix on rawcnv file
# and make sure space seperated, not tab:  more file.rawcnv | tr "\t" " " > newfile.rawcnv
# penn_id_sep needs to be uniform- cant have cases and controls; fix this again with tr

# filter by amplification/deletion before running the script
args = commandArgs(trailingOnly =  TRUE)
if( length(args) < 6 )
{
  print("USAGE:: ")
  print("runCNVRAssoc.R RAWCNVFILE IDPREFIX OUTDIR MANPLOTNAME RESNAME DUPDEL")
  print("RAWCNVFILE = the .rawcnv file to be used, that has been split into only dup or del")
  print("IDPREFIX = the prefix before each subject id in the raw CNV file")
  print("OUTDIR = output directory")
  print("MANPLOTNAME = filename for the manhattan plot")
  print("RESNAME = filename for the results file")
  print("DUPDEL = run only on duplications or deletions (TRUE) or all data together (FALSE)")
  print("")
  print("if cnv_clean returns an error, run prepRawCNV.sh on it first")
  stop()
}

RAWCNVFILE  = args[1]
IDPREFIX    = args[2]
OUTDIR      = args[3]
MANPLOTNAME = args[4]
RESNAME     = args[5]
DUPDEL      = args[6]



library(HandyCNV)
library(rtracklayer)
library(Homo.sapiens)
library(parallel)
library(pbapply)
library(dplyr)
library(data.table)
library(qqman)

# modified version of call_cnvr that accomodates dup/del separately
source("~/TECAC_CNV/call_cnvr2.R")






cnv.dat <- cnv_clean(penncnv  = RAWCNVFILE,
                     penn_id_sep = IDPREFIX)


# -- format phenotype and covariate data --
pheno <- read.table("TECAC_CNV_PHENO.csv", header = TRUE, sep = ",")

pheno <- pheno[ pheno$BID %in%  cnv.dat$Sample_ID,]
pheno <- pheno[ match(unique(cnv.dat$Sample_ID), pheno$BID),]

evec <- read.table("tecac.evec", header = F)
evec <- evec[ evec$V1 %in% pheno$SSID,]
evec <- evec[ match(pheno$SSID, evec$V1),]
pheno$EV1 <- evec$V2
pheno$EV2 <- evec$V3
pheno$EV3 <- evec$V4
pheno$site <- substr(pheno$SSID, 1, 2)

if( all(pheno$PHENO %in% c(1,2)))
{
  pheno$PHENO <- pheno$PHENO - 1
}
# --


# convert CNVs to CNVRs; if desired output is only duplication or deletion,
# use a modified version call_cnvr
if( DUPDEL == TRUE )
{
  cnvr <- call_cnvr2(clean_cnv = "cnv_clean/penncnv_clean.cnv",
                     chr_set = 22,
                     folder = OUTDIR)
} else
{
  cnvr <- call_cnvr(clean_cnv = "cnv_clean/penncnv_clean.cnv",
                     chr_set = 22,
                     folder = OUTDIR)
}


# setup the design matrix
# this is pretty fast, dont need to  paralellize
# M is the matrix of CNVR x subject with cells 0 (double deletion),  1 (single deletionn),
# 2 (no change), 3 (duplication), 4 (double duplication)
M <- lapply(cnvr$CNVR_ID, createDesignMatrix, cnvr, cnv.dat, unique(cnv.dat$Sample_ID))
M <- do.call(cbind.data.frame, M)
colnames(M) <- cnvr$CNVR_ID


M$IID <- unique(cnv.dat$Sample_ID)


map <- read.table("IKN_TECAC_SSID_updated.txt", header = TRUE)
M$IID <- map$SSID[match(M$IID, map$BID)]

# slight correction to match genotype data
# M$IID[grep("C362", M$IID)]  <- "UK0000C362F"
# M$IID[grep("C310", M$IID)]  <- "UK0000C310F"

# assumes the naming convention of filename_fix.dup.rawcnv, probably should generalize this
write.table(M, paste0(strsplit(RAWCNVFILE, "[.]")[[1]][2], "_cnvr.mat"), col.names = TRUE, row.names = FALSE, quote = FALSE)



# run association testing
cl <- parallel::makeCluster(detectCores(), setup_strategy = "sequential")

# dont run  on the IID column
k <- dim(M)[2]
out <- pblapply(colnames(M)[1:(k-1)], getCNVRAssoc, M, pheno)
out <- do.call(rbind.data.frame, out)

stopCluster(cl)


# remove NA results and format output df
# could also threshold on total number of observations, eg. frequency >= 100
out  <- out[ !is.na(out$p),]
cnvr <- cnvr[ cnvr$CNVR_ID %in%  out$CNVR  ,]
cnvr <- cnvr[match(out$CNVR, cnvr$CNVR_ID),]
out$Chr <- cnvr$Chr
out$Start <- cnvr$Start
out$End <- cnvr$End

if( substr(MANPLOTNAME, nchar(MANPLOTNAME) - 2, nchar(MANPLOTNAME)) != "png" )
{
  MANPLOTNAME <- paste0( MANPLOTNAME, ".png")
}
png(MANPLOTNAME, width = 800, height = 500)
manhattan( out, chr = "Chr", bp = "Start", p = "p", snp = "CNVR", logp = TRUE)
dev.off()

# maps genes to intervals
call_gene(refgene = "refgene/Human_hg19.txt",
          interval =  paste0(OUTDIR, "/cnvr.txt"),
          clean_cnv = "cnv_clean/penncnv_clean.cnv",
          folder = OUTDIR)

# another way of doing the above
g <- GRanges(seqnames = paste0("chr", out$Chr), ranges = IRanges(start = out$Start, end = out$End))
gns = geneRanges(Homo.sapiens, column="SYMBOL")
symInCnv = splitByOverlap(gns, g, "SYMBOL")




# merge association statistics with CNVR data
out2 <- merge(cnvr, out, by.x = 'CNVR_ID', by.y  = 'CNVR')

# cant  select by column name? whats going on here
out2 <- out2[,-c(16:18)]


colnames(out2)[which(colnames(out2) == "Chr.x")] <- "Chr"
colnames(out2)[which(colnames(out2) == "Start.x")] <- "Start"
colnames(out2)[which(colnames(out2) == "End.x")] <- "End"

#sig.res <- out2[ which(out2$p < 1e-05), ]

write.table(out2, RESNAME, col.names = TRUE, row.names = F, quote = F)

# ===================================================================================================== #
# ===================================================================================================== #
# ===================================================================================================== #
