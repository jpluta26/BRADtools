#!/usr/bin/env Rscript

# pluta 3/13/19
# script to convert gentoype data from the CAG to plink format


args = commandArgs( trailingOnly = TRUE )
if( length(args) < 2 )
{
        print("need to provide 2 arguments: INFILE, OUTNAME")
        print("INFILE is the name of the csv file of genotype data")
        print("OUTNAME is the rootname of the PLINK files to be created")
        stop()
}

INFILE = args[1]
OUTNAME = args[2]

# --------------------------------------------------------------------------- #
# ------------------------------------- functions --------------------------- #
# --------------------------------------------------------------------------- #

# --------------------------- createPlinkFiles ------------------------ #
createPlinkFiles <- function( dat, snpdat, OUTNAME )
  # function to convert CAG genotype calls to PLINK format
  # input: dat (data.frame), the genotype data. rows are subjects, columns are snps.
  #           genotypes are character strings of the form "A:C", or "No Call" for missing data
  #        snpdat (data.frame), genotype calls are in rsid; need to convert this to chrpos
  #           snpdat is the data providing this mapping
  #        OUTNAME (string), the root name of the plink files
  # output: OUTNAME.ped/map/fam
{
  # character matrix for genotype data
  # first 6 columns are FID IID MotherID FatherID Sex Pheno
  p <- 6 + (dim(dat)[2] - 1) * 2
  ped <- matrix("0", nrow = dim(dat)[1], ncol = p)

  # convert failed calls to 0 and expand genotype calls to two characters
  for( i in 1:dim(dat)[1])
  {
    if( any( dat[i,] %in% c("Invalid", "No Call")))
    {
      dat[i,which(dat[i,] %in% c("Invalid", "No Call"))] <- "0:0"
    }

    ped[i,7:p ] <- unlist( strsplit( as.character(dat[i,2:dim(dat)[2]]), ":" ))
  }

  # iid
  ped[,2] <- dat$id
 # pheno
  ped[,6] <- "-9"

  # create map file from ref data
  rsids <- colnames(dat[2:dim(dat)[2]])
  snpdat <- snpdat[ match(snpdat$rsid, rsids), ]
  map <- data.frame( chr <- snpdat$chr, rsid <- paste(snpdat$chr, snpdat$bp, sep = ":"),  gd <- 0, bp <- snpdat$bp)

  fam <- ped[,1:6]

  PEDNAME <- paste( OUTNAME, "ped", sep = ".")
  FAMNAME <- paste( OUTNAME, "fam", sep = ".")
  MAPNAME <- paste( OUTNAME, "map", sep = ".")

  write.table( ped, PEDNAME, quote = FALSE, row.names = FALSE, col.names = FALSE, append = FALSE, sep = "\t")
  write.table( fam, FAMNAME, quote = FALSE, row.names = FALSE, col.names = FALSE, append = FALSE, sep = "\t")
  write.table( map, MAPNAME, quote = FALSE, row.names = FALSE, col.names = FALSE, append = FALSE, sep = "\t")

}
# ----------------------------------------------------------------------- #


# -------------------------------------------------------------------------- #
# -------------------------------------------------------------------------- #
# -------------------------------------------------------------------------- #



# --------------------------------- MAIN -------------------------------------- #
# reference data- chrpos mapped to rsid
snpdat <- read.table("snpdat.txt", header = TRUE, as.is = TRUE)

# the number of lines to skip is specific to this dataset
dat <- read.table( INFILE, sep = ",", skip = 115, header = TRUE, as.is = TRUE, colClasses = rep("character", 98))
colnames(dat)[2] <- "id"

# drop control runs
if( length( grep("NTC", dat$id)) != 0 )
{
        dat <- dat[ -grep("NTC", dat$id),]
}

# 48 snps but they were sampled twice; divide this into two datasets
x = substr(colnames(dat), nchar(colnames(dat)) - 1, nchar(colnames(dat)))

# screen for duplicates
dat1 <- dat[,(x == ".1")]
dat1$id <- dat$id
dat1 <- subset(dat1, select = c(49, 1:48))

dat2 <- dat[,(x != ".1")]
dat2 <- dat2[,-c(1:2)]
dat2$id <- dat$id
dat2 <- subset(dat2, select = c(49, 1:48))

# get the snp names correct
colnames(dat1) <- gsub("\\.1", "", colnames(dat1))
colnames(dat1)[1] <- "id"
colnames(dat1) <- lapply(strsplit(colnames(dat1), "_"), function(l) l[[1]])

colnames(dat2)[1] <- "id"
colnames(dat2) <- lapply(strsplit(colnames(dat2), "_"), function(l) l[[1]])

if( any(!(colnames(dat1) == colnames(dat2))))
{
  stop("colname mismatch")
}


# write out the two different plink files to compare
# eg were there two runs because of genotype fails?
OUTNAME1 = paste(OUTNAME, "1", sep = ".")
OUTNAME2 = paste(OUTNAME, "2", sep = ".")
createPlinkFiles( dat1, snpdat, OUTNAME1 )
createPlinkFiles( dat2, snpdat, OUTNAME2 )

# percent concordance by snp and sub
x <-  dat1[,2:49] == dat2[,2:49]
snp.concordance <- colSums(x) / dim(dat1)[1]
sub.concordance <- rowSums(x) / dim(dat1)[2]

snp.out <- data.frame( rsid = colnames(dat1)[2:49], concordance = snp.concordance)
sub.out <- data.frame( id = dat1$id, concordance = sub.concordance)

OUTNAMESNP <- paste(OUTNAME, "snp.concordance", sep = ".")
OUTNAMESUB <- paste(OUTNAME, "sub.concordance", sep = ".")

write.table( snp.out, OUTNAMESNP, col.names = TRUE, row.names = FALSE, quote = FALSE, append = FALSE)
write.table( sub.out, OUTNAMESUB, col.names = TRUE, row.names = FALSE, quote = FALSE, append = FALSE)
