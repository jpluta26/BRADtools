# pluta 7/6/21
# script to create gene*expression/cnv/cpg methylation matrices of TCGA-TGCT data

library(IRanges)
library(data.table)
library('biomaRt')

# intended by run in /project/knathanslab/TECAC/TCGA





# $---------------------------------------------------------$ #
# $----------------------- functions -----------------------$ #
# $---------------------------------------------------------$ #

# -------------------- checkEmptyCol ---------------------- #
checkEmptyCol <- function( col )
#  check if a vector is all 0s
# input: col (integer/numeric), vector of values
# output: boolean
{
  if( all(col == 0))
  {
    return(TRUE)
  }
  
  return(FALSE)
}
# -------------------------------------------------------- #


# -------------------- removePath ------------------------ #
removePath <- function( filename )
# remove any path from a filename; filename itself contains subject id information used for
#   column naming
# input: filename (string), name of the input file with location
# output: filename (string), filename with path removed
{
  if(  grep("[/]", filename) == 1  )
  {
    x <- strsplit(filename, "[/]")
    return( x[[1]][length(x[[1]])])
  } else
  {
    return(filename)
  }
}
# -------------------------------------------------------- #



# -------------------------------------------------------- #
get.cpgMeth <- function( GENENAME, cpgmeth, hsapiens_genes )
# find all cpg methylation islands within a CNVR and return the average
# input:  
#   GENENAME (string), name of the gene of interest
#   cpgmeth (data.frame), cpg methylation data
# output:
#   (numeric), mean of all beta weights of cpg islands in the CNVR
#     or 0 if no islands found
{
  if(!(GENENAME  %in% hsapiens_genes$hgnc_symbol))
  {
    print(paste0(GENENAME, " not found"))
    return(NA)
  }
  
  gene <- hsapiens_genes[which(hsapiens_genes$hgnc_symbol == GENENAME),]
  
  cpgmeth <- cpgmeth[ which(cpgmeth$Chromosome %in% paste0("chr", gene$chromosome_name)),]
  cpgmeth <- cpgmeth[ !is.na(cpgmeth$Beta_value),]
  
  rangeA <- IRanges( cpgmeth$Start, cpgmeth$End )
  rangeB <- IRanges( gene$start_position, gene$end_position)
  
  ind <- findOverlaps(rangeA, rangeB, type = "within")
  if( length(ind) > 0 )
  {
    return(mean( cpgmeth$Beta_value[ind@from] ) )
  } else
    return(0)
}
# -------------------------------------------------------- #



# -------------------- readFPKMfile ---------------------- #
readFPKMfile <- function( filename, sample )
# read in a file of expression data and rename the columns to TCGA id
# input:
#   filename (string), name of the file to read (including path)
#   sample (data.frame), sample information and case id mappings
# output:
#   dat (data.frame), the expression file with columns labeled with TCGA id (TCGA-##-####)
{
  dat <- read.table(filename, header = F)
  
  colnames(dat) <- c("ENSG", sample$Case.ID[match( removePath(filename), sample$File.Name)])
  return(dat)
}
# -------------------------------------------------------- #

# $---------------------------------------------------------$ #
# $------------------- end functions -----------------------$ #
# $---------------------------------------------------------$ #



# ========================================================= #
# ========================= MAIN ========================== #
# ========================================================= #
args = commandArgs(trailingOnly=TRUE)

#  perform joint meta analysis of n datasets
if( length(args) < 1 )
{
  print("need to provide 1 arguments: GENELIST")
  print("GENELIST is a textfile with a list of gene names (eg. DMRT1) ")
  print("assumes expression data is in a subdirectory named expresion")
  print("methylation data is in a subdirectory named cpg_meth")
  print("all other files are in  root directory")
  print("")
  stop()
}


GENELIST <- args[1]
gene.list  <- read.table(GENELIST, header  = FALSE)$V1

print("preparing samples...")
# 150 cases total
sample <- read.table("gdc_sample_sheet.2021-07-06.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# for now, keep only unique subjects- remove secondary primary tumors
sample <- sample[ sample$Sample.Type == "Primary Tumor",]
sample$File.Name <- substr(sample$File.Name, 1, nchar(sample$File.Name) - 3)

#  should also filter based on clinical  data, eg race; but this is missing  16 subjects
clin <- read.table("clinical.tsv", header = TRUE, sep = "\t")
clin <- clin[ clin$race == "white",]
clin <- clin[ clin$ethnicity != "hispanic or latino",]

# remove nonwhite samples; this brings total down to 108?
# sample  <- sample[ sample$Case.ID  %in% clin$case_submitter_id, ]
print("done")


print("reading in expression data...")
# TCGA  - transcriptome profiling - *FPKM.txt.gz
files  <- list.files(path  = "expression", pattern = "*.FPKM.txt")
files <- files[ files %in% sample$File.Name ]
files <- paste0("expression/", files)


# ---- get expression data -----
# read and merge all subjects expression data
for( file in files )
{
  if( file ==  files[1] )
  {
    expr <- readFPKMfile( file,  sample )
  } else 
  {
    expr <-  merge(expr, readFPKMfile( file,  sample), by = "ENSG")
  }
}

expr <- as.data.frame(expr)
print("done")


print("mapping gene name to ensembl id...")
# map ensembl symbol to gene name
hsapiens_genes <- getBM(attributes = c("ensembl_gene_id", 
                                       "hgnc_symbol", "start_position", "end_position", "chromosome_name"),
                        mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "useast.ensembl.org"))

hsapiens_genes <- hsapiens_genes[ hsapiens_genes$chromosome_name %in% as.character(seq(1:22)),]

# add gene name and standardize ensembl name
expr$ENSG <- as.character(expr$ENSG)
expr$ENSG <- unlist(lapply(strsplit(expr$ENSG, "[.]"), function(x) x[1]))
expr$GeneName <- hsapiens_genes$hgnc_symbol[match(expr$ENSG, hsapiens_genes$ensembl_gene_id)]
expr <- expr[ expr$GeneName %in% gene.list, ]
print("done")


print("reading in cpg methylation  data...")
# create an empty matrix with the same rows and columns as expr, for storing cpg.meth data
cpgmeth <- matrix(0, nrow = dim(expr)[1], ncol = dim(expr)[2])
colnames(cpgmeth) <- colnames(expr)
rownames(cpgmeth) <- rownames(expr)
cpgmeth  <- as.data.frame(cpgmeth)
cpgmeth$GeneName <- expr$GeneName
cpgmeth$ENSG <- expr$ENSG
# ---------- #


# ------ get cpg methylation data ------ #
files <- list.files("cpg_meth/")
cpgmeth <- cpgmeth[ cpgmeth$GeneName %in% gene.list,]
k <- dim(cpgmeth)[2]

# ignore first and last column (gene names)
for( j in 2:(k-1))
{
  sub.id <- colnames(cpgmeth)[j]
  if( length(grep(sub.id, files)) == 0 )
  {
    print(paste0("no file found for subject ", sub.id))
  } else {
      filename <- paste0( "cpg_meth/", files[grep(sub.id, files)])
      cpgmeth.file <- fread(filename, header = TRUE, sep = "\t")
      
        # iterate through each gene and get the mean methylation beta weight
        for( i in 1:length(cpgmeth$GeneName))
        {
          cpgmeth[i,j] <- get.cpgMeth( cpgmeth$GeneName[i], cpgmeth.file, hsapiens_genes)
        }
      }
}

# remove subjects with no methylation data
ind <- !apply(cpgmeth,2,checkEmptyCol)
cpgmeth <- cpgmeth[,ind]
expr <- expr[,ind]

if(!all(colnames(cpgmeth) == colnames(expr)))
{
  stop("subject mismatch between cpgmeth and expr")
}

print("done")


print("reading CNV data...")
# cnv data
genecnv  <- read.table("all_wgs_all_thresholded.by_genes.txt", header = TRUE,  sep = "\t")
genecnv <- genecnv[ ,-c(2:3)]

colnames(genecnv) <- substr(colnames(genecnv), 1, 12)
colnames(genecnv) <- gsub("[.]", "-", colnames(genecnv))
colnames(genecnv)[1] <- "GeneName"
print("done")

# reduce cnv, expression, and methylation data to common set of subjects
common.subs <- intersect(colnames(genecnv), intersect(colnames(expr), colnames(cpgmeth)))
genecnv <- genecnv[ genecnv$GeneName %in% gene.list,]

genecnv <- genecnv[ ,colnames(genecnv) %in% common.subs]
expr  <- expr[,colnames(expr) %in% common.subs]
cpgmeth  <- cpgmeth[,colnames(cpgmeth) %in% common.subs]
genecnv <- genecnv[,match(colnames(expr), colnames(genecnv) )]

# make sure columns are in the same  order
if(!all(colnames(genecnv) == colnames(expr)))
{
  stop("subject mismatch between cpgmeth and expr")
}

print("writing data...")
# write output to be used in cnv - expression testing
write.table(expr, "TCGAexpressionMatrix.txt", col.names =  TRUE, row.names = FALSE, quote = F)
write.table(genecnv, "TCGAcnvMatrix.txt", col.names = TRUE, row.names = FALSE, quote  = FALSE)
write.table(cpgmeth, "TCGAcpgmethMatrix.txt",  col.names = TRUE, row.names = FALSE, quote  = FALSE)
print("done")
# ========================================================= #
# ========================================================= #
# ========================================================= #
