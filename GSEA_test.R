rm(list = ls())


# ------------------------------------------------------- #
# PREPROCESSOR                                            #
# ------------------------------------------------------- #
setwd("C:/Users/jpluta/Desktop/Melanoma")
source("https://bioconductor.org/biocLite.R")

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(cpvSNP)
library(GSEABase)
library(ReportingTools)
library(BiocParallel)
library(GenomicRanges)
library(enrichplot)
library(DOSE)

# reference genes
genesHg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
 

# ------------------------------------------------------- #
# ------------------------------------------------------- #
# ------------------------------------------------------- #




# ------------------------------------------------------- #
# FUNCTIONS                                               #
# ------------------------------------------------------- #


# ------------------ make.GLOSSI.input ------------------ #
make.GLOSSI.input <- function( dat.fname )
  # function to create a GRange object that is used as the input for GLOSSI
  # meaning, it is a table of snps and p-values. SNPs are independent based on LD
  # and are named using rsID rather than chr:bp.
  # 
  # input: dat.fname - a file containing the GWAS data: snps and p-values
  #        ld.fname - file containing SNPs that are independent, according to LD
  #        ref.fname - file containing the mapping betweed rsid and chrbp
  #
  # output: a GRange object with all necessary data to run GLOSSI
{
  # these files shouldnt have trailing white space. fix this in bash: sed -i 's/[[:space:]]*$//' test-dat6.txt
  dat <- read.table(dat.fname, header = TRUE, as.is = TRUE)
  
  dat$chr <- unlist(lapply(strsplit(dat$MarkerName, ":"), function(x) x[1]))
  
  
  dat$bp  <- unlist(lapply(strsplit(dat$MarkerName, ":"), function(x) x[2]))
  dat$bp <- as.integer(dat$bp)
  
  g.arr <- data.frame(P = dat$P.value, SNP = dat$rsid, Position = dat$bp, chromosome = dat$chr)
  g.arr$Position <- as.integer(g.arr$Position)
  arrayDataGR <- suppressWarnings(createArrayData(g.arr, positionName = "Position"))
  arrayDataGR$chromosome <- droplevels(arrayDataGR$chromosome)
  arrayDataGR$SNP <- as.character(arrayDataGR$SNP)
  return( arrayDataGR )
  
}
# ------------------------------------------------------- #

# -------------------- runGLOSSI ------------------------ #
runGLOSSI <- function( snpsGSC, arrayDataGR )
# function to run GLOSSI; annotate snps with rsID and run statistical test
# input: arrayDataGR, GRange object containing p-vals and snp names
#        snpsGSC, the geneSets to run glossi on, annotated with rsID
#
# output: gRes, results of running glossi
{
  pvals <- arrayDataGR$P
  names(pvals) <- arrayDataGR$SNP
  
  gRes <- glossi(pvals, snpsGSC)
  
  return(gRes)
}
# ------------------------------------------------------- #

# -------------------- runGLOSSI ------------------------ #
runGLOSSIsequence <- function( snpsGSC, arrayDataGR )
  # with a large number of genesets, running GLOSSI may require too much memory
  # or take too long to run. in that case, run one at a time.
  # function to run GLOSSI; annotate snps with rsID and run statistical test
  # input: arrayDataGR, GRange object containing p-vals and snp names
  #        snpsGSC, the geneSets to run glossi on, annotated with rsID
  #
  # output: gRes, results of running glossi
{
  pvals <- arrayDataGR$P
  names(pvals) <- arrayDataGR$SNP
  gRes <- c()
  
  
  print("Running GLOSSI...")
  pb <- txtProgressBar( min = 0, max = length(snpsGSC), style = 3)
  
  for( i in 1:length(snpsGSC))
  {
    gRes <- c(gRes, glossi(pvals, snpsGSC[[i]]))
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  return(gRes)
}
# ------------------------------------------------------- #

# ------------------------------------------------------- #
expandGeneSet <- function( gset )
{
  g <- lapply(as.list(unique(geneIds(gset))), GeneSet)
  
  for( i in 1:length(unique(geneIds(gset))))
  {
    setName(g[[i]]) <- unique(geneIds(gset))[i]
  }
  
  
  return( GeneSetCollection(g) )
  
}
# ------------------------------------------------------- #

# ------------------- expandGeneSetCollection ----------- #
# function to take a GeneSetCollection, and convert each unique 
# gene into its own geneset
# essentially, running a gwas at the gene level
# 
# input: gsets, a full GeneSetCollection
# output: g, a GeneSetCollection where each unique gene in gsets
# is its own GeneSet
expandGeneSetCollection <- function ( gsets )
{
  g <- lapply(gsets, expandGeneSet)
  g <- lapply(unique(unlist(lapply(g, geneIds))), GeneSet)
  
  for( i in 1:length(unique(g)))
  {
    setName(g[[i]]) <- geneIds(g[[i]])
  }
  
  return( GeneSetCollection(g) )
}
# ------------------------------------------------------- #

# ------------------------------------------------------- #
# simple accessor function to get p-values from glossi results list
# need to add this because 'unlist' doesnt work on this type of object
getListPValues <- function( dat )
{
  return( dat@pValue )
}
# ------------------------------------------------------- #

# ------------------------------------------------------- #
# simple accessor function to get names from glossi results list
# need to add this because 'unlist' doesnt work on this type of object
getListNames <- function( dat )
# input: dat, an s4 object from a list
{

  return( dat@geneSetName )
}
# ------------------------------------------------------- #

# ------------------------------------------------------- #
createGeneList <- function( dat, pvals )
# input: dat, a single GeneSet of interest
#        pvals, the p-values from glossi on each gene individually
{
  g.id <- geneIds(dat)
  pvals <- sort(pvals[match(geneIds(dat), names(pvals))], decreasing = TRUE)
  return(pvals)
}
# ------------------------------------------------------- #


# ------------------------------------------------------- #
# END FUNCTIONS                                           #
# ------------------------------------------------------- #






# ------------------------------------------------------- #
# MAIN
# ------------------------------------------------------- #

# get geneset statistics
# 1. set input file
dat.fname = "bms2-glossi-input.txt"

# 2. define genesets
# gene sets
# oncogenic genesets
# gsets <- getGmt("c6.all.v6.2.entrez.gmt")

# immunicologic gene sets
# do geneset analysis first to find sig geneset
gsets <- getGmt("c7.all.v6.2.entrez.gmt")

# 3. read in snp data
arrayDataGR <- make.GLOSSI.input( dat.fname )

# 4. map snps to genes
snpsGSC <- geneToSNPList(gsets, arrayDataGR, genesHg19)

sig.thresh = 0.0001

# 5. run statistical testing
gRes1 <- runGLOSSI( snpsGSC, arrayDataGR )
sig.ind <- which(pValue( gRes1 ) < sig.thresh)

g.pvals <- unlist(lapply(gRes1, getListPValues))
names(g.pvals) <- unlist(lapply(gRes1, getListNames))
geneList <- sort(g.pvals, decreasing = T)

# analysis on genes within a gene set
# convert genesets into individual genes. eg each unique gene is a set
ref.genes.gsets <- expandGeneSetCollection( gsets )

snpsGSC <- geneToSNPList( ref.genes.gsets, arrayDataGR, genesHg19)
gRes.all <- runGLOSSIsequence( snpsGSC, arrayDataGR )
g.pvals <- unlist(lapply(gRes.all, getListPValues))
names(g.pvals) <- unlist(lapply(gRes.all, getListNames))
geneList <- createGeneList( gsets[[280]], g.pvals)


# run statistics within a geneset, from part 1
dgn <- gseDGN(geneList,
              nPerm         = 100, 
              minGSSize     = 120,
              pvalueCutoff  = 0.2, 
              pAdjustMethod = "BH",
              verbose       = FALSE)
dgn <- setReadable(dgn, 'org.Hs.eg.db')


ncg <- gseNCG(geneList,
              nPerm         = 1000, 
              minGSSize     = 120,
              pvalueCutoff  = 0.2, 
              pAdjustMethod = "BH",
              verbose       = FALSE)
ncg <- setReadable(ncg, 'org.Hs.eg.db')

do <- gseDO(geneList,
              nPerm         = 100, 
              minGSSize     = 120,
              pvalueCutoff  = 0.2, 
              pAdjustMethod = "BH",
              verbose       = FALSE)
do <- setReadable(do, 'org.Hs.eg.db')

# these only need the gene IDs, not any statistics
# explore the geneset?
# (DO ontology)
edo <- enrichDO(geneIds(gsets[[280]]), pvalueCutoff = 0.05)

# network of cancer genes
# the genes in this geneset are enriched in the following diseases from NCG
ncg <- enrichNCG(geneIds(gsets[[280]]))

# DisGeNET
dgn <- enrichDGN(geneIds(gsets[[280]]))

barplot(edo, showCategory = 10)
dotplot(dgn, showCategory = 10)




# write report
{pvals <- p.adjust( unlist(pValue(gRes1) ), method = "BH" )
report <- HTMLReport(shortName = "cpvSNP glossiRes",
                     title = "GLOSSI Results", reportDirectory = ".")
publish(chr22.gsets, report, annotation.db = "org.Hs.eg", setStats = unlist(statistic (gRes)),
        setPValues = pvals)
finish(report)


dat <- data.frame( p <- unlist(pValue(gRes1)), 
                   len <- unlist( lapply(geneIds(gsets), length) ), 
                   id <- names(gRes1))
colnames(dat) <- c("p", "len", "id")

# barplot
p1 <- ggplot( data = dat, aes( x = id, y = len, colour = p)) + geom_bar(stat = "identity") + coord_flip() +
  scale_color_gradient(low = "blue", high = "red", guide = "colourbar", limits = c(0,1), breaks = c(0,0.5,1))
p1

# dotplot
p2 <- ggplot( data = dat, aes( x = len, y = id, colour = p)) + 
  geom_point(aes(size = len)) + 
  scale_color_gradient(low = "blue", high = "red", guide = "colourbar", limits = c(0,1), breaks = c(0,0.5,1))
p2}
# ------------------------------------------------------- #
# ------------------------------------------------------- #
# ------------------------------------------------------- #
