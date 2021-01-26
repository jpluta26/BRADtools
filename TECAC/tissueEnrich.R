# pluta 1/26/21
#
# 1. modify the fuma output of tissue expression; truncate to genes discovered
# in TECAC replication analysis and put testis in the last column
#
# 2. perform tissue enrichment analysis
setwd("~/Documents/nathansonlab/tecac-manuscript/FUMA_gene2func30585/")


# load libraries 
library(ggplot2)
library(dplyr)
library(reshape2)
library(SummarizedExperiment)
library(TissueEnrich)


# ----------------------------- reorderCols --------------------------- #
reorderCols <- function( x, colPos )
  # move a specified column to last place in matrix
  
  # input: x, matrix of expression data
  #       colPos, the position of the column to be moved to last
  # output: pos, a vector of the new column order
{
  k <- dim(x)[2]
  pos <- 1:k
  pos[k] <- colPos
  pos[k - 1] <- k
  pos[((colPos + 1):(k - 1) - 1)] <- ((colPos + 1):(k - 1))
  return(pos)
}
# --------------------------------------------------------------------- #

# -------------------- getEnrichmentByScore ---------------------------- #
# calculate expression
# input: x, matrix of gene x tissue expression
# tissueSpecificType: 1 = all, 2 = tissue-enriched, 3 = tissue-enhanced, 4= group-enriched
# scores: range of scores to truncate by or NA
getEnrichmentByScore <- function(x, tissueSpecificGeneType, scores)
{
  if( sum(is.na(scores)) == 0)
  {
    x <- x[ rownames(x) %in% gene.info$gene[ gene.info$score %in% scores],]
  }
  
  
  se <- SummarizedExperiment(assays = SimpleList(as.matrix(x)),
                             rowData = row.names(x),
                             colData = colnames(x))
  
  # Ulhen et al 2015 algorithm suggests foldChangeThreshold = 5
  teGenes <- teGeneRetrieval(se, foldChangeThreshold = 5)
  genes <- assay(teGenes)
  # select only genes tissue-enriched for testis
  gene.set <- GeneSet(geneIds = genes[,1][genes[,2] == "Testis" & genes[,3] %in% c("Tissue-Enhanced", "Tissue-Enriched", "Group-Enriched")], 
                      organism = "Homo Sapiens" ,
                      geneIdType = SymbolIdentifier())
  
  out = teEnrichmentCustom(inputGenes = gene.set, 
                           tissueSpecificGenes = teGenes,
                           tissueSpecificGeneType = tissueSpecificGeneType,  
                           multiHypoCorrection = T, 
                           backgroundGenes = NULL)
  
  seEnrichmentOutput <- out[[1]]
  enrichmentOutput <- setNames(data.frame(assay(seEnrichmentOutput),
                                          row.names = rowData(seEnrichmentOutput)[,1]), 
                               colData(seEnrichmentOutput)[,1])
  enrichmentOutput$Tissue <- row.names(enrichmentOutput)
  
  # out[[2]]$Testis lists the geneset that was used
  return(list(enrichmentOutput, out[[2]]$Testis))
}
# --------------------------------------------------------------------- #


# --------------------------------------------------------------------- #
# ------------------------------ main --------------------------------- #
# --------------------------------------------------------------------- #

# gene scores from kate
gene.info <- read.table("~/Documents/nathansonlab/tecac-manuscript/genescores_012521.txt", header =T, sep = "\t")
gene.info  = gene.info[ gene.info$score >= 2,]

dup.ind <- which(duplicated(gene.info$gene))

# reassign the maximum score to all instances of repeated genes, then drop the duplicates
# eg only work with the highest score
for( i in dup.ind )
{
  gene = gene.info$gene[i]
  score <- max(gene.info$score[ which(gene.info$gene == gene)])
  gene.info$score[ gene.info$gene == gene ] <- score
}

gene.info <- gene.info[ -dup.ind,]

# data from FUMA
dat <- read.table("gtex_v8_ts_general_avg_log2TPM_exp.txt", header = T)
#notfound <- gene.info$gene[!(gene.info$gene %in% dat$symbol)]

# PACC1 = ENSG00000065600 -> set this to TMEM206
# MINP = BMERB1 = ENSG00000166780 -> set this to C16orf45
# HNF1B = ENSG00000275410.5 # not found
# LOC101927151 = ENSG00000267575 -> set this to CTC-459F4.3 
# ZNF217 = ENSG00000171940.13 # not found
# PDK3 = ENSG00000067992.16
# SUPT20HL1 = FAM48B1 = ENSG00000223731.3


gene.info$gene[which(gene.info$gene == "PACC1")] <- "TMEM206"
gene.info$gene[which(gene.info$gene == "BMERB1")] <- "C16orf45"
gene.info$gene[which(gene.info$gene == "LOC101927151")] <- "CTC-459F4.3"

# truncate to replication genes
X.genes <- c("AR", "CENPI", "TKTL1", "TEX28", "PDK3", "SUPT20HL1", "DRP2")
dat2 <- dat[ dat$symbol %in% c(gene.info$gene, X.genes),]


# matrix of expression data
x = as.matrix(dat2[,3:dim(dat2)[2]])
rownames(x) <- dat2$symbol

# reorder columns so testis is last
k <- which(colnames(x) == "Testis")
x <- x[,reorderCols(x, k)]

se <- SummarizedExperiment(assays = SimpleList(as.matrix(x)),
                           rowData = row.names(x),
                           colData = colnames(x))

# Ulhen et al 2015 algorithm suggests foldChangeThreshold = 5
# this gives expression status (enriched/enhanced/etc)
genes <- assay(teGeneRetrieval(se, foldChangeThreshold = 5))
genes <- genes[!duplicated(genes[,1]),]

# AMHR2 has group-enriched for several tissues, set to testis
genes[which(genes[,1] == "AMHR2"),2] <- "Testis"

# order by tissue and enrichment type
x <- x[order(genes[,2], genes[,3], decreasing = T),]




# ----
# recreate and modify plot from FUMA
m <- melt(x)

p1 <- ggplot(data= m, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() + 
  theme_minimal() + 
  xlab("Tissue") + 
  ylab("Gene") +
  guides(fill=guide_legend(title="Expression")) +
  theme(axis.text.x = element_text(angle=90, face = c(rep("plain", 29), "bold"), size = 10),
        axis.text = element_text(size = 8)) +
  scale_fill_gradient2( low = "yellow", 
                        mid = "orange", 
                        high = "red", midpoint=3, 
                        breaks = c(1,2,3,4,5),
                        labels = c("0-1", "2", "3", "4", "5-6")) 
p1


png(height= 740, width = 675, "expr-by-score_012521.png")
print(p1)
dev.off()
# ----


# --- enrichment analysis ----
# the chosen gene set is TEX14 DMRT1 DAZL ZNF728
# out[[2]]$Testis -> the selected genes

# -log10(p) = T
# p = 10^-T


enrich.all <- getEnrichmentByScore(x, 1, NA)
enrich.high <- getEnrichmentByScore(x, 1, gene.info$score[gene.info$score >= 3.0])
enrich.med <- getEnrichmentByScore(x, 1, c(2, 2.5))
# --------------------------------------------------------------------- #
# --------------------------------------------------------------------- #
# --------------------------------------------------------------------- #
