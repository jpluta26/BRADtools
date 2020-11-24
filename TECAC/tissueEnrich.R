# pluta 11/20/20
#
# 1. modify the fuma output of tissue expression; truncate to genes discovered
# in TECAC replication analysis and put testis in the last column
#
# 2. perform tissue enrichment analysis
setwd("~/Documents/nathansonlab/tecac-manuscript/FUMA_gene2func30585/")

library(ggplot2)
library(dplyr)
library(reshape2)
library(SummarizedExperiment)



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



# ----- constants ----- #
# from replication analysis
novel.genes <- c("PPP2R5A", "BCL2L11", "TERT", "TNXB",
"BAK1", "DEPTOR","DMRT1","PSMB7","ANAPC2","ITIH5",
"ARL14EP","RAD52","METTL7A","SP1","CYTH1","ENOSF1",
"ZNF217","SUPT20HL1","AR","CENPI" ,"TKTL1")


rep.genes <- c("PMF1","UCK2","TFCP2L1","DAZL","TFDP2",
"SSR3", "TIPARP", "GPR160", "CDKL2", "G3BP2", "USO1",
"SMARCAD1", "HPGDS", "CENPE","TERT","CLPTM1L", "CATSPER3", 
"PITX1","SPRY4","BAK1", "GGNBP1","KATNA1","MAD1L1",
"PRDM14","DMRT1","GAB2", "USP35", "NARS2","PKNOX2","ATF7IP",
"KITLG","PRTG","ZWILCH","BCAR4",
"C16orf45", "MPV17L","HEATR3","RFWD3","ZFPM1","HNF1B",
"TEX14","ZNF257","LINC01859","RPSAP58", "ZNF726",
"ZNF254","NLRP12","ZFP64","MCM3AP","AIFM3")

all.genes <- unique(c(novel.genes, rep.genes))
# ----------- 3



# data from FUMA
dat <- read.table("gtex_v8_ts_general_avg_log2TPM_exp.txt", header = T)

# truncate to replication genes
dat2 <- dat[ dat$symbol %in% all.genes,]

# matrix of expression data
x = as.matrix(dat2[,3:dim(dat2)[2]])
rownames(x) <- dat2$symbol

# reorder columns so testis is last
k <- which(colnames(x) == "Testis")
x <- x[,reorderCols(x, k)]


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
        axis.text = element_text(size = 6)) +
  scale_fill_gradient2( low = "yellow", mid = "orange", high = "red", midpoint=3)
p1


png(height= 540, width = 675, "expr-plot.png")
print(p1)
dev.off()
# ----


# --- enrichment analysis ----
se <- SummarizedExperiment(assays = SimpleList(as.matrix(x)),
                         rowData = row.names(x),
                         colData = colnames(x))

# Ulhen et al 2015 algorithm suggests foldChangeThreshold = 5, but this
# is quite restrictive. use 2
genes <- assay(teGeneRetrieval(se, foldChangeThreshold = 2))

# select only genes tissue-enriched for testis
gene.set <- GeneSet(geneIds = genes[,1][genes[,2] == "Testis" & genes[,3] == "Tissue-Enriched"], 
                    organism = "Homo Sapiens" ,
                    geneIdType = SymbolIdentifier())

out = teEnrichment(inputGenes = gene.set, rnaSeqDataset = 1,
                   tissueSpecificGeneType = 2,  # tissue enriched
                   multiHypoCorrection = T, 
                   backgroundGenes = NULL)

seEnrichmentOutput <- out[[1]]
enrichmentOutput <- setNames(data.frame(assay(seEnrichmentOutput),
                                      row.names = rowData(seEnrichmentOutput)[,1]), 
                                      colData(seEnrichmentOutput)[,1])
enrichmentOutput$Tissue <- row.names(enrichmentOutput)
enrichmentOutput
