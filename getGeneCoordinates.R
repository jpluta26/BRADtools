#!/usr/bin/env Rscript


# pluta 10/5/18
# given a gene name, this function returns the genomic coordits (geneID, start, and end)

# intended to be called from a bash script
# TODO: turn this into a pure R function
# input: GENE (string), a gene name
# output: gene.dat, a data.frame with fields geneID, chr, start, and end

args = commandArgs(trailingOnly = TRUE)

if( length(args) < 1 )
{
  stop("need to provide one argument: GENE, a string containing the gene name")
}


# ---------------------- preprocessor ---------------------------------------- #
source("https://bioconductor.org/biocLite.R")

if("TxDb.Hsapiens.UCSC.hg19.knownGene" %in% rownames(installed.packages()) == FALSE)
{
  biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
}

if("org.Hs.eg.db" %in% rownames(installed.packages()) == FALSE)
{
  biocLite("org.Hs.eg.db")
}

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
# ---------------------------------------------------------------------------- #


g = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

GENE = args[1]

gene.key <- select(org.Hs.eg.db, GENE, c("ENTREZID", "ALIAS"), keytype = "ALIAS" )
ind <- which(g$gene_id == gene.key$ENTREZID)

g.range <- g@ranges[ ind ]
chr <- g@seqnames[ ind ]@values

gene.dat <- data.frame( geneID = GENE, 
                   chr =    chr, 
                   start =  g.range@start, 
                   end =    g.range@start + g.range@width - 1)

return(gene.dat)
