convertEnsemblToSymbol <- function( ensembl.genes )
{
  library("org.Hs.eg.db") # remember to install it if you don't have it already
  return(mapIds(org.Hs.eg.db, keys = ensembl.genes, keytype = "ENSEMBL", column="SYMBOL"))
}

convertSymbolToEnsembl <- function( symbol.genes )
{
  library("org.Hs.eg.db") # remember to install it if you don't have it already
  return(mapIds(org.Hs.eg.db, keys = symbol.genes, keytype = "SYMBOL", column="ENSEMBL"))
}
