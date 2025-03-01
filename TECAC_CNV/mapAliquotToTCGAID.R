# pluta 8/12/21

# get TCGA ID from aliquot ID
library(TCGAutils)

files <- read.table("somatic.files.txt", header = FALSE)
aliquot.ids <- unlist(lapply(strsplit(files$V1, "[.]"), function(x) x[2]))

out <- UUIDtoBarcode(aliquot.ids, "aliquot_ids")
colnames(out) <- c("aliquot.id", "TCGA.id")

write.table(out, "aliquotIDmap.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)