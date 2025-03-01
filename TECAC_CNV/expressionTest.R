setwd("~/Documents/nathansonlab/CNV/TCGA")


expr <- read.table("TCGAexpressionMatrix.txt", header = TRUE, sep = " ")
colnames(expr) <- gsub("[.]", "-", colnames(expr))
# sample <- read.table("gdc_sample_sheet.2021-07-06.tsv", header = TRUE, sep = "\t")
# 
# # for now, keep only unique subjects- remove secondary primary tumors
# sample <- sample[ sample$Sample.Type == "Primary Tumor",]
# sample$File.Name <- substr(sample$File.Name, 1, nchar(sample$File.Name) - 3)
# sample$File.Name <- substr(sample$File.Name, 1, nchar(sample$File.Name) - 9)
# ind <- which(substr(sample$File.Name, 1, 1) %in% seq(1:9))
# sample$File.Name[ind] <- paste0("X", sample$File.Name[ind])

genecnv  <- read.table("all_wgs_all_thresholded.by_genes.txt", header = TRUE,  sep = "\t")
genecnv <- genecnv[ ,-c(2:3)]


colnames(genecnv) <- substr(colnames(genecnv), 1, 12)
colnames(genecnv) <- gsub("[.]", "-", colnames(genecnv))
colnames(genecnv)[1] <- "GeneName"

expr <- expr[,which(colnames(expr) %in% colnames(genecnv))]
genecnv <- genecnv[,which(colnames(genecnv)  %in% colnames(expr))]
genecnv <- genecnv[,-grep("[.]1", colnames(genecnv))]
expr <- expr[, match(colnames(genecnv), colnames(expr))]
expr <- expr[rownames(expr)  %in% rownames(genecnv),]
expr <- expr[match(rownames(genecnv), rownames(expr)),]

expr.z <- as.data.frame(scale(expr[,2:132], center = rep(0,131), scale = apply(expr[,2:132], 2, sd, na.rm =  TRUE)))
expr.z$GeneName <- expr$GeneName
colnames(expr.z) <- gsub("[.]", "-", colnames(expr.z))
expr.z <- expr.z[,which(colnames(expr.z) %in% colnames(genecnv))]
genecnv <- genecnv[,which(colnames(genecnv)  %in% colnames(expr.z))]
#genecnv <- genecnv[,-grep("[.]1", colnames(genecnv))]
expr.z <- expr.z[, match(colnames(genecnv), colnames(expr.z))]
expr.z$median <- apply(expr.z[,2:132], 1, median)

expr.z$GeneName[ expr.z$median > 2 | expr.z$median < 2]

exprCNVBoxPlot <- function( expr, genecnv, GENENAME)
{
  k <- dim(expr)[2]
  expr.dat <- as.numeric(expr[ which(expr$GeneName == GENENAME),2:k])
  cnv.dat <- as.factor(genecnv[ which(genecnv$GeneName == GENENAME),2:k])
  
  #fit <- lm( expr.dat ~ cnv.dat)
  
  dat <- data.frame(expr = expr.dat, cnv = cnv.dat)
  
  p1 <- ggplot(data = dat, aes(x = cnv, y  = log2(expr), fill = cnv)) +
    geom_boxplot() + 
    theme_minimal() +
    ggtitle(GENENAME)
  return(p1)
}

png("DMRT3-expr.png")
exprCNVBoxPlot(expr, genecnv, "DMRT3")
dev.off()

exprCNVBoxPlot(expr, genecnv, "ATF7IP")
exprCNVBoxPlot(expr, genecnv, "TECRL")

png("NIPSNAP3A-expr.png")
exprCNVBoxPlot(expr, genecnv, "NIPSNAP3A")
dev.off()

png("NIPSNAP3B-expr.png")
exprCNVBoxPlot(expr, genecnv, "NIPSNAP3B")
dev.off()

exprCNVBoxPlot(expr, genecnv, "RAD52")

png("PPP2R3X-expr.png")
exprCNVBoxPlot(expr, genecnv, "PPP2R3C")
dev.off()

png("RAD51-expr.png")
exprCNVBoxPlot(expr, genecnv, "RAD51")
dev.off()


genecnv <- read.table("TGCT.focal_score_by_genes.txt", header = TRUE, sep = "\t")
genecnv$Gene.Symbol <- unlist(lapply(strsplit(genecnv$Gene.Symbol, "[.]"), function(x) x[1]))
genecnv$Gene.ID <- hsapiens_genes$hgnc_symbol[match(genecnv$Gene.Symbol, hsapiens_genes$ensembl_gene_id)]
colnames(genecnv)[4:159] <- gsub("[.]", "-", colnames(genecnv)[4:159])

map <- read.table("aliquot.tsv", sep = "\t", header  = TRUE)
map <- map[ -grep("10A", map$aliquot_submitter_id),]


map2 <- map[ map$aliquot_id %in% colnames(genecnv),]
ind <- match( colnames(genecnv), map2$aliquot_id)
colnames(genecnv)[!is.na(ind)] <- map2$aliquot_submitter_id[ind[!is.na(ind)]]

ind <- substr( map$aliquot_id, 1, 1 ) %in% seq(1:9) 
map$aliquot_id[ind] <- paste0("X", map$aliquot_id[ind])
map2 <- map[ map$aliquot_id %in% colnames(genecnv),]
ind <- match( colnames(genecnv), map2$aliquot_id)
colnames(genecnv)[!is.na(ind)] <- map2$aliquot_submitter_id[ind[!is.na(ind)]]

sample$Sample.ID <- substr(sample$Sample.ID, 1, 12)

# there is some kind of error on TCGA, many subjects have the wrong UUID on aliquot sheet
# i  had to get these manually
# colnames(genecnv)[which(colnames(genecnv) == "X9c74aae8-d90f-4130-83ca-1d97a63eb8a4")] <- "TCGA-2G-AAEW"
# colnames(genecnv)[which(colnames(genecnv) == "X7a62e1e2-191c-41ff-8215-dc6c74375032")] <- "TCGA-2G-AAGF"
# colnames(genecnv)[which(colnames(genecnv) == "X0d7960d7-76ad-4819-b664-fce1396e3c47")] <- "TCGA-2G-AAG6"
# colnames(genecnv)[which(colnames(genecnv) == "X934af87d-ce13-48e4-ac1e-6da66c656378")] <- "TCGA-2G-AAFO"
# colnames(genecnv)[which(colnames(genecnv) == "X5f863909-0821-448c-a025-4c1adf1f38c3")] <- "TCGA-2G-AAHC"
# colnames(genecnv)[which(colnames(genecnv) == "X0c8dd81d-45de-4d66-b361-efcbd2a2fe8c")] <- "TCGA-2G-AAL7"
# colnames(genecnv)[which(colnames(genecnv) == "X7887c5ce-dfed-4d95-b3d2-07299f61faa5")] <- "TCGA-2G-AAHP"   
# colnames(genecnv)[which(colnames(genecnv) == "X3f6c8076-0de7-4d47-90c2-f8df87018f2f")] <- "TCGA-2G-AAHP"
# colnames(genecnv)[which(colnames(genecnv) == "d4ac2166-ac65-4cff-aea9-6a49cd41cce5")] <- "TCGA-2G-AAF4"
# colnames(genecnv)[which(colnames(genecnv) == "X3803123b-39df-436e-8e84-20d0d6cbcaeb")] <- "TCGA-2G-AAKO"
# colnames(genecnv)[which(colnames(genecnv) == "e742adf4-f097-4e59-80a5-a79d6a4b9b62")] <- "TCGA-2G-AAFY"
# colnames(genecnv)[which(colnames(genecnv) == "be309c65-4455-4ab5-9414-cdb457887b15")] <- "TCGA-2G-AAFI"
# colnames(genecnv)[which(colnames(genecnv) == "f48b91dd-870e-46bb-aa4c-fd131a29e6e3")] <- "TCGA-2G-AAGK"
# colnames(genecnv)[which(colnames(genecnv) == "c42b5294-4208-4f87-973d-8be031331344")] <- "TCGA-2G-AALQ"
# colnames(genecnv)[which(colnames(genecnv) == "f2adc1f2-574d-4ca1-a410-1bf87aec3702")] <- "TCGA-4K-AA1H"
# colnames(genecnv)[which(colnames(genecnv) == "f76cf8d2-97f4-4b0d-a607-4127ce2c9528")] <- "TCGA-S6-A8JX"
# colnames(genecnv)[which(colnames(genecnv) == "a6e2b829-855e-4006-877a-c8d73c8c03a4")] <- "TCGA-W4-A7U2"
# colnames(genecnv)[which(colnames(genecnv) == "a8754cf8-bef5-4541-9b65-fcf9b6e20a74")] <- "TCGA-S6-A8JW"
# colnames(genecnv)[which(colnames(genecnv) == "cc85d446-9213-4c7b-8caf-3aca76187c04")] <- "TCGA-WZ-A8D5"
# colnames(genecnv)[which(colnames(genecnv) == "cb63e5cc-0a85-45b6-8f7c-b2be9fdb7e77")] <- "TCGA-ZM-AA0D"
# colnames(genecnv)[which(colnames(genecnv) == "de61c881-1a5a-4530-86cf-ab51742c43c0")] <- "TCGA-2G-AAGS"
# colnames(genecnv)[which(colnames(genecnv) == "X49a23ebc-c3e5-4b3a-8ce0-c746bc181b0c")] <- "TCGA-VF-A8AA"
# colnames(genecnv)[which(colnames(genecnv) == "f60272df-9864-47f8-84b4-cc567947a2a7")] <- "TCGA-XE-A8H4"
# colnames(genecnv)[which(colnames(genecnv) == "dfedc902-d551-42c9-9644-da0cd6304565")] <- "TCGA-2G-AAGI"
# colnames(genecnv)[which(colnames(genecnv) == "f1598631-a64f-41c1-b06c-d2235a21e2d9")] <- "TCGA-2G-AAH2"
# colnames(genecnv)[which(colnames(genecnv) == "X7887c5ce-dfed-4d95-b3d2-07299f61faa5")] <- "TCGA-2G-AAHP"
# colnames(genecnv)[which(colnames(genecnv) == "X7bd03897-b4fc-4502-8d16-4d15ad94250f")] <- "TCGA-2G-AAKO"
# colnames(genecnv)[which(colnames(genecnv) == "X13096bae-7366-4bfa-8d13-335867138e54")] <- "TCGA-2G-AAF1"
# colnames(genecnv)[which(colnames(genecnv) == "X79ad4b6e-afe1-4667-8bc7-467e0bd33794")] <- "TCGA-2G-AAGV"
# colnames(genecnv)[which(colnames(genecnv) == "X53cf153a-ce45-48c2-8175-ca58ba7d2a11")] <- "TCGA-2G-AALZ"
# colnames(genecnv)[which(colnames(genecnv) == "X198e3b27-4f5f-4b50-9973-b38b2b027111")] <- "TCGA-2G-AAKM"
# colnames(genecnv)[which(colnames(genecnv) == "b30b29d9-c821-48ae-b14b-74dbbeb6ee50")] <- "TCGA-2G-AALG"
# colnames(genecnv)[which(colnames(genecnv) == "X0b2605c1-3151-4d08-accd-4b6f020af31b")] <- "TCGA-2G-AAFE"
# colnames(genecnv)[which(colnames(genecnv) == "X8f4dda2b-c5da-4c09-a96c-6dfeb2b2fd81")] <- "TCGA-2G-AAGX"
# colnames(genecnv)[which(colnames(genecnv) == "eda339e2-6c8e-49fe-872a-472de3330363")] <- "TCGA-2G-AALX"
# colnames(genecnv)[which(colnames(genecnv) == "X6fa97436-7040-4cbe-a2d0-50cbc44a09ce")] <- "TCGA-XE-AANR"
# colnames(genecnv)[which(colnames(genecnv) == "X338ac80e-3a00-40ff-a02b-b5fbcf65a43b")] <- "TCGA-YU-A90Q"
# colnames(genecnv)[which(colnames(genecnv) == "X16a6fc56-815f-4f7e-957f-f3003237ae3a")] <- "TCGA-W4-A7U4"
# colnames(genecnv)[which(colnames(genecnv) == "X5bf8745f-80ba-40da-b1c4-6e316121cf49")] <- "TCGA-2G-AAH0"
# colnames(genecnv)[which(colnames(genecnv) == "X327a3956-e69b-4277-8c08-977bf09c0196")] <- "TCGA-2G-AAHG"
# colnames(genecnv)[which(colnames(genecnv) == "X6b1e6996-e154-4f7b-b40a-43f2705f9d4f")] <- "TCGA-2G-AAGI"
# colnames(genecnv)[which(colnames(genecnv) == "X388c6f11-f467-4978-8d47-07d5500ce3aa")] <- "TCGA-WZ-A7V4"
# colnames(genecnv)[which(colnames(genecnv) == "X13cd303a-1d78-4457-b664-4190074c2609")] <- "TCGA-2G-AALT"
# colnames(genecnv)[which(colnames(genecnv) == "de733352-5f67-474d-b5b7-5a944d2b1d32")] <- "TCGA-2G-AALN"
# colnames(genecnv)[which(colnames(genecnv) == "ecbdbf4d-8841-49ae-a48b-f073ceec62eb")] <- "TCGA-2G-AAGN"
# colnames(genecnv)[which(colnames(genecnv) == "aa07e0a4-624b-4470-a394-107ebf7be617")] <- "TCGA-2G-AALS"
# colnames(genecnv)[which(colnames(genecnv) == "X7264f8a4-344e-45f2-873e-e576c71765f6")] <- "TCGA-XE-AAOD"
# colnames(genecnv)[which(colnames(genecnv) == "X013b1afa-62fd-4516-b6cd-d1c713a5a8a7")] <- "TCGA-X3-A8G4"
# colnames(genecnv)[which(colnames(genecnv) == "X6174d04e-5460-4d37-9861-a203649df98f")] <- "TCGA-XE-A8H1"
# colnames(genecnv)[which(colnames(genecnv) == "X8c8fbd4a-8a34-4a30-9328-9320b2ffe339")] <- "TCGA-XE-AAO3"
# colnames(genecnv)[which(colnames(genecnv) == "X271977f5-9f89-4fc6-af74-a2b65814d282")] <- "TCGA-YU-A90S"
# colnames(genecnv)[which(colnames(genecnv) == "X886b97b1-c8f2-48dd-be8d-4683535302f7")] <- "TCGA-2X-A9D6"
# colnames(genecnv)[which(colnames(genecnv) == "X4646490c-961f-4a43-93ef-df8987ef114d")] <- "TCGA-W4-A7U3"
# colnames(genecnv)[which(colnames(genecnv) == "X710ae249-e496-430f-a36e-1eb23b3f0b68")] <- "TCGA-ZM-AA0F"
# colnames(genecnv)[which(colnames(genecnv) == "X64608cec-e4ab-4a47-baae-9a685e79448c")] <- "TCGA-XY-A9T9"
# colnames(genecnv)[which(colnames(genecnv) == "X0f2146a6-c393-48b3-b3cf-e445cbb30345")] <- "TCGA-ZM-AA0E"
# colnames(genecnv)[which(colnames(genecnv) == "X706ab709-7931-43d1-9e7a-2dd8da78ce90")] <- "TCGA-VF-A8A9"

#missing <- genecnv[ ,-(grep("TC", colnames(genecnv)))]

colnames(genecnv)[4:159] <- substr(colnames(genecnv)[4:159], 1, 12)
sum(sample$Sample.ID %in%  colnames(genecnv))




genecnv <- genecnv[,-which( colnames(genecnv) == "Cytoband")]







