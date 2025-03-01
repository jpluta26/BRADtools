setwd("~/Documents/nathansonlab/CNV")
library(ggplot2)

snp.dat  <-  read.table("IKN_TECAC_reclusterx14520_dbSNP.bim", header = FALSE)
snp.dat <- snp.dat[  snp.dat$V1 > 0 & snp.dat$V1 < 24,]

p1 <- ggplot( ) + 
  geom_segment(data = snp.dat, aes(x = V4, xend = V4 + .001, y = V1,  yend = V1 + .5)) +
  theme_minimal() +
  xlab("Position") +
  ggtitle("Illumina Infinium Core 24 Snp distribution - Raw") +
  scale_y_reverse("Chromosome", breaks = seq(1:23), labels = c(seq(1:22), "X"))

png("infinium24-snpdistribution.png", width =  600, height = 400)
print(p1)
dev.off()

snp.dat  <-  read.table("case_ctl_cauc.bim", header = FALSE)
snp.dat <- snp.dat[  snp.dat$V1 > 0 & snp.dat$V1 < 24,]

p1 <- ggplot( ) + 
  geom_segment(data = snp.dat, aes(x = V4, xend = V4 + .001, y = V1,  yend = V1 + .5)) +
  theme_minimal() +
  xlab("Position") +
  ggtitle("Illumina Infinium Core 24 Snp distribution - QC") +
  scale_y_reverse("Chromosome", breaks = seq(1:23), labels = c(seq(1:22), "X"))

png("infinium24-snpdistribution-qc.png", width =  600, height = 400)
print(p1)
dev.off()
