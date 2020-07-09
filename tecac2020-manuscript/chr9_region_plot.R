setwd("C:/Users/jpluta/Desktop/TECAC-fulldata/Manuscript/figures/")

library(data.table)
library(ggplot2)

getSNPs.in.LD <- function( snpname, dat )
{
  snp.fname <- paste( snpname, "snp.ld", sep = ".")
  ld.fname <- paste(snpname, "r.ld", sep = ".")
  
  snp <- t(read.table(snp.fname, header = F))
  ld <- t(read.table(ld.fname, header = F))
  ld.mat <- data.frame(snp <- snp[,1], r2 <- (ld[,1])^2)
  colnames(ld.mat) <- c("snp", "r2")
  ld.mat$snp <- as.character(ld.mat$snp)
  ld.mat <- ld.mat[ ld.mat$r2 >= .3,]
 
  ld.mat <- ld.mat[ld.mat$snp %in% dat$SNP,]
  ind <- match(ld.mat$snp, dat$SNP)
  
  dat$r2[ind] <- ld.mat$r2
  dat$grp[ind] <- snpname
  dat$grp[which(dat$SNP == sub("_", ":", snpname))] <- snpname
  return(dat)
}



dat <- fread("9_863635.cma.cojo", header = T)
dat <- as.data.frame(dat)

dat$grp <- "x"
dat$r2 <- 0
dat <- dat[ dat$bp > 7e05 & dat$bp < 9.5e05,]

dat <- getSNPs.in.LD("9_863635", dat)
dat <- getSNPs.in.LD("9_834406", dat)
dat <- getSNPs.in.LD("9_779507", dat)
dat <- getSNPs.in.LD("9_878563", dat)
#dat <- getSNPs.in.LD("9_815079", dat)

p1 <- ggplot(data = dat, aes(x = bp, y = -log10(p), colour = dat$grp)) + 
  geom_point(alpha = 0.5) +
  theme_minimal() + 
  scale_colour_manual( values = c("red", "blue", "green", "goldenrod", "black", "lightblue"))
p1
png("DMRT-region-plot.png")
print(p1)
dev.off()

p2 <- ggplot(data = dat, aes(x = bp, y = -log10(p), colour = "blue")) +
       geom_point() +
       geom_point(data = dat, aes(x = bp, y = -log10(pC), colour = "red"))
p2