# pluta 7/9/20
# script to generate a plot with 6 reference snps, and unique color/gradient for each
# used to show LD structure of multiple independent SNPs in a region
# ultimately this was not used for the paper; here for reference

library(ggplot2)
library(gridExtra)
library(ggnewscale)

# ---- function to assign a group by lead snp ----- #
getSNPgrp <- function( snp, dat )
{
  in.r <- paste0("19_", snp, ".r.ld")
  ld.dat <- t(read.table(in.r, header = FALSE))
  
  
  in.snp <- paste0("19_", snp, ".snp.ld")
  snp.dat <- t(read.table(in.snp, header = FALSE))
  
  ld.mat <- data.frame(snpname <- snp.dat[,1], r2 <- (ld.dat[,1]^2))
  colnames(ld.mat) <- c("snpname", "r2")
  
  dat$r2 <- 0
  dat$r2 <- ld.mat$r2[ match(dat$MarkerName, ld.mat$snpname)]
  
  if( any(is.na(dat$r2)))
  {
    dat$r2[is.na(dat$r2)] <- 0
    print("NA detected in ld file")
  }
  
  dat$grp[ dat$r2 > 0.4 ] <- paste("19", as.character(snp), sep = ":")
  return(dat)
}
# ----------------------------------------------- #





# ----------------------------------------------- #
# include all snps within ld of r2.thresh
r2.thresh <- 0.2

# snps of interest (bp)
# ref snps in chr19
snps <- c(  28356614,
            22228091,
            23203496,
            24050828,
            24149545,
            28257393  )
snps <- snps[order(snps)]

# snp summary stats
dat <- read.table("newVal_8site_chr19_ggman.txt", header = TRUE)

dat$MarkerName <- as.character(dat$MarkerName)
dat$CHR <- unlist(lapply(strsplit(dat$MarkerName, ":"), function(x) x[1]))
dat$BP <- unlist(lapply(strsplit(dat$MarkerName, ":"), function(x) x[2]))
dat$BP <- as.numeric(dat$BP)

# window for the plot
dat <- dat[ dat$BP >= (snps[1] - 50000) & dat$BP <= snps[6] + 50000,]
dat$grp <- "none"

# assign each snp to a group defined by its refsnp
for( i in 1:length(snps))
{
  dat <- getSNPgrp( snps[i], dat)
}

hits <- paste("19", snps, sep = ":")
dat$grp[ match(hits, dat$MarkerName) ] <- hits

dat$grp <- as.factor(dat$grp)
dat$type <- "Replicated"
dat$type[ dat$grp == "19:28356614"] <- "Novel"
dat$type <- as.factor(dat$type)

# set up an individual data frame for each snp
d1 = getSNPgrp(22228091, dat)
d1 = d1[ d1$grp == hits[1],]

d2  = getSNPgrp(23203496, dat)
d2 = d2[ d2$grp == hits[2],]

d3 = getSNPgrp(24050828, dat)
d3 = d3[ d3$grp == hits[3],]

d4 = getSNPgrp(24149545, dat)
d4 = d4[ d4$grp  == hits[4],]

d5 = getSNPgrp(28257393, dat)
d5 = d5[ d5$grp  == hits[5],]

d6 = getSNPgrp(28356614, dat)
d6 = d6[ d6$grp  == hits[6],]




# r plot with individual color and gradient for each snp
p1 <- ggplot( data = dat, aes(x = BP, y = -log10(P.value))) +
  geom_point(alpha = 0.1, color = "black") +
  geom_point( data = d1, aes(x = BP, y = -log10(P.value), color = r2)) +
  scale_fill_gradient(low = "white", high = "blue", aesthetics= "color") +
  new_scale_color() +
  
  geom_point( data = d2, aes(x = BP, y = -log10(P.value), color = r2)) +
  scale_fill_gradient(low = "white", high = "maroon", aesthetics= "color") +
  new_scale_color() +
 
  geom_point( data = d3, aes(x = BP, y = -log10(P.value), color = r2)) +
  scale_fill_gradient(low = "white", high = "green", aesthetics= "color") +
  new_scale_color() +
  
  geom_point( data = d4, aes(x = BP, y = -log10(P.value), color = r2)) +
  scale_fill_gradient(low = "white", high = "darkorchid1", aesthetics= "color") +
  new_scale_color() +
  
  geom_point( data = d5, aes(x = BP, y = -log10(P.value), color = r2)) +
  scale_fill_gradient(low = "white", high = "dodgerblue", aesthetics= "color") +
  new_scale_color() +

  geom_point( data = d6, aes(x = BP, y = -log10(P.value), color = r2)) +
  scale_fill_gradient(low = "white", high = "red", aesthetics= "color") +
  
  xlim(xmin, xmax) +
  theme_minimal()  +
  geom_point( data = dat[ dat$MarkerName == hits[1],], aes(x = BP, y = -log10(P.value)),
              color = "gray", fill = "blue", pch=23, size = 3) +
  geom_point( data = dat[ dat$MarkerName == hits[2],], aes(x = BP, y = -log10(P.value)),
              color = "gray", fill = "maroon", pch=23, size = 3) +
  geom_point( data = dat[ dat$MarkerName == hits[3],], aes(x = BP, y = -log10(P.value)),
              color = "gray", fill = "green", pch=23, size = 3) +
  geom_point( data = dat[ dat$MarkerName == hits[4],], aes(x = BP, y = -log10(P.value)),
              color = "gray", fill = "darkorchid1", pch=23, size = 3) +
  geom_point( data = dat[ dat$MarkerName == hits[5],], aes(x = BP, y = -log10(P.value)),
              color = "gray", fill = "dodgerblue", pch=23, size = 3) +
  geom_point( data = dat[ dat$MarkerName == hits[6],], aes(x = BP, y = -log10(P.value)),
              color = "gray", fill = "red", pch=24, size = 3) +
  theme(legend.position = "none") 

p1
# ----------------------------------------------- #
