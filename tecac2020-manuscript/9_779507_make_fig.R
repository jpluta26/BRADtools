#!/usr/bin/env Rscript
setwd("C:/Users/jpluta/Desktop/TECAC-fulldata/Manuscript/figures")
# pluta 3/21/18

# script to plot conditional analysis data


library(ggplot2)


# ----- setup variables -----
refsnp <- "9:779507"
ref.p  <- 1.04e-15
ref.chr <- as.character(strsplit(refsnp, ":")[[1]][1])
ref.bp <- as.integer(strsplit(refsnp, ":")[[1]][2])

# colon doesnt work well with unix, change to underscore
snpname <- paste(ref.chr, ref.bp, sep = "_")

print(paste("Begin plotting snp: ", snpname, sep = ""))

# these files are generated from conditionalAnalysis.sh
# snp-level data
COJOFILE <- paste(snpname, "cma.cojo", sep = ".") 

# file with the names of snps in LD calculation
SNPFILE <- paste(snpname, "snp.ld", sep = ".")

# r value corresponding to each snp in SNPFILE
LDFILE <- paste(snpname, "r.ld", sep = ".")
# ------




# ----- read in and format the data ----- #
# the conditional analysis file
dat <- read.table(COJOFILE, header = TRUE, 
                  col.names = c("Chr", "SNP", "bp", "refA", "freq", "b", 
                                "se", "p", "n", "freq_geno", "bC", "bC_se", "pC"), 
                  colClasses = c("character", "character", "integer", "character", 
                                 "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric"))

# ------



# ----- calculate r2 for each snp and attach to the main df -----
# ld information is stored separately from the corresponding snps
# bind them together and convert to r^2
ld <- t(read.table(LDFILE, header=F))
snp <- t(read.table(SNPFILE, header=F))
ld.mat <- data.frame(snp <- snp[,1], r2 <- (ld[,1])^2)
colnames(ld.mat) <- c("snp", "r2")

# put r2 values in the main dataframe
dat$r2 <- 0
dat$r2 <- ld.mat$r2[match(dat$SNP, ld.mat$snp)]

# error checking, can probably drop this
if( any(is.na(dat$r2)) )
{
  dat$r2[is.na(dat$r2)] <- -9
  print("NA detected in ld file.")
}
# ------



# ----- flag each snp as dependent or independent of the reference snp -----
# a snp is dependent (dat$indp = FALSE) if the p-value of the tested snp
# goes from nominal significance (p < 0.00001) to below significance (p > 0.00001)
# after modeling with the reference snp 
# *** might need to tune this; criteria for determining independence
dat$indp <- !(dat$p < 0.00001 & dat$pC > 0.00001)

#if( any( !dat$indp) )
#{
#	stop("no indepednent snps!")
#}

# snps with r2 > 0.8 are flagged as dependent
dat$indp[dat$r2 > 0.8] <- FALSE
dat$indp[dat$r2 <= 0.4] <- TRUE
# snps that are too colinear with the reference are not run
# flag these as dependent
if( any( is.na(dat$pC) ) )
{
  dat$indp[which(is.na(dat$pC))] <- FALSE
}



write.table(dat[ dat$indp == FALSE & dat$r2 >= 0.3, ], paste(snpname, "_credSet.txt", sep = ""), col.names = T, row.names = F, quote = F, append = F)



# identify the top dependent and independent snps
# separate them from the data so they can be plotted as a 
# seperate layer



# top independent snp
ind <- which(dat$bp == 878563)
topindsnp.bp <- dat$bp[ind]
topindsnp.p <- -log10(dat$p[ind])
topindsnp.snp <- dat$SNP[ind]
dat <- dat[-ind,]
# -----

ind <- which(dat$bp == 834406)
snp2.bp <- dat$bp[ind]
snp2.p <- -log10(dat$p[ind])
snp2.snp <- dat$SNP[ind]
dat <- dat[-ind,]

ind <- which(dat$bp == 863635)
snp3.bp <- dat$bp[ind]
snp3.p <- -log10(dat$p[ind])
snp3.snp <- dat$SNP[ind]
dat <- dat[-ind,]

# ----- plot the data -----
# probably a smarter way to do this in ggplot
# n is the number of levels in the gradient color part for 
# plotting ld
n <- 20
dat$ld.f <- 0
for( i in 1:n)
{
  dat$ld.f[which(dat$r2 > ((i - 1)/10) & dat$r2 < (i/10) )] <- i/10
}

# calculate offset for text annotation
xmin <- min(c(dat$bp, topindsnp.bp, topsnp.bp, ref.bp))
xmax <- max(c(dat$bp, topindsnp.bp, topsnp.bp, ref.bp))
offset <- (xmax - xmin) / 20

# y limit is the lowest p-value
ymax <- max(c(-log10(dat$p), -log10(ref.p), topindsnp.p, topsnp.p))

# make the big ggplot
# ref snp is plotted separately as an inverted triangle
# top indp snp and top dep snp are plotted separately, in pink, with corresponding shape
p1 <- ggplot(data = dat, aes(x=bp, y=-log10(p), shape=indp, colour=dat$ld.f)) + 
  scale_color_gradient2(low="red", high="blue", mid="green", midpoint=0.5, name="LD", limits=c(0,1)) + 
  theme_set(theme_light(base_size = 28)) + 
  ylim(0, 75) +
  theme_minimal() +
  theme(plot.title = element_text(size = 28, face = "bold"), 
        legend.text = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 22, face = "bold"),
        axis.text.x = element_text(size = 16), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20), 
        axis.text.y = element_text(size = 16)) + 
  geom_hline(yintercept = -log10(1e-5), linetype="dashed", colour="blue") +
  geom_point(size=4) + 
  geom_point(x=ref.bp, y=-log10(ref.p), colour = "black", fill="deeppink", shape=21, size=4) +
  geom_point(x = topindsnp.bp, y = topindsnp.p, colour = "black", fill = "deeppink", shape = 24, size=4) +
  geom_point(x = snp2.bp, y = snp2.p, colour = "black", fill = "deeppink", shape = 24, size=4) +
  geom_point(x = snp3.bp, y = snp3.p, colour = "black", fill = "deeppink", shape = 24, size=4) +
 # geom_point(x = topsnp.bp, y = topsnp.p, colour = "black", fill = "deeppink", shape = 21, size = 4) +
#  annotate("text", x = round(topsnp.bp + offset), y = topsnp.p, hjust = -0.5, label = topsnp.snp, size = 8) +
  annotate("text", x = round(ref.bp  - 135000), y = -log10(ref.p), label = refsnp, size = 8) +
  annotate("text", x = topindsnp.bp + 155000, y = topindsnp.p, label = topindsnp.snp, size = 8) +
  annotate("text", x = snp2.bp + 155000, y = snp2.p, label = snp2.snp, size = 8) +
  annotate("text", x = snp3.bp + 155000, y = snp3.p, label = snp3.snp, size = 8) +
  #ggtitle(paste(refsnp, "Conditional Plot", sep = " ")) + 
  scale_shape_discrete(name="Independence", breaks=c("FALSE", "TRUE"), 
                       labels=c("FALSE", "TRUE"), guide=guide_legend(reverse=TRUE)) 
# -----

p1

# ----- write output files ----- #
# get the time stamp for output names
dt <- strsplit(date(), " ")[[1]][c(1,3,5)]
dt <- paste(dt[1], dt[2], dt[3], sep="_")

print(paste("Top independent snp: ", topindsnp.snp, sep=""))
print(paste("Top dependent snp: ", topsnp.snp, sep=""))
print("writing image...")
png(paste(snpname,  "cond.png", sep="_"), width=700, height=900)
print(p1)
dev.off()
print("done!")

OUTNAME <- paste(snpname, "summary.txt", sep = "_")
write.table(dat, OUTNAME, append = FALSE, quote = FALSE, col.names = TRUE, row.names = FALSE)
print("successfully completed!!")
# -----