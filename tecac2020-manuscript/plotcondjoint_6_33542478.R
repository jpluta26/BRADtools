
#!/usr/bin/env Rscript

# pluta 3/21/18

# conditional analysis conditioning on > 1 snp
setwd("C:/Users/jpluta/Desktop/TECAC-fulldata/Manuscript/figures")



library(ggplot2)

refsnpfile = "snps-pvals"

ref <- read.table(refsnpfile, header=F, col.names = c("snp", "p"), colClasses = c("character", "numeric"))

ref$bp <- as.integer(unlist(lapply(strsplit(ref$snp, ":"), function(x) x[2])))
ref.chr <- as.integer(strsplit(ref$snp[1], ":")[[1]][1]) # this is a constant



# these files are generated from conditionalAnalysis.sh
COJOFILE <- "6_33542478_2_joint.cma.cojo"

# not sure how to incorporate LD just yet
#SNPFILE <- paste(snpname, "snp.ld", sep = ".")
#LDFILE <- paste(snpname, "r.ld", sep = ".")


dat <- read.table(COJOFILE, header=T, 
                  col.names = c("Chr", "SNP", "bp", "refA", "freq", "b", 
                                "se", "p", "n", "freq_geno", "bC", "bC_se", "pC"), 
                  colClasses = c("character", "character", "integer", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                  "numeric", "numeric", "numeric"))



#ld <- t(read.table(LDFILE, header=F))
#snp <- t(read.table(SNPFILE, header=F))

# combine ld to each snp
#ld.mat <- data.frame(snp <- snp[,1], r2 <- (ld[,1])^2)
#colnames(ld.mat) <- c("snp", "r2")

# put r2 values in the main dataframe
#dat$r2 <- 0
#dat$r2 <- ld.mat$r2[match(dat$SNP, ld.mat$snp)]
#dat$r2[is.na(dat$r2)] <- 0

# might need to tune this
dat$indp <- !(dat$p < 0.00001 & dat$pC > 0.00001)

if( any( is.na(dat$pC) ) )
{
	dat$indp[which(is.na(dat$pC))] <- FALSE
}

# top dependent snp
ind <- which(dat$p == min(dat$p[!dat$indp]))
topsnp.bp <- dat$bp[ind]
topsnp.p <- -log10(dat$p[ind])
topsnp.snp <- dat$SNP[ind]

# top independent snp
ind <- which(dat$p == min(dat$p[dat$indp]))
topindsnp.bp <- dat$bp[ind]
topindsnp.p <- -log10(dat$p[ind])
topindsnp.snp <- dat$SNP[ind]


# probably a smarter way to do this in ggplot
# n is the number of levels in the gradient
#n <- 20
#dat$ld.f <- 0
#for( i in 1:n)
#{
#  dat$ld.f[which(dat$r2 > ((i - 1)/10) & dat$r2 < (i/10) )] <- i/10
#}


# make the big ggplot
#p1 <- ggplot(data = dat, aes(x=bp, y=-log10(p), shape=indp, colour=dat$ld.f)) + 
#  scale_color_gradient2(low="red", high="blue", mid="green", midpoint=0.5, name="LD") + 
#  theme_light() + 
#  theme_linedraw() +
#  geom_point() + geom_point(x=ref.bp, y=-log10(ref.p), colour = "black", fill="blue", shape=24, size=4) +
#  annotate("text", x = topsnp.bp + 51000, y = topsnp.p, label = topsnp.snp) +
#  annotate("text", x = topindsnp.bp + 51000, y = topindsnp.p, label = topindsnp.snp) +
#  ggtitle(paste(refsnp, "Conditional Plot", sep = " ")) + 
#  scale_shape_discrete(name="Independence", breaks=c("FALSE", "TRUE"), labels=c("FALSE", "TRUE"), guide=guide_legend(reverse=TRUE))


xmin <- min(c(dat$bp, topindsnp.bp, topsnp.bp, ref$bp))
xmax <- max(c(dat$bp, topindsnp.bp, topsnp.bp, ref$bp))
offset <- (xmax - xmin) / 20

ymax <- max(c(-log10(dat$p), -log10(ref$p), topindsnp.p, topsnp.p))

p1 <- ggplot( data = dat, aes(x=bp, y=-log10(p), colour=indp)) +
  labs(colour = "Independent") +
  scale_colour_manual(values= c("red", "blue")) +
      theme_minimal() +
      theme(plot.title = element_text(size = 18, face = "bold"),
	    legend.text = element_text(size = 14, face = "bold"),
            legend.title = element_text(size = 18, face = "bold"), 
            axis.text.x = element_text(size = 10),
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            axis.text.y = element_text(size = 12)) +
      ylim(0, ymax + 0.5) +
      geom_point(size=3) +
      geom_hline(yintercept = -log10(1e-5), linetype="dashed", colour="blue") +
 #     geom_point(x = topsnp.bp, y = topsnp.p, size = 4, shape = 21, colour = "black", fill = "red") +
#      annotate("text", x = round(topsnp.bp - 4 * offset), y = topsnp.p, label = topsnp.snp, size = 4, hjust = -0.01) +
#      geom_point(x = topindsnp.bp + offset, y = topindsnp.p, size = 4, shape = 21, colour = "black", fill = "blue") +
 #     annotate("text", x = round(topindsnp.bp + offset), y = topindsnp.p, label = topindsnp.snp, size = 4, hjust = -0.15) +
      ggtitle(paste(ref$snp[1], "Joint Conditional Plot", sep = " ")) +
        geom_point(x = ref$bp[1], y = -log10(ref$p[1]), colour="black", fill="green", shape = 25, size = 4) +
        annotate("text", x = ref$bp[1] + offset * .6, y = -log10(ref$p[1]), label = ref$snp[1], size= 4, hjust = -0.05) +
        geom_point(x = ref$bp[2], y = -log10(ref$p[2]), colour="black", fill="green", shape = 25, size = 4) +
        annotate("text", x = ref$bp[2] - 120000, y = -log10(ref$p[2]), label = ref$snp[2], size= 4, hjust = -0.05)
      

print(p1)
png("joint-633542478e.png", height = 600, width = 550)
print(p1)
dev.off()
# format date for output
x <- c(1,3,5)
dt <- strsplit(date(), " ")[[1]][x]
dt <- paste(dt[1], dt[2], dt[3], sep="_")

# ensure a unique output name
OUTNAME=paste(ref.chr, ref$bp[1], dt, "jt_cond.png", sep="_")

if(file.exists(OUTNAME))
{
	i=2
	repeat 
	{
		OUTNAME = paste(ref.chr, ref$bp[1], dt, "jt_cond", i, ".png", sep="_")
		if( !file.exists(OUTNAME) )
		{
			break
		} else
		  { 
			i = i + 1
		  }
	}	
}

print(paste("Top independent snp: ", topindsnp.snp, sep=""))
print(paste("Top dependent snp: ", topsnp.snp, sep=""))
print("writing image...")


#png(OUTNAME, width=1000, height=1200)
print(p1)
#dev.off()


print("done!")
print("successfully completed")
