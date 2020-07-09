
# pluta 9/18/18

library(data.table)
library(qqman)
setwd("~/Desktop/tecac-manuscript/")
# script to create manhattan plot of either a single chromosome of the whole genome
# this script greatly increases the speed of plotting by subsampling a small fraction
# of the non-significant snps

# INPUT: INFILE (text file) with two columns, MarkerName (chr:pos) or rsid and 
# P-value
#
# OUTPUT: OUTFILE (string), the name of the pdf to output

#setwd("/Users/jpluta/Desktop/TECAC-fulldata/Manuscript/figures")

# the manhattan plot function from qqman, but modified to add text and labeled top hits
manhattan1<-function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10",
                                                                               "gray60"), 
                      chrlabs = NULL, 
                      suggestiveline = -log10(1e-05),
                      genomewideline = -log10(5e-08), 
                      highlight1 = NULL, 
                      highlight2 = NULL, 
                      highlight3 = NULL, 
                      logp = TRUE,
                      
                      ...)
{
  if (!(chr %in% names(x)))
    stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x)))
    stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x)))
    stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x)))
    warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]]))
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]]))
    stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]]))
    stop(paste(p, "column should be numeric."))
  d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
  if (!is.null(x[[snp]]))
    d = transform(d, SNP = x[[snp]])
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  } else 
    {
      d$logp <- d$P
    }
  d$pos = NA
  d$index = NA
  ind = 0
  for (i in unique(d$CHR)) {
    ind = ind + 1
    d[d$CHR == i, ]$index = ind
  }
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    options(scipen = 999)
    d$pos = d$BP/1e+06
    ticks = floor(length(d$pos))/2 + 1
    xlabel = paste("Chromosome", unique(d$CHR), "position(Mb)")
    labs = ticks
  } else 
    {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      } else 
      {
        lastbase = lastbase + tail(subset(d, index ==
                                            i - 1)$BP, 1)
        d[d$index == i, ]$pos = d[d$index == i, ]$BP +
          lastbase
      }
      ticks = c(ticks, (min(d[d$CHR == i, ]$pos) + max(d[d$CHR ==
                                                           i, ]$pos))/2 + 1)
    }
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i",
                   las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0,
                                                                     ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log10))
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in%
                                            names(dotargs)]))
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  if (nchr == 1) {
    axis(1, ...)
  }
  else {
    axis(1, at = ticks, labels = labs, ...)
  }
  col = rep(col, max(d$CHR))
  if (nchr == 1) {
    with(d, points(pos, logp, pch = 20, col = col[1], ...))
  }
  else {
    icol = 1
    for (i in unique(d$index)) {
      with(d[d$index == unique(d$index)[i], ], points(pos,
                                                      logp, col = col[icol], pch = 20, ...))
      icol = icol + 1
    }
  }
  if (suggestiveline)
    abline(h = suggestiveline, col = "blue")
  if (genomewideline)
    abline(h = genomewideline, col = "red")
  if (!is.null(highlight1)) {
    if (any(!(highlight1 %in% d$SNP)))
      warning("You're trying to highlight1 SNPs that don't exist in your results.")
    d.highlight1 = d[which(d$SNP %in% highlight1), ]
    with(d.highlight1, points(pos, logp, bg = "green2", col = "black", pch = 21, cex = 1.1,
                              ...))
  }
  if (!is.null(highlight2)) {
    if (any(!(highlight2 %in% d$SNP)))
      warning("You're trying to highlight2 SNPs that don't exist in your results.")
    d.highlight2 = d[which(d$SNP %in% highlight2), ]
    with(d.highlight2, points(pos, logp, bg = "dodgerblue", col = "black", pch = 22, cex = 1.1,
                              ...))
  }
  if (!is.null(highlight3)) {
    if (any(!(highlight3 %in% d$SNP)))
      warning("You're trying to highlight2 SNPs that don't exist in your results.")
    d.highlight3 = d[which(d$SNP %in% highlight3), ]
    with(d.highlight3, points(pos, logp, bg = "red", col = "black", pch = 23, cex = 1.1,
                              ...))
  }
  
 # add highlighting of new/prior hits
  snps <- read.table("tecac-snps.csv", header = TRUE, sep = ",")
  snps$lab <- ''
  snps$lab[ snps$type == "Novel" ] <- letters[ seq( from = 1, to = 22)]
  snps$lab[ snps$type == "Replicated" ] <- seq(from = 1, to = 45)
  snps$lab[ nchar(snps$lab) == 2] <- paste0(" ", snps$lab[ nchar(snps$lab) == 2])
  
 #par(xpd = TRUE)
 # tmp <- subset(d, d$SNP %in% snps$snp)
#  snps <- snps[ snps$snp %in% tmp$SNP,]
  
  
  # add hits labeling
 # ind <- match(tmp$SNP, snps$snp )
#  snps <- snps[ match(tmp$SNP, snps$snp ), ]
#  tmp$lab <- snps$lab

  #                                    
#  tmp$pos[ tmp$SNP == "4:76550768"] <- 686864431
#  tmp$pos[ tmp$SNP == "19:22228091"] <- 2689971882
  #segments(2679971882, 10.069917, 2689971882, 11)
  #print(paste0("tmp$pos = ", tmp$pos))
  #text(tmp$pos + 35000000, tmp$logp, label = as.character(tmp$lab), cex = 1.5)
  #print(tmp)

}

INFILE  = "meta_newVal_8site_ggman.txt"
OUTFILE = "meta_newVal_8site"


# PARAMS
# the threshold for p-values; sub-sample any values with p > p.thresh
p.thresh = 0.001

# proportion of snps to keep
k.thresh = 0.1


MANOUTFILE <- paste(OUTFILE, "_man.png", sep="")
QQOUTFILE <- paste(OUTFILE, "_qq.png", sep = "")
# ---

print("reading association testing data...")

dat <- fread(INFILE, header=TRUE)
dat <- as.data.frame(dat)

# check input validity
# important for consistency with locuszoom
if( colnames(dat)[1] != "MarkerName")
{
	stop("Colnames of INFILE should be 'MarkerName' and 'P-value' ")
}


dat$MarkerName <- as.character(dat$MarkerName)
dat$CHR <- unlist(lapply(strsplit(dat$MarkerName, ":"), function(x) x[1]))
dat$BP <- unlist(lapply(strsplit(dat$MarkerName, ":"), function(x) x[2]))
dat$BP <- as.numeric(dat$BP)


if( any(dat$CHR == "X"))
{
        dat$CHR[which(dat$CHR == "X")] <- 23
}

dat$CHR <- as.numeric(dat$CHR)


colnames(dat) <- c("rsid", "p", "chr", "bp")
dat <- dat[!is.na(dat$chr),]
chrlabs <- as.character(seq(1:length(unique(dat$CHR))))
dat$p <- as.numeric(dat$p)
dat <- dat[ !is.na(dat$chr), ]
print("done!")

snps <- read.table("tecac-snps.csv", header = TRUE, sep = ",")

# original tecac-snps.csv was out of date, these were some corrections
#snps$snp <- as.character(snps$snp)
#snps$type[ snps$snp == "14:55912663"] <- "Not Replicated"
#snps$snp[ snps$snp == "6:33533625"] <- "6:33542478"
#write.table(snps, "tecac-snps.csv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ",")

# 44 replicated, 12 not

snps$type <- factor(snps$type, levels = c("Novel", "Replicated", "Not Replicated"))
p1 <- ggplot(data = snps, aes(x = seq(1:78), y= seq(1:78), fill= type,shape = type)) +
  geom_point(size = 3.5) +
  scale_color_manual(values = c(
    "Novel" = "dodgerblue",
    "Replicated" = "green2",
    "Not Replicated" = "red"), aesthetics = c("colour", "fill")) +
  scale_shape_manual(values = c(22,21,23))
png("legend.png")
p1
dev.off()
# ---

print("subsampling association data...")

# s is the sub sample
x <- dim(dat[which(dat$p > p.thresh),])[1]
s <- sample(dat$rsid[which(dat$p > p.thresh)], round(k.thresh * x), replace=F)

print("done!")


# ---

print("creating manhattan plot...")

snps$lab <- ''
snps$lab[ snps$type == "Novel" ] <- letters[ seq( from = 1, to = 22)]


#dat <- dat[ dat$chr == 5 & dat$bp > 1280000 ,]

#dat2 <- dat[ dat$rsid %in% snps$snp,]
#ind <- match(dat2$rsid, snps$snp)
#dat2$lab <- snps$lab[ind]

png(MANOUTFILE, height=600, width=1200)
manhattan1(dat[(dat$rsid %in% dat$rsid[match(s, dat$rsid)]) | (dat$p <= p.thresh),],
        chr="chr", bp="bp", p="p", snp="rsid", 
	col=c("gray", "black"), #ylim = c(0,40),
        logp=TRUE, chrlabs = chrlabs, 
        highlight1 = snps$snp[snps$type == "Replicated"], 
        highlight2 = snps$snp[snps$type == "Novel"], 
        highlight3 = snps$snp[snps$type == "Not Replicated"]) 
        

dev.off()

#png(QQOUTFILE, height=600, width=600)
#qq(dat$p)
#dev.off()

print("done!")

p1 <- ggplot(data = snps, aes(x = bp, y = chr, fill = type, color = type)) + geom_point() +
   scale_color_manual(values = c("Replicated" = "green2", "Novel" = "dodgerblue", "Not Replicated" = "red"))
