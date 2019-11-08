# pluta 8/11/16

# ------------------------------ preprocessor ------------------------------ #
library(BMS)
library(ggplot2)

setwd("C:/Users/jpluta/Desktop/bayes")
# -------------------------------------------------------------------------- #



# ----------------- create odds ratio distribution plot -------------------- #
RRplot <- function( dat , PNGname )
  # input: dat, a dataframe with two entries: $RR, the odds ratios
  # and $MAF, the minRR allele frequencies. this is taken from SNP test
  #   PNGname, name of the output file
  # output: 
  # a plot of the distribution of Odds ratios
  # returns odds ratio threshold, below which 95% of RRs are contained
{
  # convert to density and get quantiles; want 95% coverage
  d1 <- density(dat$RR)
  RR.q <- quantile(d1, probs=c(0.95), names=FALSE)
  
  # convert density to data.frame fRR ggplot
  pdf.dat <- data.frame(x=d1$x, y=d1$y)
  
  p <- ggplot() + geom_line(data=pdf.dat, aes(x=x, y=y)) + 
    ggtitle("Distribution of Odds Ratios - CHR1") +
    xlab("Odds Ratio") + 
    ylab("Probability") + 
    geom_vline(xintercept = RR.q, linetype="dashed") + 
    annotate("text", x=RR.q + 0.25, y=max(pdf.dat$y) * .8, label=round(RR.q, 2)) +
    geom_ribbon(data=pdf.dat[pdf.dat$x < RR.q,], aes(x, ymin=0, ymax=y), fill="red", alpha=0.1)
  
  png(PNGname)
  print(p)
  dev.off()
  
  return(RR.q)
}
# -------------------------------------------------------------------------- # 


# -------------------------------------------------------------------------- #
getWBF <- function(V, W, B)
  # input: B, the coefficient estimate on the allele from logistic regression
  # ie log(RR)
  #    V = var(B)
  #    W = variance of the prior on beta.hat
{
  
  #  H0: B = 0;  B ~ N(0,V)
  #  Ha: B != 0; B ~ N(0, V + W)
  # WBF is a reduced form of the ratio of these pdfs
  z.sq <- (B^2) / V
  WBF <- sqrt( (V + W) / V ) * exp( -.5 * z.sq * (W / (V+W)))
  
  # note that this can be written equivalently as:
  # pH0 <- dnorm(B, m=0, s=sqrt(V))
  # pH1 <- dnorm(B, m=0, s=sqrt(V+W))
  # WBF <- pH0/pH1;
  
  
  return(WBF)
}
# -------------------------------------------------------------------------- #


# ------------------------------- getWBFfromPval ---------------------------- #
# convert p-value directly to bayes factRR
getWBFfromPval <- function(V, B.hat, K)
{
  z.sq <- (B.hat^2) / V
  WBF <- sqrt(1 + K) * exp( -.5 * z.q * (K/(1+K)))
  return(WBF)
}
# -------------------------------------------------------------------------- #

# ----------------------------- getW.MAFind --------------------------------- # 
# Effect-MAF Independence
# assumption: the variance, W, is independent of MAF
# RRu is the upper value of relative risk, above which we believe
# that relative risks will occur with low probability. The prior probability
# of a relative risk above RRu is q
getW.MAFind <- function( dat, q)
  # input: dat, the SNP data in data, q basically the alpha level
  # output: W
{
  RRu <- RRplot(dat, "allSNP.png")
  print(paste("RRu: ", RRu))
  # W increases with RRu
  W = ( log(RRu) / qnorm(1-q) )^2
  return(W)
}
# -------------------------------------------------------------------------- # 


# -------------------------------------------------------------------------- #
# Effect-MAF Dependence

# Mlo -> the MAF of rare SNPs; this MAF has a corresponding RRlo, above which
# we observe snps with probability q. in other words, we want this to be a boundary
# on the MAF and corresponding RR. So if we conclude that for SNPs with very low MAF,
# 95% of RR will be <= RRlo
getW.MAFdep <- function(dat, q)
{
  # trunctate date by MAF threshold
  dat.R <- dat[which(dat$MAF <= MAF.rare.thresh),]
  dat.NR <- dat[which(dat$MAF > MAF.nrare.thresh),]
  
  # i think these are determined once and applied to all SNPs
  # find the 95th percentile RR for rare snps and non-rare snps
  RRlo <- RRplot(dat.R, "rareSNP.png")
  RRhi <- RRplot(dat.NR, "nonerareSNP.png")
  print(paste("RRlo=", RRlo))
  print(paste("RRhi=", RRhi))
  
  # get MAFs (approximately) corresponding to this threshold
  Mlo <- mean(dat$MAF[which(round(dat$RR,3) == round(RRlo, 3))])
  Mhi <- mean(dat$MAF[which(round(dat$RR,3) == round(RRhi, 3))])
  
  Wlo <- ( log(RRlo) / qnorm(1-q) )^2
  Whi <- ( log(RRhi) / qnorm(1-q) )^2
  
  d1 <- (log(Wlo) - log(Whi)) / (Mhi - Mlo)
  d0 <- Wlo * exp(d1 * Mlo)
  
  W <- d0 * exp(-d1 * dat$MAF)
  return(W)
}
# -------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------------- #
# -------------------------------- main -------------------------------------------- #
# ---------------------------------------------------------------------------------- #

# ---
# constants:
# snps have RR above RRhi with probability q
q = 0.05

# threshold below which MAFs are considered rare
MAF.rare.thresh <- 0.01

# threshold above which MAFs are considered non-rare
MAF.nrare.thresh <- 0.1 
# ---



# input should be:
# MAF; Beta; se(Beta); p-value
# extract this from SNPtest.out
# this will depend on snptest version
dat <- read.table("C:/Users/jpluta/Desktop/bayes/snptest_chr1.out", header=TRUE,
                  colClasses= c("character", "character", "integer", "numeric", "numeric", "numeric", "numeric", "numeric"),
                  col.names = c("CHR", "rsID",  "BP","MAF", "p", "info", "B", "B.se"))
                 

# set up some vars
# variance of log(RR)
dat$V <- dat$B.se^2
dat$RR <- exp(dat$B)

# simple QC
# some snps are clearly erroneous, with RR <0 RR huge RR that throws everything off
# so only retain those in some sensible range
# threshold on RR rather than log(RR) because RR is more interpretable
ind <- which(dat$RR > 5 | dat$RR < 0)
dat <- dat[-ind,]

# for  Effect-MAF dependence or Effect-MAF-independence, get W 
dat$W <- getW.MAFind(dat, q )


dat$WBF <- getWBF(dat$V, dat$W, dat$B )
dat$WBF_alt <- 1 / dat$WBF

out <- dat[with(dat, order(dat$WBF_alt, decreasing=TRUE)),]


# make credible set
WBF.sum <- sum(out$WBF_alt)

# .95 creidble set
thresh <- .95 * WBF.sum
i=1
sum=0

while( sum < thresh )
{
  sum = sum + out$WBF_alt[i]
  i = i + 1
}

# variance plot under additive model
MAF <- .1
n <- seq(100,10000,100)
rt.V <- ( n * MAF * (1-MAF))^(-.5)
plot(n, rt.V)
# ---------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------- #




# REMEMBER!
# log(RR) +/- 1.96*se
# beta is log(RR); convert to RR
# also create confidence interval
# this is the right way to get confidence interval of RR, i checked
# be careful about using log(RR) versus (RR)
# annotate CAREFULLY
# dat$RR.hi <- exp(dat$B + 1.96 * dat$se.B)
# dat$RR.lo <- exp(dat$B - 1.96 * dat$se.B)
