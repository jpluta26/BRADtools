
#!/usr/bin/env Rscript
# pluta 3/4/25
# script to create quantile plots from PRS scores

# --------------------------------- parse input ----------------------------- #
library(dplyr)
library(ggplot2)
suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c("-s", "--SCORES"), action="store", default = NA, type="character",
              help = "The PRS scores calculated from PLINK; assuming SUM of scores"),
  make_option(c("-p", "--PHENO"), action = "store", default = NA, type = "character",
              help = "file containing phenotype and covariate information; IDs must match those in SCORES."),
  make_option(c("-o", "--OUTNAME"), action = "store", default = NA, type = "character",
              help = "prefix attached to output")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if( (opt$help))
{
  print_help(opt_parser)
  stop()
}

SCORES = opt$SCORES
PHENO = opt$PHENO
OUT = opt$OUT
# ------------------------------------------------------------------------------- #


# ======================= functions ========================================== #

# !!!!! for now covariates and covar names are hard coded
# TODO: fix this
# ------------------- compareQuantileToMedian --------------------- #
compareQuantileToMedian <- function( dat, median, quantile )
  # input: dat (data.frame), the df containing the PRS scores, covariates, pheno type, and
  #     quantiles
  #        median (integer vector), the range of quantiles to use as the median
  #        quantile (integer vector), the range of quantiles to compare against the median
  # output: a list with the Beta coefficient, which is the logOR of quantile vs median
  #     also returns standard error of beta
{
  # truncate data to just the subjects in the median and quantile groups
  dat <- dat[ dat$quantile %in% c(median, quantile),]
  
  if( dim(dat)[1] == 0 )
  {
    stop("dat is empty after subsetting, check your median/quantile definitions")
  }
  
  # subjects in the median quantiles are coded 0, subjects in the comparison quantile are 1
  dat$PHENO.q <- 0
  dat$PHENO.q[ dat$quantile %in% quantile ] <- 1
  
  dat$Site <- as.factor(dat$Site)
  dat[,grep("EV", colnames(dat))] <- apply(dat[,grep("EV", colnames(dat))], 2, as.numeric)
  
  # same model as in tecac association
  fit <- glm( data = dat, as.integer(phenotype) ~ PHENO.q + EV1 + EV2 + EV3 + Site, family = "binomial")
  B <- summary(fit)$coefficients[2,1]
  B.se <- summary(fit)$coefficients[2,2]
  
  return( list(B, B.se) )
}
# ------------------------------------------------------------------ #

# ------------------- convertBetaToOR ------------------------------ #
convertBetaToOR <- function( B, B.se )
  # input: B (numeric), the beta weight from regression; the logOR of cases v controls
  #        B.se (numeric), standard error of beta
  # output: OR (numeric), the odds ratio
  #         CI.hi/lo (numeric), 95% confidence interval of the odds ratio
{
  OR <- exp(B)
  CI.hi <- exp( B + 1.96 * B.se )
  CI.lo <- exp( B - 1.96 * B.se )
  return( list(OR, CI.lo, CI.hi) )
}
# ------------------------------------------------------------------ #
# ================================================================== #



# ======================== MAIN ==================================== #
median.range <- c(45:55)
r.df <- data.frame(range.lo = seq(from = 1, to = 96, by = 5),
                   range.hi = seq(from = 5, to = 100, by = 5) )

print("reading in data...")
scores <- read.table(SCORES, header=TRUE)

# .sample format has an extra header line
fam <- read.table(PHENO, header=TRUE)
if( strsplit(PHENO, "[.]"[[1]][1]) == ".sample")
{
  fam <- fam[-1,  ]
}

# make sure scores and phenotypes have the same ids in the same order
# TODO: generalize the input names
fam$ID_2 <- paste(fam$ID_2, fam$ID_2, sep = "_")
fam <- fam[ fam$ID_2 %in% scores$IID,]
fam <- fam[ match(scores$IID, fam$ID_2),]
if(sum(scores$IID==fam$ID_2) != dim(scores)[1])
{
  stop("missing subjects")
}
print("done!")

print(colnames(scores))
fam$SCORE <- scores[,grep("SCORE", colnames(scores))]

print("performing quantile comparison...")
dat <- fam
dat$quantile <- ntile(dat$SCORE, 100)

# define quantiles in 5% intervals
out <- c()
for( i in 1:20 )
{
  out <- rbind(out, unlist(compareQuantileToMedian( dat, median.range, seq(from = r.df[i,1], to = r.df[i,2], by = 1))))
}

OR.out <- convertBetaToOR( out[,1], out[,2] )
OR.df <- data.frame( OR = OR.out[[1]], CI.lo = OR.out[[2]], CI.hi = OR.out[[3]])
print("done!")

print("writing output...")
write.table(OR.df, paste0(OUT, "_ordf.txt"), quote=F, row.names=F, col.names= T)


p1 <- ggplot( data = OR.df ) + 
  geom_segment(aes( x = seq(1:20), xend = seq(1:20), y = CI.lo, yend = CI.hi), color = "blue") +
  geom_point(aes(x = seq(1:20), y = OR), color = "blue") + 
  theme_minimal() + xlab("Quantiles") + ylab("Odds Ratio") +
  scale_x_continuous(breaks = seq(1:20), labels = paste(r.df[,1], r.df[,2], sep = ":") )
print(p1)


png(paste0(OUT, ".png"))
print(p1)
dev.off()
print("done!")
# ================================================================== #
