setwd("C:/Users/jpluta/Desktop/TECAC-fulldata/PRS")
library(dplyr)

# pluta 11/12/19
# 
# see: Choi et al., 2018. A guide to performing polygenic risk score analyses.

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
  
  # same model as in tecac association
  fit <- glm( data = dat, PHENO ~ PHENO.q + EV1 + EV2 + EV3 + site, family = "binomial")
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








# --------------------- main -------------------------------------- #
# covariates, phenotypes, and prs scores from plink
cov <- read.table("tecac.cov", header = TRUE, as.is = T)
PRS <- read.table("prs.scores", header = TRUE, as.is = T)

if( any(PRS$PHENO == -9))
{ PRS <- PRS[-which(PRS$PHENO == -9),] }

ind <- match(cov$IID, PRS$IID)
PRS <- PRS[ind,]

# setup data.frame
cov$PRS <- PRS$SCORESUM
cov$PHENO <- PRS$PHENO - 1

# divided into 100 quantiles and choose by range
cov$quantile <- ntile(cov$PRS, 100)


median.range <- c(45:55)

# data frame of ranges (defines 20 quantiles)
r.df <- data.frame(range.lo = seq(from = 1, to = 96, by = 5),
                       range.hi = seq(from = 5, to = 100, by = 5) )

# -- compare each quantile to the median and convert to OR
out <- c()
for( i in 1:20 )
{
  out <- rbind(out, unlist(compareQuantileToMedian( cov, median.range, seq(from = r.df[i,1], to = r.df[i,2], by = 1))))
}

OR.out <- convertBetaToOR( out[,1], out[,2] )
OR.df <- data.frame( OR = OR.out[[1]], CI.lo = OR.out[[2]], CI.hi = OR.out[[3]])
# --


# plot
p1 <- ggplot( data = OR.df ) + 
  geom_segment(aes( x = seq(1:20), xend = seq(1:20), y = OR.df$CI.lo, yend = OR.df$CI.hi), color = "blue") +
  geom_point(aes(x = seq(1:20), y = OR.df$OR), color = "blue") + 
  theme_minimal() + xlab("Odds Ratio") + ylab("Quantiles") +
  scale_x_continuous(breaks = seq(1:20), labels = paste(r.df[,1], r.df[,2], sep = ":") )
print(p1)
# ----------------------------------------------------------------------- #
