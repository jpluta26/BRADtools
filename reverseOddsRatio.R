#!/usr/bin/env Rscript
# pluta 12/9/19

# given an odds ratio with confidence interval, invert the odds ratio
# this is the same as changing the reference allele

args = commandArgs(trailingOnly = TRUE)
if( length(args) < 3)
{
  print("need to provide 3 arguments")
  print("1. OR")
  print("2. Confidence interval - Low")
  stop("3. Confidence interval - High")
}

OR <- args[1]
CI.lo <- args[2]
CI.hi <- args[3]

getSEfromCI <- function(OR, ci.lo, ci.hi)
  #  input: OR, the odds ratio estimated by logistic regression
  #   ci.lo: the low end of the OR
  #   ci.hi: the high end of the OR
  #  output: se, the standard error of the coefficient
  #   this is the standard error of the original estimated coefficient,
  #   eg logOR
{
  
  # CI in logistic regression is asymmetrical
  # take the average for most stable estimate
  se.1 <- abs( (log(ci.lo) - log(OR)) / 1.96 )
  se.2 <- abs( (log(ci.hi) - log(OR)) / 1.96 )
  se <- mean(se.1, se.2)
  
  return(se)
}

get.OR.CI <- function( logOR, se )
  # input:
  #   logOR - the mean of the coefficient estimate (logOR)
  #   se - the standard error of the coefficient estimate 
  # output:
  #   (text), the odds ratio and corresponding 95% confidence interval
{
  CI.lo <- round(exp( logOR - se * 1.96), 3)
  CI.hi <- round(exp( logOR + se * 1.96), 3)
  
  return( list(OR, CI.lo, CI.hi))
}


se <- getSEfromCI( OR, CI.lo, CI.hi )
logOR <- -log(OR)
get.OR.CI( logOR, se )
