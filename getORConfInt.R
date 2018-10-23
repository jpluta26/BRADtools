# pluta 10/22/18

# function to return the odds ratio and confidence interval
# from logOR and its standard error
get.OR.CI <- function( OR, se )
# input:
#   OR - the mean of the coefficient estimate (logOR)
#   se - the standard error of the coefficient estimate 
# output:
#   (text), the odds ratio and corresponding 95% confidence interval
{
  CI.lo <- round(exp( OR - se * 1.96), 3)
  CI.hi <- round(exp( OR + se * 1.96), 3)
  
  return( paste("OR: ", round(exp(OR),3), "CI: (", CI.lo, ", ", CI.hi, ")", sep = ""))
}
