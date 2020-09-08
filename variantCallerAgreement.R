
# pluta 9/8/20

library(psych)
dat <- read.table("~/Downloads/Somatic_calls.test_data.csv", header = T, sep = ",")

# function to compute all pairwise comparisons of callers
getAllPairwiseComp <- function( dat )
{
  out <- list()
  pairwise.comp <-  t(combn(seq(1:4), 2)) + 5
  for( i in 1:dim(pairwise.comp)[1])
  {
    x = pairwise.comp[i,1]
    y = pairwise.comp[i,2]
    tmp <- cohen.kappa( dat[, c(x,y)])$agree
    names(dimnames(tmp)) <- c(colnames(dat)[x], colnames(dat)[y])
    out[[i]] <- round(tmp,4)
  }
  
  return(out)
}

# function to check if any value is "NO_CALL"
nocall_any <- function( dat )
{
  if( any(dat == "NO_CALL"))
  {
    return(TRUE)
  }
  
  return(FALSE)
}

# all data
dat.all <- getAllPairwiseComp(dat)

# only data with lof == 1
dat.lof1 <- getAllPairwiseComp(dat[ dat$LOF.level == 1,])

# only data without NO_CALL
dat.callonly <- getAllPairwiseComp(dat[ !apply(dat[,6:9], 1, nocall_any),])




