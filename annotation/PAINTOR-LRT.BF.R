# simple function to run the likelihood ratio test between the baseline
# pluta 11/4/19

args = commandArgs(trailingOnly = TRUE)
if( length(args) < 1)
{
  stop()
}

featureName = args[1]

INFILE <- paste("BF", featureName, sep = ".")

logBF.M1 <- as.numeric(read.table(INFILE, header = F)$V1)
logBF.M0 <- as.numeric(read.table("BF.Base", header = F)$V1)

LRT = -2 * (logBF.M0 - logBF.M1)
p = 1 - pchisq(LRT, 1)

print(paste(featureName, p, sep = ","))
