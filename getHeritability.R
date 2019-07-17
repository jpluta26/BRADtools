getHeritability <- function( B, maf, lambda = 4)
# input: B, the odds ratio or log odds ratio
#        maf, the minor allele frequency
#        lambda, estimate of genetic + environmental variance. for TGCT, estimated at lambda = 4 for fathers-sons,
#             8 for siblings
# output: proportion of heritability
{
  if( !any(B < 0))
  {
    print("you are probably using OR, converting to logOR")
    B <- log(B)
  }
  
  h = sum( (B * B * 2 * maf * (1 - maf)) / log( lambda^2 ))
  return(h)
}
