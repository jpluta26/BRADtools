
# pluta 5/7/21

# set of functions for aligning snps between two datasets



# ---------------------------------------------------------------------- #
checkGenoFlip <- function(bim1, bim2)
# compare genotypes from two datasets to see if they match
# input: bim1, bim2 (data.frames)- bim files from plink
# output: boolean- TRUE if both alleles match
{
  # if the same alleles are in both studies but in reverse order, need to flip
  if( sum( c(bim1$V5, bim1$V6) %in% c(bim2$V5, bim2$V6)) == 2)
  {
    if( bim1$V5 == bim2$V6 )
    {
      return(TRUE)
    }
  } 
  return(FALSE)
}
# ---------------------------------------------------------------------- #


# ---------------------------------------------------------------------- #
flipGeno <- function( geno.vec )
# if change minor allele to major allele, also flip geno types
# 0 becomes 2, 2 becomes 0
# input: geno.vec (integer), vector of genotype values (0,1,2)
# output: geno.vec (integer), vector of genotypes with 0 and 2 reversed
{
  geno.vec[ geno.vec == 0] <- -9
  geno.vec[ geno.vec == 2] <- 0
  geno.vec[ geno.vec == -9] <- 2
  return(geno.vec)
}
# ---------------------------------------------------------------------- #

# --------------------------- flipAllele ------------------------------------------- #
flipAllele <- function( allele )
# input: allele (character), one of AGCT
# output: allele (character), one of AGCT, flipped from input
{
  
  allele <- as.character(allele)
  
  if( allele == "A")
  {
    return("T")
  } else
    if( allele == "G")
    {
      return("C")
    } else
      if( allele == "C" )
      {
        return("G")
      } else
        if( allele == "T")
        {
          return("A")
        } else
          return(allele)
}
# ---------------------------------------------------------------------------------- #


# ------------------------------------- getCommonSnps ----------------------------------------- #
getCommonSnps <- function( bim.list )
  # function to get snps present in all studies in bim.list
  # input: bim.list (list), the list of bim files
  #
  # output: common.snps (character vector), the list of common snps
{
  for( i in 1:length(bim.list))
  {
    if(i == 1)
    {
      common.snps <- bim.list[[i]]$V2
    } else
    {
      snps <- bim.list[[i]]$V2
      common.snps <- intersect(common.snps, snps)
    }
  }
  
  return(common.snps)
}
# --------------------------------------------------------------------------------------------- #


# ------------------------------- reduceBEDtoCommonSet ---------------------------------------- #
reduceBEDToCommonSet <-  function( bed.list, common.snps )
  # function to reduce bed files to common sets
  # input: bed.list (list), list of bed files to align
  #        common.snps (character vector), list of snps common to all data sets
  #
  # output: bed.list (list), list of aligned bed files
{
  for( i in 1:length(bed.list))
  {
    bed.list[[i]] <- bed.list[[i]][ ,colnames(bed.list[[i]])  %in% common.snps ]
  }
  
  # if the dimensions of bed files aren't the same, something went wrong
  n.snps <- unlist(lapply(bed.list, dim))[ seq( from = 2, to = length(bed.list) * 2, by = 2)]
  
  if( length(unique(n.snps)) !=  1)
  {
    stop("all bed files do not have equal dimensions")
  }
  
  return(bed.list)
}
# --------------------------------------------------------------------------------------------- #


# ------------------------------ reduceBIMToCommonSet ----------------------------------------- #
reduceBIMToCommonSet <-  function( bim.list, common.snps )
  # reduce BIM files to the set of common snps
  # input: bim.list (list), list of bim files
  #        common.snps (character vector), the snps common to all sets
  #
  # output: bim.list (list), the list of aligned bim files
{
  for( i in 1:length(bim.list))
  {
    bim.list[[i]] <- bim.list[[i]][ bim.list[[i]]$V2  %in% common.snps, ]
  }
  
  n.snps <- unlist(lapply(bim.list, dim))[ seq( from = 2, to = length(bim.list) * 2, by = 2)]
  
  if( length(unique(n.snps)) !=  1)
  {
    stop("all bed files do not have equal dimensions")
  }
  return(bim.list)
}
# --------------------------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------------- #
alignSnps <- function( snp, BIM1, BIM2,  geno.dat )
# align snps in BIM2 to BIM1, flip corresponding genotypes in geno.dat
# BIM1, BIM2 (data.frame), bim files from plink
# geno.dat (data.frame), genotype data imported by BEDMatrix
{
 
  
  
 
  BIM1$V2 <- as.character(BIM1$V2)
  BIM2$V2 <- as.character(BIM2$V2)
  
 
  bim1 <- BIM1[ which(BIM1$V2 == snp),]
  bim2 <- BIM2[ which(BIM2$V2 == snp),]
  geno <- geno.dat[,which(colnames(geno.dat) == snp)]
  
  # align alleles
  # case one, the ref and alt allele are flipped in direction
  # eg C T  and T C
  
  if( checkGenoFlip(bim1, bim2))
  {
   
    geno <- flipGeno(geno)
  }
  
  # case two, the ref and alt allele are flipped in both name and/or direction
  # eg C T  and G A
 
  if( sum( c(bim1$V5, bim1$V6) %in% c(bim2$V5, bim2$V6)) == 0)
  {
   
    bim2$V5 <- flipAllele(bim2$V5)
    bim2$V6 <- flipAllele(bim2$V6)
 
    if( checkGenoFlip(bim1, bim2) )
    {
     
      geno <- flipGeno(geno)
    }
  }
  
  return(geno)
}
# ---------------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------------- #
alignBIM <- function(snp, BIM1, BIM2)
  # originally intended for writing bim/bed as output to only need to do alignment once
  # havent gotten that to work
  # aligns BIM2 to BIM1
  # if this crashes its probably because a snp isnt found
{
  
  
  BIM1$V2 <- as.character(BIM1$V2)
  BIM2$V2 <- as.character(BIM2$V2)
  
  bim1 <- BIM1[ which(BIM1$V2 == snp),]
  bim2 <- BIM2[ which(BIM2$V2 == snp),]
  
  # align alleles
  # case one, the ref and alt allele are flipped in direction
  # eg C T  and T C
  if( bim1$V5 ==  bim2$V6 )
  {
    bim2$V6 <- bim1$V5
    bim2$V5 <- bim1$V6
  }
  
  # case two, the ref and alt allele are flipped in both name and/or direction
  # eg C T  and G A
  
  if( sum( c(bim1$V5, bim1$V6) %in% c(bim2$V5, bim2$V6)) == 0)
  {
    
    bim2$V5 <- flipAllele(bim2$V5)
    bim2$V6 <- flipAllele(bim2$V6)
    
    if( bim1$V5 ==  bim2$V6 )
    {
      bim2$V6 <- bim1$V5
      bim2$V5 <- bim1$V6
    }
  }
 
  return(bim2)
}
# ---------------------------------------------------------------------------------- #





# sample call:
#test = lapply(colnames(dat1[[1]][,1:10]), alignSnps, BIM1, BIM2, dat1[[1]])
#out <- matrix(unlist(test), nrow=226,ncol=10)