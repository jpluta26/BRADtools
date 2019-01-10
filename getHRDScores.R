

# john pluta & kara maxwell 
# v1.00: 11/7/2016
# v1.01: 11/8/2018
# citation: Maxwell et al. BRCA locus specific loss of heterozygosity in germline BRCA1 and BRCA2 carriers. 2017. Nat Comm 8(1):319.
# contact: jpluta@pennmedicine.upenn.edu

# ===================================================================================== #
# ================================ constants ========================================== #
# ===================================================================================== #
# the columns needed for an analysis, per standard sequenza naming conventions
seq.cols.needed = c("chromosome", "start.pos", "end.pos", "CNt", "A", "B")

# predefined data about chromosome size, centromere, and telomere location
# chromosome size, centromere and telomere locations (in hg19/GRCh37)
# TODO: add ref data for GRCH38
ref.dat = data.frame( chromosome = paste("chr", c(seq(1:22), "X", "Y"), sep=""),
                      centromere.start = c(121535434, 92326171, 90504854, 49660117, 46405641, 58830166,
                                           58054331, 43838887, 47367679, 39254935, 51644205, 34856694, 
                                           16000000, 16000000, 17000000, 35335801, 22263006, 15460898,
                                           24681782, 26369569, 11288129, 13000000, 58632012, 10104553),
                      centromere.end = c(124535434, 95326171, 93504854, 52660117, 49405641, 61830166,
                                         61054331, 46838887, 50367679, 42254935, 54644205, 37856694,
                                         19000000, 19000000, 20000000, 38335801, 25263006, 18460898,
                                         27681782, 29369569, 14288129, 16000000, 61632012, 13104553),
                      p.telomere.end = rep(10000, 24),
                      q.telomere.start = c(249240621, 243189373, 198012430, 191144276, 180905260,
                                           171105067, 159128663, 146354022, 141203431, 135524747,
                                           134996516, 133841895, 115159878, 107339540, 102521392,
                                           90344753, 81185210, 78067248, 59118983, 63015520,
                                           48119895, 51294566, 155260560, 59363566),
                      chr.size = c(249250621, 243199373,  198022430, 191154276, 180915260,
                                   171115067, 159138663, 146364022, 141213431, 135534747,
                                   135006516, 133851895, 115169878, 107349540, 102531392,
                                   90354753, 81195210, 78077248, 59128983, 63025520,
                                   48129895, 51304566, 155270560, 59373566)
                      
)
# ===================================================================================== #
# ===================================================================================== #
# ===================================================================================== #



# ------------------------------------------------------------------------------------- #
# ----------------------------------- functions --------------------------------------- #
# ------------------------------------------------------------------------------------- #


# ------------------------------- getTAI.raw ----------------------------------------- #
# function to get telomeric allelic imbalance (TAI) without including main CNt segments
getTAI.raw <- function(seq.dat, min.seg.size = 11e06)
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  #   min.seg.size, (integer), the minimum segment size required for analysis
  # output:
  #   HRD.TAIm (numeric)
{
  
  # seq.dat$s: length of the segment
  # if length of segment is greater than the minimum; and allelic imbalance is present; and
  # segment is in the telomeres; and is not on centromere...
  HRD.TAI <- sum(seq.dat$s > min.seg.size & seq.dat$AI & !seq.dat$cross.arm 
                 & (seq.dat$post.telomere | seq.dat$pre.telomere))

  return(HRD.TAI)
}
# ------------------------------------------------------------------------------- #


# ------------------------------- getTAI.norm ----------------------------------- #
# function to get telomeric allelic imbalance (TAI) without including main CNt segments
getTAI.norm <- function(seq.dat, CN.dat, min.seg.size = 11e06)
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  #   min.seg.size, (integer), the minimum segment size required for analysis
  # output:
  #   HRD.TAIm (numeric)
{
  
  # seq.dat$s: length of the segment
  # if length of segment is greater than the minimum; and allelic imbalance is present; and
  # segment is in the telomeres; and is not on centromere...

  # create an index of main.CN segments; these get removed to normalize TAI
  rm.ind <- c()
  
  for( i in 1:length(CN.dat$chromosome))
  {
    rm.ind <- c(rm.ind, which(seq.dat$chromosome == CN.dat$chromosome[i] & seq.dat$CNt == CN.dat$main.CN[i]))
  }
  
  # normalized
  HRD.TAI <- getTAI.raw(seq.dat[-rm.ind,], 11e06)
  
  return(HRD.TAI)
}
# ------------------------------------------------------------------------------- #


# ------------------------------- getNTAI --------------------------------------- #
# function to get NTAI without including main CNt segments
# non-telomeric allelic imbalance
getNTAI.raw <- function(seq.dat, min.seg.size = 11e06)
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  #   min.seg.size, (integer), the minimum segment size required for analysis
  # output:
  #   HRD.NTAIm (numeric)
{
  # if length of segment is greater than the minimum; and allelic imbalance is present; and
  # segment is not in the telomeres or centromere...
  HRD.NTAI <- sum(seq.dat$s > min.seg.size & 
                    seq.dat$AI & 
                    !seq.dat$cross.arm & 
                    !seq.dat$post.telomere & 
                    !seq.dat$pre.telomere)

  return(HRD.NTAI)
}
# --------------------------------------------------------------------------------- #


# ------------------------------- getNTAI.norm --------------------------------------- #
# function to get NTAI without including main CNt segments
# non-telomeric allelic imbalance
getNTAI.norm <- function(seq.dat, CN.dat, min.seg.size = 11e06 )
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  #   min.seg.size, (integer), the minimum segment size required for analysis
  # output:
  #   HRD.NTAIm (numeric)
{
  # create an index of main.CN segments; these get removed to normalize TAI
  rm.ind <- c()
  
  for( i in 1:length(CN.dat$chromosome))
  {
    rm.ind <- c(rm.ind, which(seq.dat$chromosome == CN.dat$chromosome[i] & seq.dat$CNt == CN.dat$main.CN[i]))
  }
  
  # normalized
  HRD.NTAI <- getNTAI.raw(seq.dat[-rm.ind,], min.seg.size)
  
  return(HRD.NTAI)
}
# --------------------------------------------------------------------------------- #


# ------------------------------- getLOH --------------------------------------- #
# function to get LOH 
getLOH <- function(seq.dat)
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  # output:
  #   HRD.LOH (numeric)
{
  # HRD-LOH calculation
  # loss of heterozygosity
  # if B == 0 & s > 15mbp & within chromosome, then HRD-LOH is TRUE
  
  HRD.LOH <- sum( (seq.dat$s > 15e06) & ((seq.dat$s / seq.dat$chr.size) < 0.9) & 
                    (seq.dat$B == 0) )
  # mdacc doesnt have seq.dat$B == 0 condition
  return(HRD.LOH)
}
# --------------------------------------------------------------------------------- #



# ---------------------------------- getLST.raw ----------------------------------- #
getLST.raw <- function(seq.dat)
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  # output:
  #   HRD.LST (numeric), large state transition
{
  
  # length of gap between segments
  seq.dat$brk.len <- 0
  seq.dat$LST <- FALSE
  
  
  n.segs <- dim(seq.dat)[1]
  
  # if the gap between 2 segments is < 3mbp, and each adjacent segment is > 10mbp, and the segment does not
  # cross the centromere, it is an LST
  for( i in 1:(n.segs - 1))
  {
    seq.dat$brk.len[i + 1] <- seq.dat$start.pos[i + 1] - seq.dat$end.pos[i]
    seq.dat$LST[i] <- seq.dat$brk.len[i] < 3e06 &
      seq.dat$cross.arm[i] == FALSE & 
      seq.dat$s[i] > 10e06 & 
      seq.dat$s[i + 1] > 10e06
  }
  
  HRD.LST <- sum(seq.dat$LST)
  return(HRD.LST)
  
}
# ------------------------------------------------------------------------------- #



# ---------------------------------- getLST.norm --------------------------------------- #
getLST.norm <- function(seq.dat, ploidy.dat)
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  # output:
  #   HRD.LST (numeric), large state transition
{
  
  # any large genomic rearrangement could be increased simply by having more 
  # chromosomes (higher ploidy), and not because of the biological process
  # adjust for this
  HRD.LSTm = getLST.raw(seq.dat) - (15.5 * ploidy.dat$ploidy.estimate[2])
  
}
# ------------------------------------------------------------------------------- #



# ------------------------------------------------------------------------------------- #
# some long segments are falsely read as multiple shorter segments
# recombine these into one segment
combineSeg <- function(seq.dat, max.brk.len)
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  #   max.brk.len (integer), the maximum size of the break (gap between two segments)
  #       combine segments where break.len < max.brk.len
  # output:
  #   HRD.LOH (numeric)  
  
{
  seq.dat$brk.len <- 0
  
  n.segs <- dim(seq.dat)[1]
  rm.ind <- c()
  
  for( i in 1:(n.segs - 1))
  {
    seq.dat$brk.len[i + 1] <- seq.dat$start.pos[i + 1] - seq.dat$end.pos[i]
    
    
    # if break length is < max.brk.len; & CnT1 == CnT2 & A1 == A2 & B1 == B2
    # then combine the two segments
    if( seq.dat$brk.len[i + 1] < max.brk.len & 
          seq.dat$brk.len[i + 1] > 0 &
          seq.dat$CNt[i + 1] == seq.dat$CNt[i] &
            seq.dat$A[i + 1] == seq.dat$A[i] &
              seq.dat$B[i + 1] == seq.dat$B[i]   )
    {
 
      seq.dat$start.pos[i + 1] <- seq.dat$start.pos[i]
      rm.ind <- c(rm.ind, i)
    }
  }
  
  if( !is.null(rm.ind) )
  {
    print(rm.ind)

    seq.dat <- seq.dat[-rm.ind,]
  }
  
  
  print(paste("Initial data has ", n.segs, " segments.", sep = ""))
  print(paste("Using break-length threshold of ", max.brk.len, sep = ""))
  print(paste("Combined data has ", dim(seq.dat)[1], " segments.", sep = ""))
  return(seq.dat)
}
# ------------------------------------------------------------------------------------- #



# ------------------------------- preprocessSeq ---------------------------------------- #
# define allelic imbalance (AI), telomere positions, segment length, cross arms
# these are used in the various HRD scores
preprocessSeq <- function( seq.dat )
# input:
#   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
# output:
#   seq.dat (data.frame), the sequencing data with more factors added
{

  # check the sequencing data input to make sure it has the data we need
  if( any(!(seq.cols.needed %in% colnames(seq.dat))))
  {
    print(paste("column", 
                seq.cols.needed[which(!(seq.cols.needed %in% colnames(seq.dat)))],
                "missing in the seq data.", sep=" "))
    print(paste("Looking for columns:", seq.cols.needed, sep=" "))
    print(paste("Found:", seq.cols.needed[seq.cols.needed %in% colnames(seq.dat)]))
    stop("Column mismatch. Exiting.")
  } 
  
  # remove rows w/ missing entries
  if( any( is.na(seq.dat$CNt) ))
  {
    print(paste("Row: ", which(is.na(seq.dat$CNt)), " contains NA; removing.", sep = ""))
    seq.dat <- subset(seq.dat, !is.na(seq.dat$CNt))
  }
  
  levels(seq.dat$chromosome) <- levels(ref.dat$chromosome)
  
  seq.dat$frac.chr <- (seq.dat$end.pos - seq.dat$start.pos) / ref.dat$chr.size[match(seq.dat$chromosome, ref.dat$chromosome)]
  
  # segment length
  seq.dat$brk.len <- 0
  seq.dat <- combineSeg(seq.dat, 3e06)
  
  # matches reference data to the correct chromosome in the subject data
  key = match(seq.dat$chromosome, ref.dat$chromosome)
  
  seq.dat$s <- seq.dat$end.pos - seq.dat$start.pos
  seq.dat$chr.size <- ref.dat$chr.size[match(seq.dat$chromosome, ref.dat$chromosome)]
  
  seq.dat$AI <- (seq.dat$A > seq.dat$B) & (seq.dat$A != 1) & (seq.dat$B != 0)
  
  seq.dat$start.arm <- seq.dat$end.arm <- rep("NA", dim(seq.dat)[1])
  
  ind = seq.dat$start.pos - ref.dat$centromere.start[key] > 1
  seq.dat$start.arm[ind]  = "p"
  seq.dat$start.arm[!ind] = "q"
  
  ind = seq.dat$end.pos - ref.dat$centromere.end[key] > 1
  seq.dat$end.arm[ind]  = "p"
  seq.dat$end.arm[!ind] = "q"
  
  seq.dat$cross.arm <- seq.dat$start.arm != seq.dat$end.arm
  
  
  # define end telomeres
  seq.dat$post.telomere <- (seq.dat$start.pos - ref.dat$p.telomere.end[key] <= 1000)
  seq.dat$pre.telomere  <- (ref.dat$q.telomere.start[key] - seq.dat$end.pos <= 1000)
  # ---
  
  return(seq.dat)
}
# ------------------------------------------------------------------------------------- #



# ----------------------------------- getCNT ------------------------------------------ #
# setup CNt values
getCNt <- function( seq.dat )
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  # out:
  #   CN.dat (data.frame), copy number data
{
  main.CN <- rep(0, length(unique(seq.dat$chromosome)))
  
  ct1 <- 0
  for( chr in unique(seq.dat$chromosome))
  {
    ct1 <- ct1 + 1
    
    # reduce data to one chromosome at a time
    dat <- seq.dat[seq.dat$chromosome == chr,]
    
    # remove any NAs or this will crash
    dat <- dat[!is.na(dat$CNt),]
    CNt.vals <- rep(0, length(unique(dat$CNt)))
    ct2 <- 0
    
    # store CNt values per CNt
    for( CNt in unique(dat$CNt) )
    {
      ct2 <- ct2 + 1
      CNt.vals[ct2] <- sum(dat$frac.chr[dat$CNt == CNt])
    }
    
    # record values + associated CNt
    out <- data.frame(CNt=unique(dat$CNt), CNt.vals=CNt.vals)
    
    # get one main.CN per chromosome
    main.CN[ct1] <- out$CNt[which(CNt.vals == max(CNt.vals))]
    
  }
  
  CN.dat <- data.frame(chromosome=unique(seq.dat$chromosome), main.CN=main.CN)
  return(CN.dat)
}
# ------------------------------------------------------------------------------- #



# ------------------------------------- hrd.stats ------------------------------------- #
# hrd.stats is a function to compute the three HRD metrics (HRD-LOH, HRD-NTAI, and
# HRD-LST), as well as total HRD and mean HRD. this simply wraps up the output in one
# dataframe

#
# input: seq.dat, (data.frame) with chromosome, start.pos, end.pos, CNt, alleleA, alleleB;
#         ploidy.dat (data.frame), the ploidy data
#         CN.dat (data.frame), copy number data
#         min.seg.size (integer), minimum segment size used in TAI calculations
#         scaleTotal (boolean), rescale HRD total to 0-100
# output: out, a data.frame with HRD metrics
hrd.stats <- function(seq.dat, ploidy.dat, CN.dat, min.seg.size = 11e06, scaleTotal = FALSE)
{
  
  seq.dat <- preprocessSeq(seq.dat)
  
  # raw data
  HRD.NTAIr <- getNTAI.raw( seq.dat, min.seg.size )
  HRD.TAIr  <- getTAI.raw(  seq.dat, min.seg.size )
  HRD.LSTr  <- getLST.raw(  seq.dat )
  HRD.LOH   <- getLOH(  seq.dat )
  
  HRD.NTAIm <- getNTAI.norm( seq.dat, CN.dat, min.seg.size )
  HRD.TAIm  <- getTAI.norm(  seq.dat, CN.dat, min.seg.size )
  HRD.LSTm  <- getLST.norm( seq.dat, ploidy.dat )
  
  out = data.frame(HRD.LOH   = HRD.LOH, 
                   HRD.TAIr  = HRD.TAIr,
                   HRD.TAIRm = HRD.TAIm,
                   HRD.NTAIr = HRD.NTAIr,
                   HRD.NTAIm = HRD.NTAIm,
                   HRD.LSTr  = HRD.LSTr,
                   HRD.LSTm  = HRD.LSTm )
  
  out$HRD.TOTAL <- out$HRD.LOH + out$HRD.LSTr + out$HRD.NTAIm
  
  # rescale HRD total to 0-100
  if( scaleTotal == TRUE )
  {
    library(scales)
    out$HRD.TOTAL <- round( rescale(out$HRD.TOTAL, to = c(0,100), from = range(out$HRD.TOTAL)) )
  }
  
  return(out)
  
}
# ------------------------------------------------------------------------------- #

# re: min seg length-
# timms et al 2014. -"HRD-TAI score was defined as the number of regions with allelic imbalances that extend to
# one of the subtelomeres but do not cross the centromere. A region was counted only if it
# encompassed a certain minimum number of SNPs (on average approximately 1.8 Mb). We tested for 
# association of HRD-TAI score with BRCA1, BRCA2, and RAD51C deficiency in three datasets of 609 ovarian
# tumors (data not shown) and found the association to be more significant if the cutoff size of
# HRD-TAI regions was increased to 11MB. Therefore, a modified HRD-TAIm score was defined as the
# number of regions with allelic imbalances that extend to one of the subtelomeres, do not cross the 
# centromere, and are longer than 11Mb." 
#



# ------------------------------------------------------------------------------------- #
# ----------------------------------- end functions ----------------------------------- #
# ------------------------------------------------------------------------------------- #

# example usage
# seq.dat <- read.table("sub01_seqz_segments.txt", header = TRUE)
# ploidy.dat <- read.table("sub01.seqz_confints_CP.txt", header = TRUE)
# sub.id <- "sub01"
# seq.dat <- preprocessSeq(seq.dat)
# CN.dat <- getCNt(seq.dat)
# hrd.dat <- round(hrd.stats(seq.dat, ploidy.dat, CN.dat), 3)
# hrd.dat$ID <- sub.id

