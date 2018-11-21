

# john pluta & kara maxwell 
# v1.00: 11/7/2016
# v1.01: 11/8/2018
# citation: Maxwell et al. BRCA locus specific loss of heterozygosity in germline BRCA1 and BRCA2 carriers. 2017. Nat Comm 8(1):319.

# example usage:
# sub.id = "subX"
# ploidy.dat <- read.table("653-Brca2Br32-T.seqz_confints_CP.txt", header = T)
# seq.dat <- read.table("653-Brca2Br32-T.seqz_segments.txt",  header = T)

# hrd.dat = getHRD.Data( sub.id, seq.dat, ploidy.dat )


# ===================================================================================== #
# ================================ constants ========================================== #
# ===================================================================================== #
# the columns needed for an analysis, per standard sequenza naming conventions - you can change this
# to whatever colnames you use)
seq.cols.needed = c("chromosome", "start.pos", "end.pos", "CNt", "A", "B")

# predefined data about chromosome size, centromere, and telomere location
# chromosome size, centromere and telomere locations (in hg19/GRCh37)
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


# ------------------------------- getTAI --------------------------------------- #
# function to get TAI without including main CNt segments
# telomeric allelic imbalance
getTAI <- function(seq.dat)
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  # output:
  #   HRD.TAIm (numeric)
{
  
  # seq.dat$s: length of the segment
  HRD.TAI <- sum(seq.dat$s > 11000000 & seq.dat$AI & !seq.dat$cross.arm 
                 & (seq.dat$post.telomere | seq.dat$pre.telomere))
  return(HRD.TAI)
}
# ------------------------------------------------------------------------------- #


# ------------------------------- getNtAI --------------------------------------- #
# function to get NtAI without including main CNt segments
# non-telomeric allelic imbalance
getNtAI <- function(seq.dat)
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  # output:
  #   HRD.NtAIm (numeric)
{
  HRD.NtAI <- sum(seq.dat$s > 11000000 & seq.dat$AI & !seq.dat$cross.arm & !seq.dat$post.telomere & !seq.dat$pre.telomere)
  return(HRD.NtAI)
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
  # if B == 0 & s > 15000mbp & within chromosome, then HRD-LOH is TRUE
  
  HRD.LOH <- sum( (seq.dat$s > 15000000) & ((seq.dat$s / seq.dat$chr.size) < .9) & 
                    (seq.dat$B == 0) )
  
  return(HRD.LOH)
}
# --------------------------------------------------------------------------------- #


# to do: try combining segments
# if break length is < x; & CnT1 == CnT2 & A1 == A2 & B1 == B2
# then combine the two segments

combineSeg <- function(seq.dat, brk.len.thresh)
{
  seq.dat$brk.len <- 0
  
  n.segs <- dim(seq.dat)[1]
  rm.ind <- c()
  
  for( i in 1:(n.segs - 1))
  {
    seq.dat$brk.len[i + 1] <- seq.dat$start.pos[i + 1] - seq.dat$end.pos[i]
    
    if( seq.dat$brk.len[i + 1] < brk.len.thresh & 
          seq.dat$brk.len[i + 1] > 0 &
          seq.dat$CNt[i + 1] == seq.dat$CNt[i] &
            seq.dat$A[i + 1] == seq.dat$A[i] &
              seq.dat$B[i + 1] == seq.dat$B[i]   )
    {
   # nope, wrong.
      seq.dat$start.pos[i + 1] <- seq.dat$start.pos[i]
      seq.dat$s[i + 1] <- seq.dat$end.pos[i + 1] - seq.dat$start.pos[i + 1]
      rm.ind <- c(rm.ind, i)
    }
  }
  
  if( !is.null(rm.ind) )
  {
    print(rm.ind)
    seq.dat <- seq.dat[-rm.ind,]
  }
  
  
 # print(paste("Initial data has ", n.segs, " segments.", sep = ""))
#  print(paste("Using break-length threshold of ", brk.len.thresh, sep = ""))
#  print(paste("Combined data has ", dim(seq.dat)[1], " segments.", sep = ""))
  return(seq.dat)
}

# ---------------------------------- getLST --------------------------------------- #
getLST <- function(seq.dat)
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  # output:
  #   HRD.LST (numeric)
{
  # HRD-LST
  # large state transition
  seq.dat$brk.len <- 0
  seq.dat$LST <- FALSE
  
  seq.dat <- combineSeg(seq.dat, 180000)
  n.segs <- dim(seq.dat)[1]
  
  for( i in 1:(n.segs - 1))
  {
    seq.dat$brk.len[i + 1] <- seq.dat$start.pos[i + 1] - seq.dat$end.pos[i]
    seq.dat$LST[i] <- seq.dat$brk.len[i] < 3e06 &
                        seq.dat$cross.arm[i] == FALSE & 
                          seq.dat$s[i] > 10e06 & 
                            seq.dat$s[i + 1] > 10e06
  }
  
  write.table(seq.dat, "tmp-seq-dat.txt", col.names = TRUE, append = FALSE, row.names = FALSE)
  HRD.LST <- sum(seq.dat$LST)
  return(HRD.LST)
  
}
# ------------------------------------------------------------------------------- #

preprocessSeq <- function( seq.dat )
{
  # matches reference data to the correct chromosome in the subject data
  key = match(seq.dat$chromosome, ref.dat$chromosome)
  
  # ---
  # define allelic imbalance (AI), telomere positions, segment length, cross arms
  # these are used in multiple functions, attach to seq.dat here
  
  # segment length
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
  
  seq.dat$post.telomere <- (seq.dat$start.pos - ref.dat$p.telomere.end[key] <= 1000)
  seq.dat$pre.telomere  <- (ref.dat$q.telomere.start[key] - seq.dat$end.pos <= 1000)
  # ---
  
  return(seq.dat)
}


# ------------------------------------- hrd.stats ------------------------------------- #
# hrd.stats is a function to compute the three HRD metrics (HRD-LOH, HRD-NtAI, and
# HRD-LST), as well as total HRD and mean HRD.
# these calculations are based on kara's work
#
# input: seq.dat, a data.frame with chromosome, start.pos, end.pos, CNt, alleleA, alleleB
# output: out, a data.frame with HRD metrics
hrd.stats <- function(seq.dat, ploidy.dat, CN.dat)
{
  
  seq.dat <- preprocessSeq(seq.dat)
  
  # raw data
  HRD.NtAIr <- getNtAI(seq.dat)
  HRD.TAIr  <- getTAI( seq.dat )
  HRD.LSTr  <- getLST( seq.dat )
  HRD.LOH   <- getLOH( seq.dat )
  
  # create an index of main.CN segments; these get removed to normalize TAI
  rm.ind <- c()
  
  for( i in 1:length(CN.dat$chromosome))
  {
    rm.ind <- c(rm.ind, which(seq.dat$chromosome == CN.dat$chromosome[i] & seq.dat$CNt == CN.dat$main.CN[i]))
  }
  
  # normalized
  HRD.NtAIm <- getNtAI(seq.dat[-rm.ind,])
  HRD.TAIm <- getTAI(seq.dat[-rm.ind,])
  
  # any large genomic rearrangement could be increased simply by having more 
  # chromosomes (higher ploidy), and not because of the biological process
  # adjust for this
  HRD.LSTm = HRD.LSTr - (15.5 * ploidy.dat$ploidy.estimate[2])
  
  
  
  out = data.frame(HRD.LOH   = HRD.LOH, 
                   HRD.TAIr  = HRD.TAIr,
                   HRD.TAIRm = HRD.TAIm,
                   HRD.NtAIr = HRD.NtAIr,
                   HRD.NtAIm = HRD.NtAIm,
                   HRD.LSTr  = HRD.LSTr,
                   HRD.LSTm  = HRD.LSTm )
  
  
  return(out)
  
}
# ------------------------------------------------------------------------------- #



# ---------------------------------- getHRD.Data ---------------------------------- #
getHRD.Data <- function( sub.id, seq.dat, ploidy.dat )
  # setup data structures for HRD analysis
  # input: sub.id (character)
  #        seq.dat, a data.frame with columns: 
  #           chromosome, start.pos, end.pos, CNt, alleleA, alleleB
  #        ploidy.dat, a data.frame with columns:
  #           cellularity, ploidy.estimate, ploidy.mean.cn
  # output: hrd.dat, a data.frame with all of the calculated HRD scores
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
  
  levels(seq.dat$chromosome) <- levels(ref.dat$chromosome)
  
  
  seq.dat$frac.chr <-  (seq.dat$end.pos - seq.dat$start.pos) / ref.dat$chr.size[match(seq.dat$chromosome, ref.dat$chromosome)]
  
  
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
  
  # compute HRD stats, and record along with id
  
  hrd.dat = round(hrd.stats(seq.dat, ploidy.dat, CN.dat),3)
  hrd.dat$ID <- as.character(sub.id)
  
  return(hrd.dat)
}
# ------------------------------------------------------------------------------------- #


# get number of segments per break-length threshold
mb.thresh <- seq(from = 10000, to = 3e06, by = 10000)
out <- rep(0, length(mb.thresh) )
ct <- 0
for( i in mb.thresh)
{
  ct <- ct + 1
  x <- combineSeg(seq.dat, i)
  out[ct] <- dim(x)[1]
  
}

png("n.segs by break length.png")
plot(mb.thresh, out, type = "l")
dev.off()


# ------------------------------------------------------------------------------------- #
# ----------------------------------- end functions ----------------------------------- #
# ------------------------------------------------------------------------------------- #
