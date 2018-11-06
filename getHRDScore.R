

# BRADtools update, pluta 10/6/18
# v1.0: john pluta & kara maxwell 11/7/2016


# example usage:
# sub.id = "subX"
# ploidy.dat <- read.table("653-Brca2Br32-T.seqz_confints_CP.txt", header = T)
# seq.dat <- read.table("653-Brca2Br32-T.seqz_segments.txt",  header = T)

# hrd.dat = getHRD.Data( sub.id, seq.dat, ploidy.dat )
options(error=dump.frames)



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

# ------------------------------------- hrd.stats ------------------------------------- #
# hrd.stats is a function to compute the three HRD metrics (HRD-LOH, HRD-NtAI, and
# HRD-LST), as well as total HRD and mean HRD.
# these calculations are based on kara's work
#
# input: seq.dat, a data.frame with chromosome, start.pos, end.pos, CNt, alleleA, alleleB
# output: out, a data.frame with HRD metrics
hrd.stats <- function(seq.dat, ploidy.dat, CN.dat)
{
  # matches reference data to the correct chromosome in the subject data
  key = match(seq.dat$chromosome, ref.dat$chromosome)
  
  # check the sequencing data input to make sure it has the data we need
  if( any(!(seq.cols.needed %in% colnames(seq.dat))))
  {
    print(paste("column", 
                seq.cols.needed[which(!is.true(seq.cols.needed %in% colnames(seq.dat)))],
                "missing in the seq data.", sep=" "))
    print(paste("Looking for columns:", seq.cols.needed, sep=" "))
    print(paste("Found:", seq.cols.needed[seq.cols.needed %in% colnames(seq.dat)]))
    stop("Column mismatch. Exiting.")
  } 
  
  seq.dat$s <- seq.dat$end.pos - seq.dat$start.pos
  seq.dat$chr.size <- ref.dat$chr.size[match(seq.dat$chromosome, ref.dat$chromosome)]
  
  # HRD-LOH calculation
  # loss of heterozygosity
  # if B == 0 & s > 15000mbp & within chromosome, then HRD-LOH is TRUE
  seq.dat$LOH <- seq.dat$B == 0
  seq.dat$sub.chr <- (seq.dat$s / seq.dat$chr.size) < .9
  seq.dat$HRD.LOH <- (seq.dat$s > 15000000) & seq.dat$sub.chr & seq.dat$LOH
  
  HRD.LOH = sum(seq.dat$HRD.LOH)
  
  # HRD-NtAI calculation
  # non-transcriptome allelic imbalance
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
  
  seq.dat$HRD.NtAIr <- (seq.dat$s > 11000000 & seq.dat$AI & !seq.dat$cross.arm & !seq.dat$post.telomere & !seq.dat$pre.telomere)
  HRD.NtAIr = sum(seq.dat$HRD.NtAIr)
  
  # create an index of main.CN segments; these get removed
  rm.ind <- c()
  
  for( i in 1:length(CN.dat$chromosome))
  {
    rm.ind <- c(rm.ind, which(seq.dat$chromosome == CN.dat$chromosome[i] & seq.dat$CNt == CN.dat$main.CN[i]))
  }
  
  HRD.NtAIm <- getNtAIm(seq.dat[-rm.ind,])
  
  
  # HRD-LST
  seq.dat$brk.on.arm <- paste(seq.dat$chromosome, seq.dat$start.arm, seq.dat$end.arm, sep="")
  n.segs <- length(seq.dat$brk.on.arm[seq.dat$s > 3000000])
  n.breaks <- n.segs - length(unique(seq.dat$brk.on.arm)) + 1
  HRD.LSTm = n.breaks - (15.5 * ploidy.dat$ploidy.estimate[2])
  HRD.LSTr = n.breaks
  
  out = data.frame(HRD.LOH=HRD.LOH, 
                   HRD.NtAIr=HRD.NtAIr,
                   HRD.NtAIm=HRD.NtAIm,
                   HRD.LSTm=HRD.LSTm, 
                   HRD.LSTr=HRD.LSTr,
                   HRD.TOTALrr=sum(HRD.LOH, HRD.NtAIr, HRD.LSTr), 
                   HRD.MEANrr=mean(c(HRD.LOH, HRD.NtAIr, HRD.LSTr)),
                   HRD.TOTALrm=sum(HRD.LOH, HRD.NtAIm, HRD.LSTr),
                   HRD.MEANrm =mean(c(HRD.LOH, HRD.NtAIm, HRD.LSTr)),
                   HRD.TOTALmr=sum(HRD.LOH, HRD.NtAIr, HRD.LSTm), 
                   HRD.MEANmr=mean(c(HRD.LOH, HRD.NtAIr, HRD.LSTm)),
                   HRD.TOTALmm=sum(HRD.LOH, HRD.NtAIm, HRD.LSTm),
                   HRD.MEANmm =mean(c(HRD.LOH, HRD.NtAIm, HRD.LSTm)))
  
  
  
  return(out)
  
}
# ------------------------------------------------------------------------------- #


# ------------------------------- getNtAIm --------------------------------------- #
# function to get NtAI without including main CNt segments
# non-telomeric allelic imbalance
getNtAIm <- function(seq.dat)
  # input:
  #   seq.dat (data.frame), the sequencing data (eg, .seqz_segments.txt)
  # output:
  #   HRD.NtAIm (numeric)
{
  seq.dat$HRD.NtAIm <- (seq.dat$s > 11000000 & seq.dat$AI & !seq.dat$cross.arm & !seq.dat$post.telomere & !seq.dat$pre.telomere)
  HRD.NtAIm <- sum(seq.dat$HRD.NtAIm)
  return(HRD.NtAIm)
}
# --------------------------------------------------------------------------------- #


# ---------------------------------- getHRD.Data ---------------------------------- #
getHRD.Data <- function( sub.id, seq.dat, ploidy.dat )
  # script to setup data structures for HRD analysis
  # input: sub.id (character)
  #        seq.dat, a data.frame with columns: 
  #           chromosome, start.pos, end.pos, CNt, alleleA, alleleB
  #        ploidy.dat, a data.frame with columns:
  #           cellularity, ploidy.estimate, ploidy.mean.cn
  # output: hrd.dat, a data.frame with all of the calculated HRD scores
{
  
  levels(seq.dat$chromosome) <- levels(ref.dat$chromosome)
  hrd.dat <- data.frame(ID=sub.id,
                        HRD.LOH     = 0, 
                        HRD.NtAIr   = 0,
                        HRD.NtAIm   = 0,
                        HRD.LSTm    = 0,
                        HRD.LSTr    = 0,
                        HRD.TOTALrr = 0,
                        HRD.MEANrr  = 0,
                        HRD.TOTALrm = 0,
                        HRD.MEANrm  = 0,
                        HRD.TOTALmr = 0, 
                        HRD.MEANmr  = 0,
                        HRD.TOTALmm = 0,
                        HRD.MEANmm  = 0)
  
  hrd.dat$ID <- as.character( hrd.dat$ID )
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
  m <- dim(hrd.dat)[2]
  hrd.dat[,2:m] = round(hrd.stats(seq.dat, ploidy.dat, CN.dat),3)
  
  return(hrd.dat)
}
# ------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------- #
# ----------------------------------- end functions ----------------------------------- #
# ------------------------------------------------------------------------------------- #
