#!/usr/bin/env Rscript

# john pluta & kara maxwell 11/7/2016

rm(list=ls())

options(error=dump.frames)


args = commandArgs(trailingOnly = TRUE)
if(length(args) < 3)
{
	print("need to provide 3 arguments:")
	print("1: group.dir = group-level directory containing individual directories for each subject's sequenza data")
	print("2: out.dir   = output directory")
	print("3: gene.file  = a text file describing the genes of intetest, in the format GeneID chr hg19_start hg19_end")
	stop()
}


# group-level directory containing individual directories for each subject's sequenza data
group.dir = args[1]

# output directory
out.dir   = args[2]

# genes of interest
# this is a text file with columns:
# GeneID CHR hg19_start hg19_end
# the last two columns are the start and end physical coordinates of the gene
#
# example:
# GeneID  chr hg19_start  hg19_end
# AKT1  14  105235687 105262080
# AKT2  19  40736224  40791302
gene.file  = args[3]




# ===================================================================================== #
# ================================ constants ========================================== #
# ===================================================================================== #
# the columns needed for an analysis, per standard sequenza naming conventions
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





# ===================================================================================== #
# ================================== functions ======================================== #
# ===================================================================================== #
# ------------------------------------- hrd.stats ------------------------------------- #
# hrd.stats is a function to compute the three HRD metrics (HRD-LOH, HRD-NtAI, and
# HRD-LST), as well as total HRD and mean HRD.
# these calculations are based on kara's work
#
# input: 
# output:
hrd.stats <- function(seq.dat, ploidy.dat, CN.dat)
{
  # matches reference data to the correct chromosome in the subject data
  key = match(seq.dat$chromosome, ref.dat$chromosome)
  
  # check the sequenza input to make sure it has the data we need
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
  # if B == 0 & s > 15000mbp & within chromosome, then HRD-LOH is TRUE
  seq.dat$LOH <- seq.dat$B == 0
  seq.dat$sub.chr <- (seq.dat$s / seq.dat$chr.size) < .9
  seq.dat$HRD.LOH <- (seq.dat$s > 15000000) & seq.dat$sub.chr & seq.dat$LOH
  
  HRD.LOH = sum(seq.dat$HRD.LOH)
  
  # HRD-NtAI calculation
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
getNtAIm <- function(seq.dat)
# input:
#   seq.dat (data.frame), the sequencing data (.seqz_segments.txt)
# output:
#   HRD.NtAIm (numeric)
{
  seq.dat$HRD.NtAIm <- (seq.dat$s > 11000000 & seq.dat$AI & !seq.dat$cross.arm & !seq.dat$post.telomere & !seq.dat$pre.telomere)
  HRD.NtAIm <- sum(seq.dat$HRD.NtAIm)
  return(HRD.NtAIm)
}
# ------------------------------------------------------------------------------- #


# -------------------------------- getPTENcopyStatus ---------------------------------- #
getPTENcopyStatus <- function( group.dir )
# input: 
#   group.dir (string), the group level directory
# output:
#   grp.dat (data.frame), alleles and copy number fo reach subejct
{
  # list the subdirectories, eg the individual subject directories, within
  # the overall group directory
  sub.dirs = list.dirs( group.dir )
  
  sub.dirs = sub.dirs[2:length(sub.dirs)]
  n <- length(sub.dirs)
  grp.hrd = rep(0, n)
  
  sub.ids <- rep("0", n)
  sub.out <- rep("0", n)
 
  print("Gathering PTEN Copy data...")
  pb <- txtProgressBar(min = 0, max=n, style=3 )
  
  for(i in 1:n)
  {
    sub.dir = sub.dirs[i]
    sub.id = strsplit(sub.dir, "/")[[1]][lengths(strsplit(sub.dir, "/"))]
    seq.file = paste(sub.dir, "/", sub.id, ".seqz_segments.txt", sep="")
    seq.dat = read.table(seq.file, header=TRUE)
    
    # per kara's specifications...
    dat <- seq.dat[seq.dat$chromosome == "chr10" &
    seq.dat$start.pos < 87917521 &
    seq.dat$end.pos   > 87917521,]
    
    # if the data is missing
    if( dim(dat)[1] == 0)
    {
      sub.out[i] <- "NA"
    } else
    {
      sub.out[i] <- paste(dat$CNt, dat$A, dat$B, sep="_")
      sub.ids[i] <- sub.id
    }
 
    setTxtProgressBar(pb, i)
  }
  print("Done!")
  grp.dat <- data.frame(sub.id=sub.ids, GCS=sub.out)
  close(pb)
  return(grp.dat)
}
# ------------------------------------------------------------------------------- #


# -------------------------------- getHRDData ---------------------------------- #
# function to compute various HRD measures from sequenza data
getHRDData <- function( group.dir )
# input: 
#   group.dir (string), the directory containing the subject-level subdirectories of sequenza data
# output:
#   hrd.dat (data.frame), data.frame of HRD data
{
  # list the subdirectories, eg the individual subject directories, within
  # the overall group directory
  sub.dirs = list.dirs( group.dir )

  # sub.dirs[1] = root dir
  # dont include this
  sub.dirs = sub.dirs[2:length(sub.dirs)]
  
  n = length(sub.dirs)
  sub.ids <- rep("0", n)
  
  # data frame to store the stats
  hrd.dat = data.frame(ID=sub.ids,
                       HRD.LOH=rep(0,n), 
                   HRD.NtAIr=rep(0,n),
                   HRD.NtAIm=rep(0,n),
                   HRD.LSTm=rep(0,n),
                   HRD.LSTr=rep(0,n),
                   HRD.TOTALrr=rep(0,n),
                   HRD.MEANrr=rep(0,n),
                   HRD.TOTALrm=rep(0,n),
                   HRD.MEANrm =rep(0,n),
                   HRD.TOTALmr=rep(0,n), 
                   HRD.MEANmr=rep(0,n),
                   HRD.TOTALmm=rep(0,n),
                   HRD.MEANmm =rep(0,n))
  
  hrd.dat$ID <- as.character(hrd.dat$ID)
  
  print("Gathering HRD statistics...")
 
  pb <- txtProgressBar(min=0, max=n, style=3)
  
  for(i in 1:length(sub.dirs))
  {
    sub.dir = sub.dirs[i]
    sub.id = strsplit(sub.dir, "/")[[1]][lengths(strsplit(sub.dir, "/"))]
    hrd.dat$ID[i] <- sub.id
    
    
    # read in seq and ploidy data
    seq.file = paste(sub.dir, "/", sub.id, ".seqz_segments.txt", sep="")
    ploidy.file = paste(sub.dir, "/", sub.id, ".seqz_confints_CP.txt", sep="")
    seq.dat = read.table(seq.file, header=TRUE)
    ploidy.dat = read.table(ploidy.file, header=TRUE)
    levels(seq.dat$chromosome) <- levels(ref.dat$chromosome)
    
    
    # calculate main CN
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
    hrd.dat[i,2:m] = round(hrd.stats(seq.dat, ploidy.dat, CN.dat),3)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  print("Done!")
  return(hrd.dat)
}
# ------------------------------------------------------------------------------- #


# -------------------------------- getPloidyData -------------------------------- #
getPloidyData <- function( group.dir )
# input: 
#   group.dir (string), the directory containing the subject-level subdirectories of sequenza data
# output:
#   ploidy.stats.dat (data.frame), ploidy data
{
  # list the subdirectories, eg the individual subject directories, within
  # the overall group directory
  sub.dirs = list.dirs( group.dir )
  
  # sub.dirs[1] = "/Users/jpluta/Desktop/sequenza"
  # dont include this
  sub.dirs = sub.dirs[2:length(sub.dirs)]
  
  n = length(sub.dirs)
  sub.ids <- rep("0", n)
  
  # data frame to store the stats
  ploidy.stats.dat = data.frame(ID=sub.ids,
                       cellularity1=rep(0,n), 
                       cellularity2=rep(0,n),
                       cellularity3=rep(0,n),
                       estimate1=rep(0,n),
                       estimate2=rep(0,n),
                       estimate3=rep(0,n),
                       mean.cn=rep(0,n))
                       
  
  ploidy.stats.dat$ID <- as.character(ploidy.stats.dat$ID)
  
  print("Gathering ploidy data...")
  
  pb <- txtProgressBar(min=0, max=n, style=3)
  
  for(i in 1:length(sub.dirs))
  {
    sub.dir = sub.dirs[i]
    sub.id = strsplit(sub.dir, "/")[[1]][lengths(strsplit(sub.dir, "/"))]
    ploidy.stats.dat$ID[i] <- sub.id
    
    # read ploidy data
    
    ploidy.file = paste(sub.dir, "/", sub.id, ".seqz_confints_CP.txt", sep="")
    ploidy.dat = read.table(ploidy.file, header=TRUE)
   
    ploidy.stats.dat[i,2:8] <- c(ploidy.dat$cellularity, ploidy.dat$ploidy.estimate, ploidy.dat$ploidy.mean.cn[1])
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  print("Done!")
  return(ploidy.stats.dat)
}
# ------------------------------------------------------------------------------- #



# ---------------------------------getGeneData---------------------------------------- #
# generic function to get gene data, works a little differently than p10
getGeneData <- function( group.dir, chr, start, end )
# input: 
#   group.dir = the directory containing the individual subjects data
#   chr = the chromosome of interest
#   start = start position (in bp) of the chromosome
#   end = ending position (in bp) of the chromosome
# output:
#   grp.dat = group-level subject data, genomic copy status for each 
#     subject (copy number and allele)
{
  # list the subdirectories, eg the individual subject directories, within
  # the overall group directory
  sub.dirs = list.dirs( group.dir )

  # the first directory listed is the directory itself; skip this
  sub.dirs = sub.dirs[2:length(sub.dirs)]
  n <- length(sub.dirs)
  grp.hrd = rep(0, n)

  sub.ids <- rep("0", n)
  sub.out <- rep("0", n)


  pb <- txtProgressBar(min = 0, max=n, style=3 )

  # calculate midpoint of the gene
  midpoint <- round(  (end - start)/2,  0) + start

  
  for(i in 1:n)
  {
    sub.dir = sub.dirs[i]
    sub.id = strsplit(sub.dir, "/")[[1]][lengths(strsplit(sub.dir, "/"))]
    seq.file = paste(sub.dir, "/", sub.id, ".seqz_segments.txt", sep="")
    seq.dat = read.table(seq.file, header=TRUE)
  
    # adjust chr naming convention (from eg 7 to chr7)
    # need the info from seq.dat, but only want to run this once
    if(substr(seq.dat$chromosome[1], 1, 3) == "chr" & nchar(chr) == 1)
    {
        chr = paste("chr", chr, sep="")
    }
   
    # per kara's specifications...
    dat <- seq.dat[as.character(seq.dat$chromosome) == chr &
                   seq.dat$start.pos < midpoint &
                   seq.dat$end.pos   > midpoint,]
  
    # if the data is missing
    if( dim(dat)[1] == 0)
    {
      sub.out[i] <- "NA"
      
    } else
    {
      sub.out[i] <- paste(dat$CNt, dat$A, dat$B, sep="_")
      sub.ids[i] <- sub.id
    }
  
    setTxtProgressBar(pb, i)
  }
  
  if( sum(sub.out == "NA") == n )
  {

    print("WARNING: no data found!")
  }
  
  grp.dat <- data.frame(sub.id=sub.ids, GCS=sub.out)
  close(pb)
  return(grp.dat)
}
# ------------------------------------------------------------------------------- #
# ===================================================================================== #
# ===================================================================================== #
# ===================================================================================== #



# ===================================================================================== #
# ===================================== MAIN ========================================== #
# ===================================================================================== #


setwd(out.dir)


# r is raw, m is corrected.

# TOTAL = LOH + LST +NtAI
# TOTALrr = LOH + LSTraw + NtAIraw
# TOTALrm = LOH + LSTraw + NtAI modified
# we used mr in our analysis (LSTmodified, NtAI raw)

# input genes
gen.dat <- read.table(gene.file, header=TRUE,
                      colClasses= c("character", "character", "integer", "integer"), 
                      col.names = c("Gene", "Chr", "Start", "End"))

for( i in 1:dim(gen.dat)[1])
{
  print(paste("Writing data for gene: ", gen.dat$Gene[i], "...", sep=""))
  out.dat <- getGeneData( group.dir, gen.dat$Chr[i], gen.dat$Start[i], gen.dat$End[i])
  write.table(out.dat, paste(gen.dat$Gene[i], "-genomic-copy-status.txt", sep=""), quote=FALSE, append=FALSE, row.names=FALSE, col.names = TRUE)
  print("Done!")
}

grp.hrd = getHRDData(group.dir)

# p10 copy stats
PTEN.dat <- getPTENcopyStatus(group.dir)

# ploidy statistics
ploidy.stats <- getPloidyData(group.dir)

# output
print("Writing data ...")
write.table(PTEN.dat, "PTEN-genomic-copy-status.txt", quote=FALSE, append=FALSE, row.names=FALSE, col.names = TRUE)
write.table(grp.hrd, "Group-HRD.txt", quote=FALSE,append=FALSE,row.names=FALSE,col.names=TRUE)
write.table(ploidy.stats, "ploidy-stats.txt", quote=FALSE, append=FALSE, row.names=FALSE, col.names = TRUE)
print("Done")
# ===================================================================================== #
# ===================================================================================== #
# ===================================================================================== #