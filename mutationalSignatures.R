
# john pluta 8/7/18
# jpluta@pennmedicine.upenn.edu

# ------------------------------------------------------------------ #
# ---------------------- preprocessor ------------------------------ #
rm(list=ls())

# uses bioconductor
source("http://bioconductor.org/biocLite.R")

if("ggdendro" %in% rownames(installed.packages()) == FALSE)
{
  install.packages("ggdendro")
}

if("SomaticSignatures" %in% rownames(installed.packages()) == FALSE)
{
  biocLite("SomaticSignatures")
}

if("deconstructSigs" %in% rownames(installed.packages()) == FALSE)
{
  biocLite("deconstructSigs")
}

if("SomaticCancerAlterations" %in% rownames(installed.packages()) == FALSE)
{
  biocLite("SomaticCancerAlterations")
}

if("BSgenome.Hsapiens.UCSC.hg38" %in% rownames(installed.packages()) == FALSE)
{
  biocLite("BSgenome.Hsapiens.UCSC.hg38")
}

if("BSgenome.Hsapiens.UCSC.hg19" %in% rownames(installed.packages()) == FALSE)
{
  biocLite("BSgenome.Hsapiens.UCSC.hg19")
}

if("rtracklayer" %in% rownames(installed.packages()) == FALSE)
{
  biocLite("rtracklayer")
}

if("ggplot2" %in% rownames(installed.packages()) == FALSE)
{
  install.packages("ggplot2")
}

# load libraries
library(SomaticSignatures)
library(SomaticCancerAlterations)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(deconstructSigs)
library(rtracklayer)
library(ggdendro)
# ------------------------------------------------------------------ #
# ------------------------------------------------------------------ #






# ------------------------------------------------------------------ #
# ---------------------------- functions --------------------------- #
# ------------------------------------------------------------------ #


# -------------------------- convertAFtoRLE ------------------------------- #
convertAFtoRLE <- function( vr )
# if certain fields in a vrange object are NOT in RLE format, convert them
# this function works on the AF (allele frequency) field
  # input: vr, a VRange object
  # output: vr, a VRange o bject
{
  l <- rep(vr@elementMetadata@listData$AF[1], length(vr@elementMetadata@listData$AF))
  vr@elementMetadata@listData$AF <- NumericList(split(l, seq(1:length(l))), compress=TRUE)
  names(vr@elementMetadata@listData$AF@partitioning) <- NULL
  return(vr)
}
# ------------------------------------------------------------------ #


# -------------------------- convertACtoRLE ------------------------------- #
convertACtoRLE <- function( vr )
# if certain fields in a vrange object are NOT in RLE format, convert them
# this function works on the AC (allele count) field
  # input: vr, a VRange object
  # output: vr, a VRange o bject
{
  l <- rep(vr@elementMetadata@listData$AC[1], length(vr@elementMetadata@listData$AC))
  vr@elementMetadata@listData$AC <- IntegerList(split(l, seq(1:length(l))), compress=TRUE)
  names(vr@elementMetadata@listData$AC@partitioning) <- NULL
  return(vr)
}
# ------------------------------------------------------------------ #



# ---------------------- createSomSigPlots_NMF --------------------- #
# get mutational signatures via NMF decomposition
# create all plots
#
# input: grp_mm, a group level motif matrix
#        name (string) that is prepnded to plots
#        n_sigs (integer), number of signatures to use. determined by assessNumberSignatures
#         or the getNSignatures function
# 
# output: a series of 7 plots in a pdf; text output of the data used in the plots
createSomSigPlots <- function(grp_motifs, name, n_sigs)
{
  
  grp_mm <- SomaticSignatures::motifMatrix(do.call(c, grp_motifs), group="grp", normalize=TRUE)
  
  if( any(grp_mm == 0))
  {
    stop("grp_mm cannot have 0 valued or NA entries")
  }
  
  grp.p0 <- SomaticSignatures::plotMutationSpectrum(do.call(c, grp_motifs), group="grp") + ggtitle("Mutation Spectrum by Group")
  
  
  # create the MutationalSignatures object
  print("creating mutational sigs object...")
  grp_nmf = SomaticSignatures::identifySignatures(grp_mm, n_sigs, nmfDecomposition)
  print("done")
  
  # standard plot
  print("Plotting signatures...")
  grp.p1 <- SomaticSignatures::plotSignatures(grp_nmf) + ggtitle(paste(name,"Barchart",sep="-"))

  
  # heatmap
  grp.p2 <- SomaticSignatures::plotSignatureMap(grp_nmf) + ggtitle(paste(name,"Heatmap",sep="-"))
  
  # observed specturm
  grp.p3 <- SomaticSignatures::plotObservedSpectrum(grp_nmf) + ggtitle(paste(name,"Observed Spectrum",sep="-"))
  
  # fitted spectrum
  grp.p4 <- SomaticSignatures::plotFittedSpectrum(grp_nmf) + ggtitle(paste(name,"Fitted Spectrum",sep="-"))
  
  # sample map
  grp.p5 <- SomaticSignatures::plotSampleMap(grp_nmf) + ggtitle(paste(name,"Sample Map",sep="")) 
  
  # samples
  grp.p6 <- SomaticSignatures::plotSamples(grp_nmf) + ggtitle(paste(name,"Samples",sep=""))
  
  # write data thats used in the plots, for further analysis
  write.table(grp_nmf@samples, paste(name,"sigdat.txt", sep="-"), row.names = TRUE, col.names=TRUE, append=FALSE,
              quote=FALSE)
  
  
  # write everything to pdf
  pdf(paste(name, "plots.pdf", sep="_"))
  print(grp.p0)
  print(grp.p1)
  print(grp.p2)
  print(grp.p3)
  print(grp.p4)
  print(grp.p5)
  print(grp.p6)
  dev.off()
  
}
# ------------------------------------------------------------------ #



# -------------- convert_to_deconstructSigs_format ----------------- #
convert_to_deconstructSigs_format <- function( dat, gref )
  # input: dat, a VRange object
  #        gref (BSGenome object), reference genome
  # output: a deconstructSigs format object
{
  
  # get values from VRanges object, convert values as needed
  # expand RLE objects where needed
  len <- dat@seqnames@lengths
  val <- as.character(dat@seqnames@values)
  chr <- rep(val, len)
  samp <- rep(as.character(dat@sampleNames@values), dat@sampleNames@lengths)
  
  # build data.frame- input names are specifically what deconstructSig takes, dont
  # change these.
  dat <- data.frame(Sample=as.character(samp),
                    chr=as.character(chr),
                    pos=as.integer(dat@ranges@start),
                    ref=dat@ref,
                    alt=dat@alt)
  
  # do the actual conversion and return
  sigs.input <- deconstructSigs::mut.to.sigs.input(mut.ref = dat, sample.id="Sample", chr="chr",
                                                   pos="pos", ref="ref", alt="alt", bsg=gref)
  return(sigs.input)
}
# ------------------------------------------------------------------ #


# --------------------- run_deconstructSigs_grp -------------------- #
runDeconstructSigs <- function( grp_motifs, refBuild )
  # input: grp_motifs, group of vcf files
  #        refBuild (string), a stringscribing the reference build (one of hg19 or hg38)
  # output: group level statistics- weights per signature per subject,
  #   count of number of times each signature appears, and the error
  #   for each subject. also outputs plots for each subject
  #   does not return a value, but rather prints output to file
{
  #
  # convert group VRange object to deconstructSigs format
  # set reference genome
  if( refBuild == "hg19" )
  {
    gref <- BSgenome.Hsapiens.UCSC.hg19
  } else
    if( refBuild == "hg38" )
    {
      gref <- BSgenome.Hsapiens.UCSC.hg38
    } else
    {
      print(paste(refBuild, " is not a valid reference genome.", sep =""))
      stop("Select from hg19 or hg38")
    }
  
  # alexandrov data
  data(signatures21, package="SomaticSignatures")
  
  dat <- convert_to_deconstructSigs_format(grp_motifs, gref)
  
  err <- rep(0, dim(dat)[1])      # the error term for each individual; the difference between observed and reference signature
  sig.ct <- rep(0, 27)           # a count of which signatures had non-zero weight. totaled over all subjects.
  w <- c()                       # matrix to contain the weights assigned to each signature, for each subject
   
  
  print("Running deconstruct sigs...")
  pb <- txtProgressBar( min = 0, max =dim(dat)[1], style = 3)
  
  # TODO: put this in lapply
  # iterate through each subject
  for( i in 1:dim(dat)[1])
  {
    id = row.names(dat)[i]    # subject id
    
    
    
    
    sub.dat = deconstructSigs::whichSignatures(tumor.ref = dat,
                                               signatures.ref = signatures.nature2013,
                                               sample.id = id, 
                                               contexts.needed = TRUE,
                                               tri.counts.method = "exome")
    
    
    # iteratively add the weights
    w <- rbind(w, as.numeric(sub.dat$weights))
    
    # write plots to pdf
    pdf(paste(id, "plots.pdf", sep="_"))
    deconstructSigs::plotSignatures(sub.dat, sub=id)
    makePie(sub.dat, sub=id, add.color=NULL)
    dev.off()
    
    # error
    # the difference between observed signature and composite of alexandrov signatures
    err[i] = sqrt( sum( sub.dat$diff^2))
    
    # product is the weight of each signature multiplied by the tumor weights
    # weights > 0 are the signatures present in the tumor. 
    sig.ct = sig.ct + as.integer(sub.dat$weights > 0)
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  print("done!")
  
  if( is.null( grp_motifs$grp) )
  {
    grp_motifs$grp <- "grp"
  }
  
  print("writing count file")
  
  # write signature counts to file
  sig.ct = cbind(colnames(sub.dat$weights), sig.ct)
  colnames(sig.ct) = c("Signature", "Count")
  write.table(sig.ct, paste(as.character(grp_motifs$grp[1]), "sig_ct.txt", sep="_"), append=FALSE, quote=FALSE,
              row.names = FALSE, col.names=TRUE)
  print("done")
  print("writing weights")
  # write the group level matrix of weights
  colnames(w) = colnames(sub.dat$weights)
  row.names(w) = row.names(dat)
 
  
  write.table(w, paste(grp_motifs$grp[1], "weights.txt", sep="_"), append=FALSE, quote=FALSE,
              row.names=TRUE, col.names=TRUE)
  print("done")
  print("writing errors")
  
  
  # write all the errors
  err = as.matrix(err)
  row.names(err) = row.names(dat)
  colnames(err) = c("error")
  write.table(err, paste(grp_motifs$grp[1], "error.txt", sep="_"), append=FALSE, quote=FALSE,
              row.names=TRUE, col.names=TRUE)
  
  print("done!")
}
# ------------------------------------------------------------------ #





# --------------------------- checkSeqLevels ----------------------- #
# check the names of the sequences in a VR object, and make sure
# chromosomes start with "chr" to match the reference genome
#
# input: seqlevels, a character vector with the names of sequence levels
#   from the VR object
# output: seqlevels, renamed to prepend "chr"
checkSeqLevels <- function( seqlevels )
{
  if( nchar(seqlevels[1]) == 1)
  {
    seqlevels <- paste("chr", seqlevels, sep = "")
  }
  
  return(seqlevels)
}
# ------------------------------------------------------------------ #


# -------------------------- readVCF.qc ----------------------------- #
# read a VCF file into a VRange object and do some basic qc. add mutation
#   alteration and context.
# 
# input: fname, the name of the vcf file to read in
#        gref, the reference genome object
# output: vr_A, the qc'd VRange object with mutation alteration and context added
readVCF.qc <- function( fname, gref )
{
  # ignore the "duplicate keys in header" error
  dat <- suppressWarnings(  readVcfAsVRanges( fname, seqinfo(scanVcfHeader(fname)) )  )
  
  # truncate to viable sequences
  dat <- keepSeqlevels(dat, seqlevels(dat)[1:24])
  genome(dat) <- gref@provider_version
  idx <- ref(dat) %in% DNA_BASES & alt(dat) %in% DNA_BASES
  dat <- dat[idx]
  seqlevels(dat) <- checkSeqLevels(seqlevels(dat))
  seqlengths(dat) <- seqlengths(gref)[seqlevels(gref) %in% seqlevels(dat)]
  
  
  # this adds the columns 'alteration' and 'context' to the vcf
  vr_A <- SomaticSignatures::mutationContext(dat, gref)

  # account for various encoding in the VCF files
  if( class(vr_A@elementMetadata@listData$AF) == "numeric" )
  {
    vr_A <- convertAFtoRLE(vr_A)
    
  }
  
  if( class(vr_A@elementMetadata@listData$AC) == "integer" )
  {
    vr_A <- convertACtoRLE(vr_A)
  }
  
  return(vr_A)
}
# ------------------------------------------------------------------ #


# ------------------------------ getVCFmotif --------------------------------- #
# create motif matrix for a given group
# input:
#   VcfFileDir, the directory contains the vcf files
#   refBuild, a string describing the genome reference build (one of either hg19 or hg38)
# output:
#   a vr_A, a VRange object containing all of the vcf data
getVCFmotif <- function( VcfFileDir, refBuild)
{
  
  
  # set reference genome
  if( refBuild == "hg19" )
  {
    gref <- BSgenome.Hsapiens.UCSC.hg19
  } else
    if( refBuild == "hg38" )
    {
      gref <- BSgenome.Hsapiens.UCSC.hg38
    } else
      {
        print(paste(refBuild, " is not a valid reference genome.", sep =""))
        stop("Select from hg19 or hg38")
      }
  
  
  
  print("Reading VCF files...")
  
  fnames <- paste(VcfFileDir, list.files(VcfFileDir, pattern="*.vcf$"), sep = "/")
  
  if( any(is.na(fnames)) | length(fnames) == 0)
  {
    stop(print(paste("ERROR: No VCF files found in directory:", VcfFileDir, sep = " ")))
  } 
  
  # do some basic qc
  vr_All <- lapply(as.list(fnames), readVCF.qc, gref = gref)
 
  
  # convert to group-level VRange object (from VRange list)
  vr_A <- do.call(c, vr_All)
  
  vr_A$grp <- unlist(strsplit(VcfFileDir, "/"))[length(unlist(strsplit(VcfFileDir, "/")))]
  sampleNames(vr_A) = droplevels(sampleNames(vr_A))
  
  print("done!")
  print(paste(length(vr_A@sampleNames@values), "subjects loaded.", sep=" "))
  return(vr_A)
}
# ------------------------------------------------------------------ #



getNSignatures <- function( grp_mm, maxNSig = ncol(grp_mm), nReplicates = 5)
{
  gof_nmf <- SomaticSignatures::assessNumberSignatures(grp_mm, seq(from=1, to=2, by=1), nReplicates = 5)
  grp.p2 <- SomaticSignatures::plotNumberSignatures(gof_nmf)
  png("nSigs.png")
  print(grp.p2)
  dev.off()
}


# ------------------------------------------------------------------ #
# ------------------------------------------------------------------ #
# ------------------------------ end functions --------------------- #
# ------------------------------------------------------------------ #









# ------------------------------------------------------------------ #
# ------------------------------ MAIN ------------------------------ #
# ------------------------------------------------------------------ #

# --- user defined --- #
outdir <- "/Users/jpluta/Desktop/forBrad/vcf/"
refBuild <- "hg19"

# program assumes vcf files are grouped in separate directories within the main vcf directory
VcfFileDir <- "/Users/jpluta/Desktop/forBrad/vcf/"
group.names <- c("grp1", "grp2")

# you can use the getNSignatures function to determine this
n.sigs <- 2



groups <- paste(VcfFileDir, group.names, sep = "/")
grp_motifs <- lapply( groups, getVCFmotif, refBuild)

for( i in 1:length(grp_motifs))
{
  grp_motifs[[i]]$grp <- group.names[i]
}



# run deconstruct sigs
lapply( grp_motifs, runDeconstructSigs, refBuild)

createSomSigPlots( grp_motifs, outdir, n.sigs)
# ------------------------------------------------------------------ #
# ------------------------------------------------------------------ #
# ------------------------------------------------------------------ #
