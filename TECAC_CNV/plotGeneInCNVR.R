
#  pluta  9/8/21
# v1.0

# map of gene and CNVR overlap

# example annotation data
# annot.dat <- read.table("call_gene/interval_annotation.txt", header = TRUE, sep = "\t")

library(ggplot2)

# ----------------------------- plotGenesInCNVR ---------------------------------- #
plotGenesInCNVR <- function(cnvr.id, cnvr.dat, annot.dat, probe.dat = NULL)
#  input: cnvr.id (string), id of the CNVR, e.g "CNVR_190"
#         cnvr.dat (data.frame), df of CNVR information (chr and start/end); the output of runCNVRAssoc (cnvr.txt)
#         annot.dat (data.frame), from the interval_annotation.txt file created by HandyCNV; the output of call_gene
#         probe.dat (data.frame), the bim file containing snp locations for the chip used to generate cnvs
#  output: p1  (ggplot object), a plot of the CNVR and overlapping genes

{
  if(!(cnvr.id %in% cnvr.dat$CNVR_ID))
  {
    stop(paste0("CNVRID ", cnvr.id, " not found in data."))
  }
  
  if(!(all( c("Chr", "Start", "End") %in% colnames(cnvr.dat))))
  {
    print("column name mismatch")
    stop("needs column Chr, Start, and End")
  }
  
  cnvr.dat <- data.frame( chr = cnvr.dat$Chr[ cnvr.dat$CNVR_ID == cnvr.id ],
                          start = cnvr.dat$Start[ cnvr.dat$CNVR_ID == cnvr.id ],
                          end = cnvr.dat$End[ cnvr.dat$CNVR_ID == cnvr.id],
                          id = cnvr.id)
  cnvr.size <- cnvr.dat$end - cnvr.dat$start
  
  plot.label <- paste0( "chr", cnvr.dat$chr, ":", cnvr.dat$start, "-", cnvr.dat$end)
  
  # add a check to make sure  chrs match
  gene.dat <- annot.dat[ which(annot.dat$ID == cnvr.dat$id), ]
  gene.dat <- gene.dat[ !duplicated(gene.dat$name2), ]
  
  if( cnvr.dat$chr != unique(gene.dat$Chr) )
  {
    print(paste0("cnvr.dat is on chr", cnvr.dat$chr, " but gene data is mapped to chr", gene.dat$Chr[1], "."))
    stop("are you using the correct CNVR_ID?")
  }
  
  
  df <- data.frame( start = gene.dat$Start, end = gene.dat$End, name = gene.dat$name2)
  df$position <- seq(1:length(df$start)) + 1
  
  write.table(df$name, paste0(cnvr.id, "_genes.txt"), col.names = FALSE,  quote  = FALSE, row.names = FALSE, append = FALSE)
  xmin <- min(c(cnvr.dat$start, df$start))  - (cnvr.size/4)
  xmax <- max(c(cnvr.dat$end, df$end)) + 10000
  
  p1 <- ggplot(data = df) + 
    geom_segment(aes(x = start, xend = end, y = position, yend = position), size = 4, color = "blue") +
    geom_segment(data = cnvr.dat, aes(x = start, xend = end, y = 1,  yend = 1), size = 4, color = "red") +
    geom_label(data = df, aes(x = start, y = position, label = name), size = 4, nudge_x = -nchar(df$name) * 0.02 * cnvr.size) +
    geom_label(data = cnvr.dat, aes(x  = start +  (end - start)/2, y = 1, label = id), nudge_y = max(df$position) * 0.025) +
    geom_vline(xintercept = cnvr.dat$start, linetype = "dashed", color = "red") +
    geom_vline(xintercept = cnvr.dat$end, linetype = "dashed", color = "red") +
    xlim( xmin, xmax) + 
    ggtitle(plot.label) +
    theme_minimal() 
  
  # annotate plot with probe location
  if( !is.null(probe.dat))
  {
    probe.dat <- probe.dat[  probe.dat$chr == cnvr.dat$chr &
                               probe.dat$bp >= xmin &
                               probe.dat$bp <=  xmax, ]
    probe.dat$pos   <- 1
    p1 <- p1 + geom_point(data = probe.dat, aes(bp, pos),
                          stat= "identity", shape = 108, size = 3)
  }
  
  return(p1)
}
# ------------------------------------------------------------------------------ #




# ============================= main =========================== #
args = commandArgs(trailingOnly = TRUE)

if( length(args) < 3)
{
  print("need to provide 3 arguments: ")
  print("")
  print(paste0("you provided ", length(args), " arguments:"))
  print(paste0("CNVRRESFILE = ", args[1]))
  print("the cnvr_res.txt file, produced by runCNVRAssoc.R")
  print("the interval_annotation.txt file, produced by runCNVRAssoc.R")
  print(paste0("ANNOTFILE = ", args[2]))
  print("the annotation fi")
  print(paste0("CNVRIDFILE = ", args[3]))
  print("simple text list of CNVR ids of interest, in the format CNVR_###")
  print("")
  stop()
}

CNVRRESFILE = args[1]
ANNOTFILE = args[2]
CNVRIDFILE = args[3]

# load the probe data
probe.dat <- read.table("IKN_TECAC_reclusterx14520.bim", header  = FALSE)
probe.dat <-  probe.dat[ probe.dat$V1 > 0 & probe.dat$V1 < 23,]
colnames(probe.dat) <- c("chr", "rsid", "gd", "bp", "a1", "a2")

# list of association results
cnvr.dat <- read.table(CNVRRESFILE, header = TRUE, sep = " ")
annot.dat  <- read.table(ANNOTFILE, header = TRUE, sep = "\t")
cnvr.ids <- read.table(CNVRIDFILE, header = FALSE)$V1

for( cnvr.id in cnvr.ids )
{
  p1 <- plotGenesInCNVR(cnvr.id, cnvr.dat, annot.dat, probe.dat)
  OUTNAME <- paste0(cnvr.id, "_plot.png")
  png(OUTNAME)
  print(p1)
  dev.off()
}
