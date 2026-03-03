

#!/usr/bin/env Rscript

# pluta 3/3/26

# ------------------------ install required packages ------------------------- #
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

install_if_missing <- function(pkgs) {
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_install) > 0) {
    BiocManager::install(to_install, ask = FALSE, update = FALSE)
  }
}

bioc_pkgs <- c(
  "GenomicRanges",
  "IRanges"
)

install_if_missing(bioc_pkgs)
# ---------------------------------------------------------------------------- #



library(optparse)
library(GenomicRanges)
library(IRanges)


# -------------------------------- parse input ------------------------------- #
option_list <- list(
  make_option(c("-s", "--SNPLIST"),
              type = "character",
              help = "list of SNPs in chr:pos format")
)

parser <- OptionParser(
  usage = "%prog [options]",
  option_list = option_list
)

args <- parse_args(parser)

if(length(commandArgs(trailingOnly = TRUE)) == 0)
{
  print_help(parser)
  quit(status = 1)
}

if( is.null(args$SNPLIST)) 
{
  print("Return the list of PMBB imputed chunk files containing the SNPs in SNPLIST.")
  print_help(parser)
  quit(status = 1)
}
# ---------------------------------------------------------------------------- #




# ---------------------------------------------------------------------------- #
# ------                              MAIN                              ------ #
# ---------------------------------------------------------------------------- #
impute.dat.path <- "/static/PMBB/PMBB-Release-2024-3.0/Imputed/chunked_bed_files/"
files <- list.files(impute.dat.path)
files <- files[ grep("*.bed", files)]

snps <- read.table(args$SNPLIST, header = FALSE)

chunk.desc <- unlist(lapply(strsplit(files, "[.]"), function(x) x[3]))

bp <- as.integer(unlist(lapply(strsplit(snps$V1, ":"), function(x) x[2])))
chr <- unlist(lapply(strsplit(snps$V1, ":"), function(x) x[1]))

# snps need to start with chr to match pmbb
if(substr(snps$V1[1],1,1) != "c")
{
  chr <- paste0("chr", chr)
}

# create the GRange objects
snps_gr <- GRanges( seqnames = chr, 
                    ranges = IRanges(start = bp, end = bp),
                    snp_id = snps$V1)
chunks_gr <- GRanges( seqnames = unlist(lapply(strsplit(chunk.desc, "_"), function(x) x[1])),
                      ranges = IRanges(start = as.integer(unlist(lapply(strsplit(chunk.desc, "_"), function(x) x[3]))),
                                       end = as.integer(unlist(lapply(strsplit(chunk.desc, "_"), function(x) x[4])))),
                      file = files)


# find overlaps and write data
hits <- findOverlaps(snps_gr, chunks_gr, type="within")

overlap.df <- data.frame(
  snp  = mcols(snps_gr)$snp_id[queryHits(hits)],
  chr  = as.character(seqnames(snps_gr))[queryHits(hits)],
  bp   = start(snps_gr)[queryHits(hits)],
  file = mcols(chunks_gr)$file[subjectHits(hits)],
  stringsAsFactors = FALSE
)


files_to_pull <- unique(overlap.df$file)
files_to_pull <- substr(files_to_pull,1,nchar(files_to_pull)-4)
write.table(files_to_pull, "files_to_pull.txt", col.names=FALSE, row.names = FALSE, quote = FALSE)
# ---------------------------------------------------------------------------- #
