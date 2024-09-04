
# --------
library(data.table)

# pluta 11/19/19
# given a list of snps in chr:pos format, the script will find the corresponding rsID
# input: SNPLST, name of the file containing the list of snps
# output: OUTFILE, a file with two columns, the variant name and rsID

# this will only work on the hpc

library(data.table)
INFILE  = args[1]
OUTFILE = args[2]

args = commandArgs(trailingOnly = TRUE)


if(length(args) < 2)
{
  print('need to provide arguments: ')
  print('INFILE, eg. denmark_rsids.txt')
  print('OUTFILE, a string of the output file name')

  stop()
}

REF="/project/knathans_tecac/REF/HRC/hg19/HRC-chrbp-rsid-map.txt"
ref <- as.data.frame(fread(REF, header = T))

dat <- as.data.frame(read.table(INFILE, header = T))

dat <- dat[ dat$rsID %in% ref$ID,]
ref <- ref[ ref$ID %in% dat$rsID,]
ref <- ref[-which(duplicated(ref$ID)),]
colnames(ref) <- c("MarkerName", "ID")

dat <- dat[ match(ref$ID, dat$rsID),]

dat$rsID <- ref$MarkerName
dat$rsID <- as.character(dat$rsID)
dat$CHR <- unlist(lapply(strsplit(dat$rsID, ":"), function(x) x[1]))
dat$BP <- unlist(lapply(strsplit(dat$rsID, ":"), function(x) x[2]))

write.table(dat, OUTFILE, col.names=T, row.names=F, quote=F)
