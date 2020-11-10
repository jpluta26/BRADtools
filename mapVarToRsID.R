library(data.table)

# pluta 11/19/19
# given a list of snps in chr:pos format, the script will find the corresponding rsID
# input: SNPLST, name of the file containing the list of snps
# output: OUTFILE, a file with two columns, the variant name and rsID

# this will only work on the hpc

args = commandArgs(trailingOnly = TRUE)
if(length(args) < 2)
{
  print('need to provide arguments: ')
  print('SNPLST, a text file of snps in chr:pos format')
  print('OUTFILE, a string of the output file name')
  stop()
}

INFILE  = args[1]
OUTFILE = args[2]
dat <- as.data.frame(read.table(INFILE, header = F))

dat$MarkerName <- as.character(dat$V1)
dat$CHR <- unlist(lapply(strsplit(dat$MarkerName, ":"), function(x) x[1]))
dat$BP <- unlist(lapply(strsplit(dat$MarkerName, ":"), function(x) x[2]))

out <- c()
			
for( i in unique(dat$CHR))
{
  # change this line if running else where
	REFFILE <- paste("/project/knathanslab/REF/HHRC/HRC-chrbp-rsid-map-chr", i, ".txt", sep = "")
	ref <- as.data.frame(fread(REFFILE, header = F))

	tmp <- dat[ dat$CHR == i,]
	ind <- match(tmp$MarkerName, ref$V1)
	tmp$rsid <- ref$V2[ind]
	out <- rbind(tmp[ ,colnames(tmp) %in% c("MarkerName", "rsid")], out)
}


			
write.table(out, OUTFILE, col.names = T, row.names = F, quote = F, append = F)
