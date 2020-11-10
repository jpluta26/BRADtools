

#!/usr/bin/env Rscript
# pluta 9/5/19

# script used for annotation
# PAINTOR needs a non-singular LD matrix
# if two snps have LD of 1, the matrix will be singular. need to remove any of these
# snps from the credible set and the ld matrix
args = commandArgs(trailingOnly = TRUE)
if( length(args) < 1)
{
	stop('need to provide arguments: INLDFILE, the ld matrix from the previous step')
}

# the LD matrix from the previous step 
INLDFILE=args[1]
REFSNP <- gsub("_", ":", strsplit(INLDFILE, "[.]")[[1]][1])
ld.dat <- read.table(INLDFILE, header = F)
snps <- read.table("snplst1", header= F, as.is = T)
colnames(ld.dat) <- snps$V1
row.names(ld.dat) <- snps$V1

# get the list of snps to remove
# search the ld matrix starting from one cell to the right of each diagonal
# (bc diagonal entries will always be 1)
rm.lst <- c()
nsnp = dim(ld.dat)[1]
for( i in 1:(nsnp - 1))
{
		for( j in (i + 1):nsnp)
		{
				if( ld.dat[i,j] == 1)
				{
					if( colnames(ld.dat)[j] == REFSNP)
					{ 
						rm.lst <- c(rm.lst, rownames(ld.dat)[i])
					} else
					 {					
						rm.lst <- c(rm.lst, colnames(ld.dat)[j])
					 }
				}
		}	
}

# the new ld matrix
new.ld.dat <- ld.dat[!(rownames(ld.dat) %in% rm.lst), !(colnames(ld.dat) %in% rm.lst)]

n = dim(new.ld.dat)[1]

# quick error check
for( i in 1:n)
{
	if( new.ld.dat[i,i] != 1)
	{
			stop("not all diagonal entries are 1, something went wrong")
	}
}

# write the new ld matrix and the list of snps it includes
# if there was no change, the ensuing script will detect that
write.table(new.ld.dat, INLDFILE, col.names = F, row.names = F, quote = F, append = F, sep = " ")
write.table(colnames(new.ld.dat), "snplst2", col.names = F, row.names = F, quote = F, append = F)
