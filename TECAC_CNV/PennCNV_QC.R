#This script is used to read in penncnv output (.qcsum), and conduct probe
#based QC (part1) and sample-based QC (part2)
#This script also generates graphs summarizing the QC results based on call size and distributiions. Therefore it's important to specify an identifier for all graphs.


#-------Load Library---------
library(tidyverse)
library(stringr)
library(data.table)
library(ggplot2)

#--------Reference Data--------
grch37.ref.dat = data.frame( chromosome = c(seq(1:22), "X", "Y"),
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

#-------------------Functions------------
# ------------------------- in.telomere ---------------------------- #
in.telomere <- function( start.pos, end.pos, q.tel.start, OFFSET )
  # input: start.pos (integer), the starting position of the segment
  # end.pos (integer), the ending position of the segment
  # q.tel.start (integer), the starting position of the q-arm telomere (from ref data)
  # OFFSET (integer), the offset in bp. exclude all snps with in OFFSET of the telomere.
  #
  # output: boolean value
{
  return(start.pos <= 10000 + OFFSET | end.pos >= q.tel.start - OFFSET)
}

# --------------------------- in.centromere ------------------------- #
in.centromere <- function (start.pos, end.pos, ref.cent.start, ref.cent.end)
{
  return(start.pos < ref.cent.start & end.pos > ref.cent.end |
           start.pos > ref.cent.start & start.pos < ref.cent.end |
           end.pos > ref.cent.start & end.pos < ref.cent.end)
}

#----------------------Probe_based QC------------------
#Read in the .rawcnv file with all penncnv calls.
args = commandArgs(trailingOnly =  TRUE)
if( length(args) < 5 )
{
  print("USAGE:: ")
  print("PennCNV_QC.R RAWCNVFILE OUTFILE QCSUM PHENOFILE OUTNAME")
  print("RAWCNVFILE = path of the .rawcnv file to be QC")
  print("OUTFILE = path of the filtered .rawcnv file")
  print("QCSUM = path to the .qcsum file")
  print("PHENOFILE = path to the phenotype file")
  print("OUTNAME = prefix used in all QC graphs as identifier")
  print("")
  stop()
}
RAWCNVFILE  = args[1]
OUTFILE = args[2]
QCSUM =args[3]
PHENOFILE = args[4]
OUTNAME = args[5]

dat <- as_tibble(read.table(RAWCNVFILE, header = F,
                  strip.white = T))
start_n <- dim(dat)[1]

#Create columns storing chromosome, cnv start position and cnv end position
#,respectively
dat$chr <- str_split(dat$V1, ":", simplify = T)[,1]
dat$chr <- as.numeric(gsub("chr", "", dat$chr))
dat$start <- str_split(dat$V1, "-", simplify = T)[,1]
dat$start <- as.numeric(str_split(dat$start, ":", simplify = T)[,2])
dat$end <- as.numeric(str_split(dat$V1, "-", simplify = T)[,2])

#Filter samples in telomere or in centromere
dat$in.telomere <- F
dat$in.centromere <- F
OFFSET <- 1000000 # 1mbp by default

for (i in c(1:22)) {
  chr = i
  ref.dat <- grch37.ref.dat[grch37.ref.dat$chromosome == chr,]
  dat$in.centromere[dat$chr == chr] <- mapply(in.centromere, dat$start[dat$chr == chr], dat$end[dat$chr == chr],
                                              ref.dat$centromere.start - OFFSET, ref.dat$centromere.end + OFFSET)
  dat$in.telomere[dat$chr == chr] <- mapply(in.telomere, dat$start[dat$chr == chr], dat$end[dat$chr == chr],
                                            ref.dat$q.telomere.start, OFFSET)
}

dat <- dat[ !dat$in.centromere & !dat$in.telomere, ]

#Filter calls in HLA region on chr6
# remove cnvs in HLA
# https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37
# chr6, 28477797-33448354
dat$inHLA <- dat$start >= 28477797 & dat$end <= 33448354 & dat$chr == 6
dat <- dat[ !dat$inHLA, ]

dat <- dat %>% select(!contains("in"))

#Filter out calls containing less than 5 snps. 
dat$n <- as.numeric(str_split(dat$V2, "=", simplify = T)[,2])
dat <- dat[dat$n >= 5,]
print("Part 1 QC complete.")
print(paste("Starting with", start_n, "snps", sep = " "))
print(paste("Removed", start_n - dim(dat)[1], "snps", sep = " "))
print(paste("A total of", dim(dat)[1], "calls were kept.", sep = " "))
#---------------------------Sample-based QC (Part2) ----------------------
#Read in the qcsum file containing sample metrics.
start_n <- dim(dat)[1]
qcsum <- as_tibble(fread(QCSUM))
qcsum$LRR_SD <- as.numeric(qcsum$LRR_SD)
qcsum$BAF_drift <- as.numeric(qcsum$BAF_drift)
qcsum <- qcsum[!is.na(qcsum$LRR_SD) &!is.na(qcsum$BAF_drift),]
qcsum$sample <- mapply(function(i) str_split(basename(i), ".txt", simplify = T)[,1],
                       qcsum$File)
#Get the 95th percentile of LRR_SD and BAF_Drift
lrr_threshold <- quantile(qcsum$LRR_SD, 0.95)
baf_threshold <- quantile(qcsum$BAF_drift, 0.95)

#Filter samples beyond the thresholds
qcsum <- qcsum[qcsum$LRR_SD < lrr_threshold & qcsum$BAF_drift < baf_threshold,]
dat <- dat %>% dplyr::filter(V5 %in% qcsum$File)
#Filter out samples with more than 50 calls
dat$sample <- str_split(basename(dat$V5), ".txt", simplify = T)[,1]
#Read in phenotype data
pheno <- as_tibble(fread(PHENOFILE))
dat$pheno <- as.character(mapply(function(i) pheno$PHENO[pheno$BID == i],
                                 dat$sample))
dat$pheno <- mapply(function(i) gsub("1","control", i), dat$pheno)
dat$pheno <- mapply(function(i) gsub("2", "case", i), dat$pheno)
#Get number of calls in each sample
sample_counts <- count(dat, sample)
sample_counts$pheno <- as.character(mapply(function(i) pheno$PHENO[pheno$BID == i],
                                           sample_counts$sample))
sample_counts$pheno <- mapply(function(i) gsub("1","control", i), sample_counts$pheno)
sample_counts$pheno <- mapply(function(i) gsub("2", "case", i), sample_counts$pheno)

#Filter out samples with more than 50 calls
#Remove samples with more than 50 calls. (About 98th percentile)
sample_counts <- sample_counts[sample_counts$n <= 50,]
dat <- dat[dat$sample %in% sample_counts$sample,]


print("Part 2 QC complete.")
print(paste("Starting with", start_n, "snps", sep = " "))
print(paste("Removed", start_n - dim(dat)[1], "snps", sep = " "))
print(paste("A total of", dim(dat)[1], "calls were kept.", sep = " "))
print(paste(dim(sample_counts)[1], "samples are kept after QC.", sep = " "))

dat_write <- dat %>% select(contains("V"))
write.table(dat_write, OUTFILE,
            col.names = F, row.names = F, sep = "\t", quote = F)

#-------------------------------Plotting-----------------------------
#First, we can examine number of calls on each chromosome.
dat$chr <- as.character(dat$chr)
order <- as.character(c(1:22))
dat <- dat %>% mutate(chr = factor(chr, levels = order))
dat$chr <- as.numeric(dat$chr)
chrPlotName <- paste0( OUTNAME, "_chrDistribution.png")
png(chrPlotName, width = 800, height = 500)
ggplot(dat, aes(chr))+
  geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]))+
  scale_y_continuous(labels=scales::percent)+
  ylab("Percentage in each group")+
  facet_grid(pheno~.)+
  labs(title="Distribution of calls on each chromosome.")
dev.off()
  

#We can also examine the distribution of #of snps in a call.
snpsPlotName <- paste0(OUTNAME, "_snpsDistribution.png")
png(snpsPlotName, width = 1200, height = 800)
ggplot(dat, aes(n)) + 
  geom_bar(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]), color = "Black", fill = "white", alpha = 0.05) +
  scale_y_continuous(labels=scales::percent)+
  facet_grid(pheno~.)+
  xlab("Number of snps in each call")+
  ylab("Percentage in total")+
  labs(title="Number of SNPs in each CNV call.")+
  xlim(10,200)
dev.off()
 


#Next, check the distribution of # of bps in a call
dat$length <- dat$end - dat$start

bpPlotName <- paste0(OUTNAME, "_bpDistribution.png")
png(bpPlotName, width = 1200, height = 800)
ggplot(dat, aes(x=pheno, y=length)) + 
  geom_violin()+
  geom_boxplot(width=0.1)+
  scale_y_log10()+
  #facet_grid(pheno~.)+
  xlab("Phenotype")+
  ylab("Length (log10)")+
  labs(title="Length of each CNV call as in bp.")
dev.off()
 
#Plot out the distribution of number of calls
#in a sample
ncallsPlotName <- paste0(OUTNAME, "_ncallsDistribution.png")
png(ncallsPlotName, width = 800, height = 500)
ggplot(sample_counts, aes(n))+ 
  stat_count(aes(y=(..count..)/tapply(..count.., ..PANEL.., sum)[..PANEL..]),color = "Black", fill = "white", alpha = 0.5) + 
  scale_y_continuous(labels=scales::percent)+
  facet_grid(pheno~.)+
  labs(title="Number of calls in each sample",
       x = "Number of calls",
       y= "Percentage in total")
dev.off()
 

