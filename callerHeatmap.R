# pluta 12/3/19
# create a heat map comparing callers


#!/usr/bin/env Rscript

library(data.table)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
if( length(args) < 2)
{
  print("need to provide 2 arguments: INFILENAME OUTFILENAME")
  print("INFILENAME = the caller data, eg somatic.high_confidence.hg19_multianno.report.tsv")
  stop("OUTFILENAME = the ROOT of the output image and text file, eg TEST for TEST.png and TEST-stats.txt")
} else {
  INFILE <- args[1]
  OUTFILE <- args[2]
}
dat <- read.table(INFILENAME, header = TRUE, sep = "\t", as.is = TRUE)
#dat <- read.table("somatic.high_confidence.hg19_multianno.report.tsv", header = TRUE, sep = "\t", as.is = TRUE)

caller.list <- c("MUTECT2", "STRELKA2", "VARDICTJAVA", "VARSCAN2")

dat <- dat[ ,colnames(dat) %in% caller.list]

# recode anything thats not PASS or . to FAIL
for( i in caller.list )
{
  dat[[i]][ !(dat[[i]] %in% c("PASS", ".")) ] <- "FAIL"
}

dat <- as.data.frame(apply(dat, 2, as.factor))


# statistically compare the call rates of each pair of callers
# use a fisher test- might be a smarter way to do this
comp.mat <- matrix(0, nrow = 4, ncol = 4)

for( i in 1:length(colnames(dat)) )
{
  for( j in 1:length(colnames(dat)) )
  {
    comp.mat[i,j] <- fisher.test(
                table( dat[[ colnames(dat)[i] ]] , 
                 dat[[ colnames(dat)[j] ]] ) )$p.value
  }
}

diag(comp.mat) <- 1
rownames(comp.mat) <- caller.list 
colnames(comp.mat) <- caller.list

# format data for plotting
# -convert to long format
dat <- melt(dat, measure.vars = colnames(dat))
n <- dim(dat)[1] / 4
dat$x <- rep( seq(1:n), 4)
dat$value <- as.factor(dat$value)
levels(dat$value) <- c("NO CALL", "FAIL", "PASS")

p1 <- ggplot(data = dat, aes(x = x, y = variable, fill = value)) + 
  geom_tile() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()
         ) +
  guides(fill = guide_legend(title= "Call Status")) +
  ylab("Caller") +
  scale_fill_manual(labels = c("NO CALL", "FAIL", "PASS"), values = c("blue", "red", "green"))

png(paste0(OUTFILENAME, ".png"))
print(p1)
dev.off()

write.table(comp.mat, paste0(OUTFILENAME, "-stats.txt"), row.names =T, col.names = T, append = F, quote = F)
