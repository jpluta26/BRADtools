# pluta 2/27/25

# script to create a file with PRS weights and effect alleles corresponding to the genotype data.
# the output of this script is the input for the plink --score option

library(dplyr)
suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c("-a", "--ASSOC"), action="store", default = NA, type="character",
              help = "the association statistics containing weight and direction. output of METAL"),
  make_option(c("-f", "--FRQ"), action = "store", default = NA, type = "character",
              help = "the MAF (.frq) file from plink for the genotype data"),
  make_option(c("-o", "--OUT"), action = "store", default = NA, type = "character",
              help = "output file name")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if( (opt$help))
{
  print_help(opt_parser)
  stop()
}

ASSOC = opt$ASSOC
FRQ = opt$FRQ
OUT = opt$OUT

# Read the data
frq <- read.table(FRQ, header = TRUE)
dat <- read.table(ASSOC, header = TRUE)

# Filter and match SNPs
common.snps <- intersect(dat$MarkerName, frq$SNP)

# output snps that are missing from either ref or weights
write.table(frq$SNP[ !(frq$SNP %in% common.snps)], "frq-snps-missing.txt", col.names=F,quote=F,row.names=F)
write.table(dat$MarkerName[ !(dat$MarkerName %in% common.snps)], "dat-snps-missing.txt", col.names=F,quote=F,row.names=F)

frq <- frq[ frq$SNP %in% common.snps,]
dat <- dat[ dat$MarkerName %in% common.snps,]
dat <- dat %>% 
  filter(MarkerName %in% frq$SNP) %>%
  arrange(match(MarkerName, frq$SNP))

print(paste0(dim(dat)[1], " SNPs found in genotype data."))

# Ensure the order is correctly aligned
stopifnot(all(frq$SNP == dat$MarkerName))


# figure out the effect and whether positive effect is tied to major or minor allele
# this uses the effect allele so all scores will be positive
dat <- dat %>%
  mutate(EffectAllele = ifelse((Freq1 > 0.5 & Effect < 0) | (Freq1 < 0.5 & Effect > 0), frq$A1, frq$A2),
         Effect = abs(Effect))


out <- dat %>% select(MarkerName, EffectAllele, Effect)
write.table(out, OUT, row.names = F, quote = F, col.names = T)

# TODO: might as well invoke plink here
