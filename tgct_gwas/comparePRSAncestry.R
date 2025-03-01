# pluta 4/8/21
# compare PRS scores across phenotypes and ancestry groups

setwd("~/nonwhitePRS")
library(ggplot2)

# ------
# prepare datasets-  attach phenotype data to PRS scores for a given ancestry group
prepData <- function(PRSFILE, COVFILE, ANCESTRY)
# PRSFILE (string), name of text file containing the PRS scoring from PLINK
# COVFILE (string), name of text file containing phenotype information
# ANCESTRY (string), ancestry of the group, used to label plots
{
  dat <- read.table(PRSFILE, header = true())
  dat$Ancestry <- ANCESTRY
  cov <- read.table(COVFILE, header  = TRUE)
  cov  <- cov[ match(dat$IID, cov$IID),]
  dat$PHENO <- as.integer(as.factor(cov$pheno))
  
  # recode to 0 = control, 1 = cases
  dat$PHENO[ dat$PHENO  == 2 ] <- 0
  return(dat)
}
# ------

# cov file is a little different for europeans
dat  <- read.table("euro.prs.scores", header = T)
dat$Ancestry <- "european"
cov <- read.table("case_ctl_cauc.fam", header = F)
ind <- match(dat$IID, cov$V2)
dat$PHENO <- cov$V6 - 1

afr  <- prepData("afr.prs.profile", "afr.cov", "african")
hisp <- prepData("hisp.prs.profile", "hisp.cov", "hispanic")
asian <- prepData("asian.prs.profile", "asians.cov", "asian")

# subsample european subjects so the frequencies are comparable to the other groups
set.seed(2352)
dat.sub <- dat[ sample.int(10608, 500),]

all  <- rbind(dat.sub, afr, hisp, asian)


# create 4 panel plot comparing histograms of cases and controls for each ancestry group
p1 <- ggplot(data = all, aes( x = SCORESUM, fill = as.factor(PHENO))) +
  geom_histogram(color="#e9ecef", alpha = 0.6, bins = 20, position = "identity") +
  theme_minimal() +
  labs(fill = "Phenotype") +
  ggtitle("PRS Scores - All")

p4 <- p1 + facet_wrap( ~  Ancestry, nrow = 2,  ncol  =  2)  + geom_vline(xintercept =  mean(all$SCORESUM))

print(p4)


