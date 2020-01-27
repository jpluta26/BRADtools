library(data.table)
# from:http://genometoolbox.blogspot.com/2014/08/how-to-calculate-genomic-inflation.html
# http://rstudio-pubs-static.s3.amazonaws.com/9743_8a5f7ba3aa724d4b8270c621fdf6d06e.html
# pvals is a vector of p-values from the meta analysis, after filtering
# exclude multi-allelic variants, excess heterogeneity, and snps that did not appear
# in all sites
dat <- as.data.frame(fread("pvals", header = F))

chisq <- qchisq(1 - dat$V1, 1)
lgc <- median(chisq)/ qchisq(0.5, 1)

# lgc inflates with large samples
# rescale to a study of 1000 cases, 1000 ctrls
lgc.rescale <- 1 + (lgc - 1) * ((( 1 / n.cases ) + ( 1 / n.ctrls)) / (( 1/1000 ) + ( 1/1000)))
