library(RIdeogram)

# script to create ideogram and annotate with novel/replicated/not replicated hits

data(human_karyotpe, package = "RIdeogram")
data(gene_density, package = "RIdeogram")

dat <- read.table("/Users/johnpluta/Desktop/tgct_hits_track.csv", header = T, stringsAsFactors = F, sep = ",")

# colnames have to be exactly right or ideogram crashes
# need to use box/circle/triangle as shape parameter
colnames(dat) <- c("Type",	"Shape",	"Chr",	"Start",	"End",	"color")
ideogram(karyotype = human_karyotype, 
         overlaid = gene_density, 
         label = dat, 
         label_type = "marker")
convertSVG("chromosome.svg", device = "png")
