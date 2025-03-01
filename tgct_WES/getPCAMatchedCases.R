
# pluta 8/4/21
# determine caucasian v non-caucasian with k-means clustering and plot
# need to manually update the WES db after each set of subjects

library(ggplot2)
library(dplyr)
library(broom)


# --------------------------- getDist -------------------------------- #
# function to get distance between two points on a cartesian plane
getDist <- function(x1, y1, x2, y2)
{
  d = sqrt(  ((x2 - x1)^2) + ((y2 - y1)^2)    )
  return(d)
}
# -------------------------------------------------------------------- #





# ------------------------- main ------------------------------------- #
args = commandArgs(trailingOnly=TRUE)
if( length(args) < 2 )
{
  print("need to provide 2 arguments: INFILE OUTFILE")
  print("INFILE: the list of case IDs to match with controls")
  print("OUTFILE: name of the output file")
  stop()
}

INFILE  = args[1]
OUTFILE = args[2]

cagmap  <-  read.table("tecac_id_map.txt", header = TRUE, sep  =   ",")

# subjects already selected for WES
wes <- read.table("WES_subjects_db.csv", header =  T, sep  = ",")

# eigenvectors from genotype data
ev = paste("EV", seq(1:10), sep="")
evec.dat <- read.table("~/Documents/nathansonlab/tecac-manuscript/PCA/merged.evec", header=FALSE, col.names= c("IID", ev, "Pheno"))
evec.dat.bak <- evec.dat # save a copy for plotting
evec.dat$PRE  <- substr(as.character(evec.dat$IID), 1, 2)

# truncate by qc subjects
fam  <-  read.table("~/Documents/nathansonlab/tecac-manuscript/case_ctl_cauc.fam", header = F)
evec.dat <- evec.dat[ evec.dat$IID %in% fam$V2,]

# no biospecimens for MDA
evec.dat <- evec.dat[ -grep("MDA", evec.dat$IID),]

# exclude cases already controlled
evec.dat <- evec.dat[ !(evec.dat$IID %in% wes$COLLABORATOR.ID.JOHNS.ID[wes$CASE_CTRL.TEXT == "CTRL"]),]

# controls  with famhx of tgct- dont use these
famhx.ctrls <- read.table("famhx_ctrls.txt", header = F)
evec.dat <- evec.dat[ !(evec.dat$IID %in% famhx.ctrls),]
ctrls <- c()
cases <- c()
nogeno <- c()


# list of cases to match
sub.list <- read.table(INFILE, header = FALSE)

# iterate through cases and  pick matching control based on PCA
for( sub.id in sub.list$V1 )
{
  if( !(sub.id %in% evec.dat$IID))
  {
    nogeno <- c(nogeno, sub.id)
  } else
  {
    PRE <- substr(sub.id,1,2)
    ind <-  which(evec.dat$IID == sub.id)
    
    # caculate the distance of every subject from sub.id
    evec.dat$sub.dist <- getDist(evec.dat$EV1[ind], evec.dat$EV2[ind], evec.dat$EV1, evec.dat$EV2)
    tmp <- evec.dat[ evec.dat$Pheno == "Control",]
    
    # take first control with minimum distance from case (multiple subjects may have a distance of 0)
    min.d.sub <- tmp$IID[which(tmp$sub.dist == min(tmp$sub.dist))][1]
    ctrls <- c(ctrls, min.d.sub)
    cases <- c(cases, sub.id)
    
    # remove case and selected control from data
    evec.dat <- evec.dat[ -ind, ]
    evec.dat <- evec.dat[ -which(evec.dat$IID  == min.d.sub),]
  }
 
}


out.df  <- data.frame( CASES = cases, 
                         CAGID.CASES = cagmap$CAG.ID[match(cases, cagmap$TECAC.ID)],
                         MATCHEDCTRL  =  ctrls, 
                         CAGID.CTRLS = cagmap$CAG.ID[match(ctrls, cagmap$TECAC.ID)])

nogeno.df <- data.frame( IID = nogeno, CAGID = cagmap$CAG.ID[match(nogeno, cagmap$TECAC.ID)])


if(any(duplicated(out.df$MATCHEDCTRL)))
{
  stop("duplicate controls found")
}

if(any(duplicated(out.df$CASES)))
{
  stop("duplicate cases found")
}

if(any(duplicated(out.df$CAGID.CASES)))
{
  stop("duplicate CAGID cases found")
}

if(any(duplicated(out.df$CAGID.out)))
{
  stop("duplicate CAGID controls found")
}



evec.dat <- evec.dat.bak[ evec.dat.bak$Pheno != 3,]

# plot of cases and controls, to check distance
p1 <- ggplot(data = evec.dat[ evec.dat$IID %in% c(cases, ctrls),], 
             aes(x = EV1, y = EV2, color = as.factor(Pheno))) +
             theme_minimal() +
             scale_color_manual(labels = c("Case", "Control"), values = c("red", "blue")) +
             geom_point(alpha = 0.5, size=3)
png("matchedctrls.png")
print(p1)
dev.off()

tmp <- evec.dat[ evec.dat$IID %in% c(cases, ctrls),]

write.table(out.df, OUTFILE, col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(nogeno.df, "nogenosubs.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
# -------------------------------------------------------------------- #



