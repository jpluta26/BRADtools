
# pluta 4/13/21
# determine caucasian v non-caucasian with k-means clustering and plot


library(ggplot2)
library(dplyr)
library(broom)


# function to get distance between two points on a cartesian plane
getDist <- function(x1, y1, x2, y2)
{
  d = sqrt(  ((x2 - x1)^2) + ((y2 - y1)^2)    )
  return(d)
}

set.seed(413153)

# hapmap demographics
HAP.dem = read.table("relationships_w_pops_051208.txt", header=TRUE)

# eigenvectors
ev = paste("EV", seq(1:10), sep="")

evec.dat = read.table("merged.evec", header=FALSE, col.names= c("IID", ev, "Pheno"))
levels(evec.dat$Pheno) <- c("Parent", "Case", "Control")
evec.dat$PRE = substr(evec.dat$IID, 1, 2)

# caucasians are subjects within n.sd of the center of the caucasian cluster
n.sd <-  6

# perform clustering on EV1 and EV2
# 4 clusters for the 4 major genomic ancestry groups (caucasian, asian, hispanic african)
ev.cls <- kmeans( evec.dat[,2:3], 4, nstart= 20)

# largest cluster
l = which(ev.cls$size == max(ev.cls$size))

# the center of each cluster
ev.cls$cluster <- as.factor(ev.cls$cluster)
evec.dat$cluster = ev.cls$cluster

# data frame of xy coordinates for center of each cluster
centers <- data.frame(x=ev.cls$centers[,1], y=ev.cls$centers[,2], cluster = c(1:4))

# cluster center
cls1.cnt = ev.cls$centers[l,]

# compute distance of each point from center of the largest cluster (europeans)
evec.dat$dist = getDist(centers[l,1], centers[l,2], dat$EV1, dat$EV2)
d.sd   <- sd(evec.dat$dist[evec.dat$cluster == l])


# radius of the inclusion region is # of standard deviations of distance from
# the center of the region
r <- d.sd * n.sd
xc <- centers[l,1]
yc <- centers[l,2]

# check that this aligns with cluster order
levels(ev.cls$cluster) <- c("Hispanic", "Caucasian", "Asian", "African" )

p1 <- ggplot(evec.dat, aes(x=EV1, y=EV2, color= ev.cls$cluster)) +  
geom_point(alpha=0.3) +
  
  geom_point(data=evec.dat[which(evec.dat$dist <= n.sd * d.sd),], aes(x=EV1, y=EV2), col="red") +
  annotate("path", x = xc + r*cos(seq(0,2*pi, length.out=100)),
                   y = yc + r*sin(seq(0,2*pi, length.out=100))) +
  geom_point(data=centers, aes(x=x, y=y), col = "black") +
  geom_point(data = evec.dat[evec.dat$PRE == "NA",], aes(x=EV1, y=EV2), col="yellow") +
  ggtitle("PCA") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1), title = "Ancestry")) +
  theme_minimal() 

png("k-means-PCA.png", height=800, width=800)
print(p1)
dev.off()
























