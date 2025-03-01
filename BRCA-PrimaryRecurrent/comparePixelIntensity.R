# pluta  4/22/21
# script to compare median pixel intensity for different genes 

library(openxlsx)
library(ggplot2)
library(afex)



#  ================================== functions ============================= #
#  --------------------------------- load data  -----------------------  #
#  function  to  read in data from .xlsx files and truncate to the fields of interest
loadData  <- function( sub.list )
  # input:
  # sub.list (string): root of the .xlsx file containing patient information
  #
  # output:
  # dat (data.frame): data from the xlsx file truncated to the information of interest
{
  dat <- c()
  for(  i in sub.list)
  {
    if( file.exists(paste0(i, ".xlsx")))
    {
      tmp <- read.xlsx(paste0(i, ".xlsx"), sheet = 1, colNames = T)
      tmp <- tmp[,colnames(tmp) %in% c("subject", "RAD51.(RAD51)", "CTLA-4.(CTLA-4)", "PD-1.(PD-1)")]
      
      colnames(tmp)[which(colnames(tmp) == "RAD51.(RAD51)")] <- "RAD51"
      colnames(tmp)[which(colnames(tmp) == "CTLA-4.(CTLA-4)")] <- "CTLA4"
      if( "PD-1.(PD-1)" %in% colnames(tmp))
      {
        colnames(tmp)[which(colnames(tmp) == "PD-1.(PD-1)")] <- "PD1"
      }
    } else
    {
      # some files were saved as text and require different manipulation
      # text files had a column named "Region #" which will create an offset of 1
      # eg column 1 becomes row names- need manually remove the pound sign before running
      tmp <- read.table(paste0(i, ".txt"), header =  T, sep = "\t")
      tmp <- tmp[,colnames(tmp) %in% c("subject", "RAD51..RAD51.", "CTLA.4..CTLA.4.", "PD.1..PD.1.")]
      
      colnames(tmp)[which(colnames(tmp) == "RAD51..RAD51.")] <- "RAD51"
      colnames(tmp)[which(colnames(tmp) == "CTLA.4..CTLA.4.")] <- "CTLA4"
      if( "PD.1..PD.1." %in% colnames(tmp))
      {
        colnames(tmp)[which(colnames(tmp) == "PD.1..PD.1.")] <- "PD1"
      }
    }
    
    tmp$subject <- i
    dat <-  rbind(dat,  tmp)
  }
  
  
  n <- nchar(dat$subject)
  dat$subject <-  as.factor(substr(dat$subject, n - 1, n))
  return(dat)
}
#  --------------------------------------------------------------------- #


#  ------------------------------ makePlot ---------------------------------- #
makeHistPlot <- function( dat, GENENAME, plot.title, log.transform  = FALSE )
  # input:
  #  dat (data.frame): data  frame containing RAD51 pixel intensity and subject id
  #  plot.title (string): string that is the title of the plot
  #  log.transform (binary): log transform the data or not
  #
  # output:
  #  p1 (ggplot object): plot of the data
{
  if( log.transform == TRUE)
  {
    dat[GENENAME] <- log(dat[GENENAME])
  }
  
  if(!(GENENAME %in% colnames(dat)))
  {
    stop(paste0(GENENAME, " not found"))
  }
  p1 <- ggplot(data = dat, aes( x = unlist(dat[GENENAME]), fill = subject)) +
    geom_density(color="#e9ecef", alpha = 0.4) +
    theme_minimal() + 
    xlab(paste0(GENENAME, " Pixel Intensity")) +
    labs(fill = "subject") +
    ggtitle(plot.title) 
  
  if( log.transform == FALSE)
  {
    p1 <- p1 + xlim(0.0,  0.075) 
  }
  return(p1)
}
#  --------------------------------------------------------------------- #

#  ---------------------------- test.means ----------------------------- #
# get mean pixel intensity per subject, use as input to the wilcox test
test.means <- function( p.dat, r.dat, GENENAME)
  # input: 
  #   p.dat (data.frame): data.frame containing pixel intensity data per gene for 
  #     primary tumors
  #   r.dat (data.frame): data.frame containing pixel intensity data per gene for
  #    recurrent tumors
  #  GENENAME (string): gene name corresponding to a column of each df
  #
  # output:
  #  results of the wilcox test
{
  p.means <- aggregate(unlist(p.dat[GENENAME]), list(p.dat$subject), mean)$x
  r.means <- aggregate(unlist(r.dat[GENENAME]), list(r.dat$subject), mean)$x
  return(wilcox.test( p.means, r.means))
}
#  --------------------------------------------------------------------- #

#  --------------------------------------------------------------------- #
# get the median log-intensity for each subject, use as input to the wilcox test
test.medians <- function( p.dat, r.dat, GENENAME)
  # input: 
  #   p.dat (data.frame): data.frame containing pixel intensity data per gene for 
  #     primary tumors
  #   r.dat (data.frame): data.frame containing pixel intensity data per gene for
  #    recurrent tumors
  #  GENENAME (string): gene name corresponding to a column of each df
  #
  # output:
  #  results of the wilcox test
{
  p.medians <- aggregate(unlist(log(p.dat[GENENAME])), list(p.dat$subject), median)$x
  r.medians <- aggregate(unlist(log(r.dat[GENENAME])), list(r.dat$subject), median)$x
  return(wilcox.test( p.medians, r.medians))
}
#  --------------------------------------------------------------------- #

#  --------------------------------------------------------------------- #
#  compare primary versus recurrent pixel intensity using a mixed model
test.mixed.model <- function(p.dat, r.dat, GENENAME)
  # input: 
  #   p.dat (data.frame): data.frame containing pixel intensity data per gene for 
  #     primary tumors
  #   r.dat (data.frame): data.frame containing pixel intensity data per gene for
  #    recurrent tumors
  #  GENENAME (string): gene name corresponding to a column of each df
  #
  # output:
  #  results of the wilcox test
{
  p.dat$subject <- "primary"
  r.dat$subject <- "recurrant"
  all.dat <-  rbind(p.dat, r.dat)
  all.dat$subject <- as.factor(all.dat$subject)
  all.dat$subject <-  as.factor(all.dat$subject)
  
  fit <- lmer(unlist(all.dat[GENENAME])  ~  subject + (1|subject),  data =  all.dat)
  return(fit)
}
#  --------------------------------------------------------------------- #

#  --------------------------------------------------------------------- #
logTransformAndNormalize <- function( x )
{
  x <- log(x)
  if( any(x %in% c("Inf", "-Inf")))
  {
    x <- x[-which(x %in% c("Inf", "-Inf", "NaN"))]
  }
  return( (x - mean(x)) / sd(x) )
}
#  --------------------------------------------------------------------- #

#  --------------------------------------------------------------------- #
getValidVals <- function( x )
# log transformation might introduce Inf, -Inf, NaN- need to  drop these
# input:
#   x (numeric), vector of numeric values
{
  ind <- rep(TRUE, length(x))
  if( any(log(x) %in% c("Inf", "-Inf", "NaN")))
  {
    ind[which(log(x)  %in% c("Inf", "-Inf", "NaN"))] <- FALSE
  }
  
  return(ind)
}
#  --------------------------------------------------------------------- #


# -------------------------- makeCompPlot ------------------------------------------- #
makeCompPlot <- function( dat, GENENAME, plot.title )
# plot boxplot of each subject, colored by group (primary/recurrent)
# p-value is from wilcoxon test comparing median-log
# input:
#   dat (data.frame), data frame containing subject, group, and pixel intennsity per gene
#   GENENAME (string), column name of the gene of interest
#   plot.title (string), the title of the plot
#
# output:
#   p (ggplot.obj), the plot
{

  out <- test.medians(dat[dat$group == "primary",],
                      dat[dat$group == "recurrent",],
                      GENENAME)
  
  p <- ggplot(data = dat, aes(x = subject, y =  log(unlist(dat[GENENAME])), fill = subject, color = group)) +
    geom_boxplot() +
    scale_color_manual(values = c("gray", "black"))  +
 #   scale_fill_manual(values = c("blue", "blue1", "blue2", "blue3", "red", "red1", "red2", "red3")) +
    theme_minimal() +
    ylab( paste0("log(", GENENAME, ")")) +
    labs(title = plot.title,
         subtitle = paste0("p = ", round(out$p.value, 4)))
  p
}
# ------------------------------------------------------------------------------ #
#  ================================================================================ #


# ======================================== main ======================================= #

# --------- TMA1

# load data for each group and do some simple preprocessing
setwd("~/Documents/nathansonlab/forDana/TMA1")
p.list <- c("P1", "P2", "P3")
p.list <- paste0("TMA1_primary_CTLA-4/", p.list)
p.dat <- loadData( p.list )

r.list <- c("R1", "R2", "R3", "R4")
r.list <-  paste0("TMA1_recurrence_CTLA-4/", r.list)
r.dat <- loadData( r.list )

p.dat$group <- "primary"
r.dat$group <- "recurrent"

dat <- rbind(p.dat, r.dat)

p1 <- makeCompPlot( dat, "CTLA4", "CTLA4 Pixel Intensity")
png("TMA1 - CTLA.png")
print(p1)
dev.off()

p.list <- c("P1", "P2", "P3")
p.list <- paste0("TMA1_primary_PD-1/", p.list)
p.dat <- loadData( p.list )

r.list <- c("R1", "R2", "R3", "R4")
r.list <-  paste0("TMA1_recurrence_PD-1/", r.list)
r.dat <- loadData( r.list )

p.dat$group <- "primary"
r.dat$group <- "recurrent"

dat <- rbind(p.dat, r.dat)
p2 <- makeCompPlot( dat, "PD1", "PD1 Pixel Intensity")
png("TMA1 - PD1.png")
print(p2)
dev.off()


p.list <- c("P1", "P2", "P3")
p.list <- paste0("TMA1_BRCA2Isoform_1/", p.list)
p.dat <- loadData( p.list )

r.list <- c("R1", "R2", "R3", "R4")
r.list <-  paste0("TMA1_BRCA2Isoform_3/", r.list)
r.dat <- loadData( r.list )

p.dat$group <- "primary"
r.dat$group <- "recurrent"

dat <- rbind(p.dat, r.dat)
p3 <- makeCompPlot( dat, "RAD51", "RAD51 Pixel Intensity")
png("TMA1 - RAD51.png")
print(p3)
dev.off()

# ------- TMA2
setwd("~/Documents/nathansonlab/forDana/TMA2")
# primary tumors
p.list <- paste0("primary/", c("p1", "p2", "p3", "p4"))
p.dat <- loadData( p.list )


# recurrent tumors
r.list<- paste0("recurrent/", c("r1", "r2", "r3", "r4"))
r.dat  <- loadData( r.list )

p.dat$group <- "primary"
r.dat$group <- "recurrent"

dat <- rbind(p.dat, r.dat)

p1 <- makeCompPlot( dat, "CTLA4", "CTLA4 Pixel Intensity")
p2 <- makeCompPlot( dat, "RAD51", "RAD51 Pixel Intensity")
p3 <- makeCompPlot( dat, "PD1", "PD1 Pixel Intensity")

png("TMA2 - CTLA4.png")
print(p1)
dev.off()

png("TMA2 - RAD51.png")
print(p2)
dev.off()

png("TMA2 - PD1.png")
print(p3)
dev.off()


# ----- TMA3
setwd("~/Documents/nathansonlab/forDana/TMA3")
p.list <- paste0("primary/p", c("0", "4", "5", "7"))
p.dat <- loadData( p.list )

r.list <- paste0("recurrent/r", c("0", "4", "5", "7"))
r.dat <- loadData( r.list )

p.dat$group <- "primary"
r.dat$group <- "recurrent"

dat <- rbind(p.dat, r.dat)


p1 <- makeCompPlot( dat, "CTLA4", "CTLA4 Pixel Intensity")
png("TMA3 - CTLA4.png")
print(p1)
dev.off()

p2 <- makeCompPlot( dat, "RAD51", "RAD51 Pixel Intensity")
png("TMA2 - RAD51.png")
print(p2)
dev.off()

# TMA3 - PD1
p.list <- paste0("primary_PD1/", c("p2", "p5", "p6", "p9"))
p.dat <- loadData( p.list )
r.list <- paste0("recurrence_PD1/", c("r3","r4","r7","r8"))
r.dat <- loadData( r.list )
p.dat$group <- "primary"
r.dat$group <- "recurrent"

dat <- rbind(p.dat, r.dat)

p1 <- makeCompPlot( dat, "PD1", "PD1 Pixel Intensity")
png("TMA3 - PD1.png")
print(p1)
dev.off()
# ===================================================================================== #
