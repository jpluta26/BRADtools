# -------------------------------------------------------------------------- #
# plots two roc curves, and statistically tests for a difference between them
compare2ROC <- function( roc.dat1, roc.dat2, title, grp1.name, grp2.name)
# input: roc.dat1, roc.dat2 (list), list objects that are output from function "getROCstats"
#        title (character), the title of the plot
#        grp1.name, grp2.name (character), the name of the groups, added to the legend
{
  compROC = roc.test(roc.dat2[[3]], roc.dat1[[3]], paired = TRUE, method = "delong", alternative = "greater")
  
  # for display purposes, can change this to be more exact
  p <- round(compROC$p.value, 3)
  if( p < 0.001 )
  {
    p.annot <- "p < 0.001"
  } else
  {
    p.annot <- paste("p = ", p, sep ="")
  }
  
  # probably should change this to a ggplot object and pass, rather than write the plot in the function
  # summary plot of the test
  png(paste(title, "png", sep = "."))
  text1 <- paste(grp1.name, " AUC = ", round(roc.dat1[[2]][[1]][1], 2), sep = "")
  text2 <- paste(grp2.name, " AUC = ", round( roc.dat2[[2]][[1]][1], 2), sep = "")
  plot(roc.dat1[[1]], col = "blue", main = title)
  plot(roc.dat2[[1]], col = "red", add = TRUE)
  mtext(p.annot, side = 3)
  legend(x = 0.5, y = 0.2, cex = .9, legend = c(grp1.name, grp2.name), col = c("blue", "red"), lwd = 3)
  dev.off()
}  
# -------------------------------------------------------------------------- #
