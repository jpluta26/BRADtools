library(ggplot2)



# ------------------------------ get.GIF ------------------------------ #
get.GIF <- function( p.vec )
# compute genomic inflation factor
# input: p.vec (numeric), vector of p-values
# output: GIF
{
  median( qchisq(1- p.vec,1))/qchisq(0.5,1)
}
# --------------------------------------------------------------------- #


# -------------------------- gg.qq ------------------------------------ #
gg.qq <- function( p.vec, x.pos , y.pos  )
# recreation of qq from qqman, but using ggplot
# input:
#   p.vec (numeric), vector of p-values
#.  x.pos, y.pos (numeric), x/y coordinates for the position of the GIF text
# output:
#   p1, a ggplot object
{

  if( (!is.numeric(p.vec)))
      stop("input must be numeric")
  
  p.vec <- p.vec[!is.na(p.vec) & !is.nan(p.vec) & !is.null(p.vec) & is.finite(p.vec) &
                 p.vec < 1 & p.vec > 0]
  o = -log10(sort(p.vec, decreasing = FALSE))
  e = -log10(ppoints(length(p.vec)))
  
  
  txt = paste("lambda", round(get.GIF(p.vec),2), sep = " == ")
  p1 <- ggplot(data = data.frame(o = o, e = e), aes(x = e, y = o)) + geom_point() +
    theme_minimal() +
    xlab("Expected") +
    ylab("Observed") + 
    geom_abline(slope=1, intercept=0, color="red") +
     annotate("text", x = x.pos, y = y.pos, 
              label = txt, 
              parse = T,
              size = unit(10, "pt"))
  return(p1)
    
}
# --------------------------------------------------------------------- #
