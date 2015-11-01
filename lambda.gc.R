lambda.gc <- function(p,plot=FALSE,col = palette()[4], lcol = palette()[2], ...) {
  ## for any vector p containing p-values
  ## Convert p-value to Chi-sq and calculate genomic controls
  ## Modify from gcontrol2() from library("gap") version 1.1.3
  ## Bhoom Suktitipat
  ## 08/01/2012
  p <- p[!is.na(p)]
  n <- length(p)
  x2obs <- qchisq(p, 1, lower.tail = FALSE)
  x2exp <- qchisq(1:n/n, 1, lower.tail = FALSE)
  lambda <- median(x2obs)/median(x2exp)
  if (plot) qqunif(p, col = col, lcol = lcol, ...)
  return(lambda)
}
