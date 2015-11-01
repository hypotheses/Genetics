# Rank based inverse normal using Rankit.
# Based on Anthony J. Bishara and James B. Hittner, 
# Testing the Significance of a Correlation With Nonnormal Data: 
# Comparison of Pearson, Spearman, Transformation,and Resampling Approaches. 
# Psychological Methods, 2012, 17(3), 399-417.
RINfun <- function(yorig){
  yranksrank(yorig)
  tempp(yranks-.5)/(length(yranks))
  return(qnorm(tempp))
}