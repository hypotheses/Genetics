##Bhoom Suktitipat, MD,PhD
##September 7, 2012
## calculate minor allele frequency of genotype coded as 0,1,2 following additive model
## as the number of minor allele

recodeA.freq <- function(x) {
  tx=table(x)
  if (length(tx < 3)) {
    if( is.na(tx["2"]) ) tx["2"]=0
    if( is.na(tx["1"]) ) tx["1"]=0
  } 
  maf = (2*tx["2"]+tx["1"])/(2*sum(tx))
  maf
}

## Examples
# geno <- list()
# P = c(.01,.05,.1,.2,.3,.4)
# for (i in 1:length(P)) {
#   p=P[i]
#   geno[[i]] <- sample(c(2,1,0),10000,replace=TRUE,prob=c(p^2,2*p*(1-p),(1-p)^2))
# }
# geno <- do.call(cbind,geno)
# apply(geno,2,recodeA.freq)
