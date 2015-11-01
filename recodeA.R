geno012 <- function(genotype,sep.allele="/") {
  ## recode a matrix of character allele of size n x 2 to 
  ## a vector of number of minor alleles vector coded as 0,1, or 2
  ## required input as a matrix of genotypes with 2 columns/marker
  if ( ncol(genotype) %% 2 ) stop("Number of alleles is not a denomination of 2")
  genotype <- as.matrix(genotype)
  genotype <- sub(sep.allele,"",genotype)
  nsnp=ncol(genotype)/2
  recoded <- matrix(ncol=nsnp,nrow=nrow(genotype))
  for (i in 1:nsnp) {
    alleles <- as.matrix(genotype[,(2*(i-1)+1:2)])
    allele <- names(sta <- sort(table(alleles))) ## find minor allele by sorting by allele frequency
    alleles.f <- as.factor(alleles) ## convert to factor for conversion
    alleles.n <- matrix(c(0,1)[alleles.f],ncol=2) ## convert alllele to number of minor allele
    recoded[,i] <- apply(alleles.n,1,sum) ## sum number of minor allele in the genotype
  }
  return(recoded)
}