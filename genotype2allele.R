geno2allele <- function(x,allele.sep="",snp.name,allele.name=paste(".a",1:2,sep="")) {
  ## input genotype in one column format e.g. AA, AG, GG
  ## convert to ?[allele.sep]?
  x <- gsub(allele.sep,"",x) ## remove allele separator
  allele1 <- substr(x,1,1)
  allele2 <- substr(x,2,2)
  #   if (length(snp.name)) {
  #     alleleName <- paste(snp.name,allele.name,sep="")
  #   } else {
  #     alleleName <- paste("SNP",allele.name,sep="")
  #   }
  allele <- data.frame(allele1,allele2)
  names(allele) <- c("a1","a2")
  return(allele)
}

# Example of converting a data frame of genotype to a data frame of alleles
# geno <- list()
# for (i in 3:9) {
#   geno[[names(pon123)[i]]] <- geno2allele(pon123[,i],snp.name=names(pon123[i]))
# }
# head(geno[[1]])
# geno <- do.call(cbind,geno)
# genotype <- cbind(pedigree,geno)