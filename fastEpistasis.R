snp.assoc <- functioN(x1,x2) {
  ## test for association between two snps using allelic association test based on Chi-squared 
  ## This is the same function implemented for --fast-epistasis --case-only in PLINK
  tx <- table(x1,x2)
  tx2 <- tx[c(1,4,7,2,5,8,3,6,9)]
  tx.a <- matrix(c(4*tx2[1]+2*tx2[2]+2*tx2[4]+tx2[5],
                   4*tx2[3]+2*tx2[2]+2*tx2[6]+tx2[5],
                   4*tx2[7]+2*tx2[8]+2*tx2[4]+tx2[5],
                   4*tx2[9]+2*tx2[8]+2*tx2[6]+tx2[5]),
                 c(2,2),byrow=TRUE)
  chisq.test(tx.a)
}
x1 <- sample(c("1/1","1/3","3/3"),size=3000,replace=TRUE,prob=c(.25,5,.25))
x2 <- sample(c("2/2","2/3","3/3"),size=3000,replace=TRUE,prob=c(.25,5,.25))
