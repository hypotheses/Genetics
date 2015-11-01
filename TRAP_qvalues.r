## Calculate q-values of TRAP GWAS RESULTS
PVALUE="z:/SIM1/SET1C/gwasResult/PHE1C1_GWAS.txt"
batch::parseCommandArgs()
## INSTALL QVALUE IF NOT PRESENT
## source("http://bioconductor.org/biocLite.R")
## biocLite("qvalue")

library(qvalue)
FDR<-numeric(100)
for (REP in 1:100) {
  PVALUE=paste("z:/SIM1/SET1C/gwasResult/PHE1C",REP,"_GWAS.txt",sep="")
  message("READING TRAP GWAS RESULT ",PVALUE)
  data = read.table(PVALUE,header=FALSE,col.names=c("CHR","SNP","BP","TILE","P"),as.is=TRUE)
  save(data,file=paste("z:/SIM1/SET1C/gwasResult/PHE1C",REP,"_GWAS.RData",sep=""),compress=TRUE)
  qobj <- qvalue(data$P)
  ##hist(qobj$qvalues); hist(data$P)
  ##qsummary(qobj,cuts=1e-4)
  FDR[REP] <- sum(qobj$qvalues<1e-4)
  rm(data);gc(reset=TRUE)
}