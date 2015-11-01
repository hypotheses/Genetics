################################################################################
## PROGRAM: read.mldose.r
## BY: Bhoom Suktitipat, MD, PhD
## DATE: 03/03/10
################################################################################
## GOAL: read MACH mldose data and return a data frame containing genotype
## OPTION:
##       - mldosefile: genotype file from MACH "mldose file"
##       - mlinfofile: mlinfo file
##       - pedfile: pedigree structure file
##       - sep.gen: default separator in genfile " "
##       - sep.ped: default seprator in pedfile "," 
## OUTPUT: return genotype data frame 
################################################################################# 
read.mldose <- function(mldosefile,mlinfofile,pedfile,sep.ped,sep.gen) {
  message("Reading in White Split Imputed genotype data from mldose file\nUpdated 02/02/2011 ")
  ped.dat <- read.table(mldosefile, header = FALSE, sep = sep.gen, as.is=T)
  ped.dat[,1] <- matrix(unlist(strsplit(ped.dat[,1],"->")),byrow=T,ncol=2)[,2]
  ped.dat[,2] <- NULL
  message("Reading mlinfo files")
  mlinfo <- read.table(mlinfofile,header=T,as.is=T)
  snp.names <- paste(gsub("X","",mlinfo[,1]),mlinfo[,2],sep="_")
  names(ped.dat) <- c("id",snp.names)
  pedigree <- read.table(pedfile, header = TRUE, sep = sep.ped)
  genotype <- merge(pedigree, ped.dat, by = "id")
  names(genotype)[1:5] <- c("IID","FID","PAT","MAT","SEX")
  print("Done reading in imputed genotype data")
  return(genotype)
}
