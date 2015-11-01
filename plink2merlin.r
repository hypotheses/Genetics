################################################################################
## PROGRAM: plink2merlin.r
## BY: Bhoom Suktitipat, MD, PhD
## DATE: 03/03/10
################################################################################
## GOAL: Convert PLINK Export in pedigree format to Melrin format with or
##    without phenotype file to merge with Phenotype file hsould have the first
##    two column as PEDID and PID following by PHENOTYPE DATA
## OPTION:
##       - pedfile: pedigree file in PLINK format from PLINK --recode with default
##                  space delimitted
##       - mapfile: default map file from PLINK --recode option
##       - phefile: phenotype file to be included in MERLIN the first two column =
##                  FID and IID
##       - outfile: output filename for MERLIN format file
################################################################################# 
plink2merlin <- function(pedfile,mapfile,phefile,outfile="plink2merlin") {
  if (missing(pedfile)) stop("Please enter pedigree file")
  if (missing(mapfile)) stop("Please input map file in the four columns format CHR SNP CM")
  cat("Reading map",mapfile,"\r")
  map <- read.table(mapfile,col.names=c("chr","snp","cm","pos"),as.is=T)
  if (!ncol(map)==4) stop("Please input map file in the four columns format CHR SNP CM")
  mysnp<-paste(rep(map[,2],each=2),c("_1","_2"),sep="")
  cat("Reading pedigree",pedfile,"\r")
  ped <- read.table(pedfile,col.names=c("famid","id","fa","mo","sex","phe",mysnp),as.is=T)
  expectedcol <- (nrow(map)*2)+6
  if ( !expectedcol == ncol(ped) ) {
    warning(paste("Number of expected column",expectedcol,"not equal number of column found",ncol(pedigree)))
  }
  if (!missing(phefile)) {
    cat("Reading phenotype file",phefile,"\r")
    phen <- read.table(phefile,header=T,as.is=T)
    ped[,names(phen)[-2:-1]] <- phen[match(ped$id,phen$id),-c(1:2)]
    cat("Saving new pedigree",paste(outfile,".ped",sep=""),"\r")
    write.table(ped,paste(outfile,".ped",sep=""),sep="\t",col.names=F,row.names=F,na="x",quote=F)
    map$cm <- map$pos/1e6
    cat("Saving new map",paste(outfile,".map",sep=""),"\r")
    write.table(map[,-4],paste(outfile,".map",sep=""),sep="\t",col.names=F,row.names=F,quote=F)
    dat <- data.frame("key"=c("S",rep("M",nrow(map)),rep("T",ncol(phen)-2)),"name"=c(names(ped)[6],map$snp,names(phen)[-2:-1]),stringsAsFactors=F)
    cat("Saving new datfile",paste(outfile,".dat",sep=""),"\r")
    write.table(dat,paste(outfile,"merlin.dat",sep=""),sep="\t",col.names=F,row.names=F,quote=F)
  }
  else {
    cat("Saving new pedigree",paste(outfile,".ped",sep=""),"\r")
    write.table(ped,paste(outfile,".ped",sep=""),sep="\t",col.names=F,row.names=F,na="x",quote=F)
    map$cm <- map$pos/1e6
    cat("Saving new map",paste(outfile,".map",sep=""),"\r")
    write.table(map[,-4],paste(outfile,".map",sep=""),sep="\t",col.names=F,row.names=F,quote=F)
    cat("Saving new dat",paste(outfile,".dat",sep=""),"\r")
    dat <- data.frame("key"=c("S",rep("M",nrow(map))),"name"=c(names(ped)[6],map$snp),stringsAsFactors=F)
    write.table(dat,paste(outfile,".dat",sep=""),sep="\t",col.names=F,row.names=F,quote=F)
  }
}
