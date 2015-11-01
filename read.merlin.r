################################################################################
## PROGRAM: read.merlin.R
## BY: Bhoom Suktitipat, MD,PhD
## DATE: Tue Dec 20 15:38:18 EST 2011
## Bug Fix 1/21/2013 -- fix read in data from T to logical "TRUE"
################################################################################
## GOAL: read pedigree file from Merlin into R 
## INPUT:
##       - pedfile: name of pedigree file. This is a linkage format file
##       - datfile: dat type file: a two column file format with 1st column
##                  describing the data type "M/S/C/T", and 2nd column
##                  for variable names
## OUTPUT: create a list of class "readMerlin" for writing using write.merlin()
##         readMerlin contains "$pedigree" and "$dat"
################################################################################
read.merlin <- function(pedfile,datfile,force=1,...) {
  dat=read.table(file=datfile,header=F,col.names=c("TYPE","NAME"),as.is=T,...)
  if ( is.logical(dat[,1]) ) {
    dat[dat[,1],1]="T"
  }
  expectedcol<-sum(dat[,1]=="M")*2+sum(dat[,1]=="S2")*2+sum(dat[,1]=="S")+sum(dat[,1]=="T")+sum(dat[,1]=="C")+5
  pedigree=read.table(file=pedfile,header=F,as.is=T,na.strings="x",...)
  ## Warn if the number of column in pedigree file is different from expected
  if ( !expectedcol == ncol(pedigree) ) {
    warning(paste("Number of expected column",expectedcol,"not equal number of column found",ncol(pedigree)))
    if (force==1) {
      message("Trailing spaces and extra-columns will be removed")
      pedigree <- pedigree[,1:expectedcol]
    }
  }
  datname<-character()
  for (i in 1:nrow(dat) ) {
    if (dat[i,1]=="M") {
      datname<-append(datname,paste(rep(dat[i,2],each=2),c("_1","_2"),sep=""))
    }
    else if (dat[i,1]!="E") {
      datname<-append(datname,dat[i,2])
    }
  }
  names(pedigree)<-c("PID","IID","PAT","MAT","SEX",datname)
  merlin <- list(ped=pedigree,dat=dat)
  class(merlin) <- "readMerlin"
  return(merlin)
}
