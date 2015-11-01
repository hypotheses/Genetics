################################################################################
## PROGRAM: read.plink.r
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
##       - sep: if missing assuming the space delimitted format as default
##              change to "\t" if --recode --tab was used to create the PLINK file
## OUTPUT: return a list of class readPlink
##       - $ped: genotype data frame
##       - $map: map read from mapfile
##       - $dat: data type
################################################################################# 
## READ PLINK EXPORT PEDIGREE WITH MAP FILE
read.plink <- function(pedfile,mapfile,sep) {
  if (missing(pedfile)) stop("Please enter pedigree file")
  if (missing(mapfile)) stop("Please input map file in the four columns format CHR SNP CM")
  cat("Reading map",mapfile,"\r")
  map <- read.table(mapfile,col.names=c("chr","snp","cm","pos"),as.is=T,header=FALSE)
  if (!ncol(map)==4) stop("Please input map file in the four columns format CHR SNP CM")
  cat("Reading pedigree",pedfile,"\r")
  ## DEFAULT : SPACE DELIMITED FILE
  if (missing(sep)) sep=" "
  if ( sep == " ") {
    sep=" "
    ped <- as.data.frame(matrix(scan(pedfile,
                                     what=c(PID=0,IID=0,PAT=0,MAT=0,SEX=0,PHE=0,
                                       rep(character(),6+nrow(map)*2)),sep=" "),
                                ncol=6+nrow(map)*2,byrow=TRUE),stringsAsFactors=FALSE)
    mysnp<-paste(rep(map[,2],each=2),c("_1","_2"),sep="")
    names(ped) <- c("PID","IID","PAT","MAT","SEX","PHE",mysnp)
    expectedcol <- (nrow(map)*2)+6
    if ( !expectedcol == ncol(ped) ) {
      warning(paste("Number of expected column",expectedcol,"not equal number of column found",ncol(pedigree)))
    }
    map$cm <- map$pos/1e6
    dat <- data.frame("TYPE"=c("S",rep("M",nrow(map))),"NAME"=c(names(ped)[6],map$snp),stringsAsFactors=F)
  }
  ## TAB DELIMITED FILE
  if ( sep == "\t") {
    ped <- as.data.frame(matrix(scan(pedfile,
                                     what=c(PID=0,IID=0,PAT=0,MAT=0,SEX=0,PHE=0,
                                       rep(character(),6+nrow(map)*2)),sep=sep),
                                ncol=6+nrow(map),byrow=TRUE),stringsAsFactors=FALSE)
    mysnp<-map[,2]
    names(ped) <- c("PID","IID","PAT","MAT","SEX","PHE",mysnp)
    expectedcol <- (nrow(map))+6
    if ( !expectedcol == ncol(ped) ) {
      warning(paste("Number of expected column",expectedcol,"not equal number of column found",ncol(pedigree)))
    }
  }
  map$cm <- map$pos/1e6
  dat <- data.frame("TYPE"=c("S",rep("M",nrow(map))),"NAME"=c(names(ped)[6],map$snp),stringsAsFactors=F)
  plink <- list("ped"=ped,"map"=map,"dat"=dat)
  class(plink) <- "readPlink"
  return(plink)
}
