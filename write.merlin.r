################################################################################
## PROGRAM: write.merlin.R
## BY: Bhoom Suktitipat, MD,PhD
## DATE: Tue Dec 20 15:38:18 EST 2011
################################################################################
## GOAL: write pedigree file for Merlin (LINKAGE FORMAT) with DAT file from
##       1) object from read.merlin
##       2) when giving both pedfile and dat file, will write out pedfile&datfile
##       3) when giving only pedfile and "type", will create datfile coresponding
##          to "type"
## INPUT: 
##       - x: a list of class readMerlin from read.merlin()
##       - prefix: the prefix of pedfile and datfile you want to create
##       - ped: name of pedigree file. This is a linkage format file
##       - dat: dat type file: a two column file format with 1st column
##                  describing the data type "M/S/C/T", and 2nd column
##                  for variable names
##       - map: Four columns format with "CHR, SNP, CM, BP"
##              OPTIONAL: If map specified, it will be writen out.
##       - type: NULL (default) use data from .dat
##               when specified as a character string same as the first column of
##               .dat
## OUTPUT: write ${prefix}.ped file and ${prefix}.dat file to the current directory
##         readMerlin contains "$pedigree" and "$dat"
################################################################################
write.merlin <- function(x,prefix,ped,dat,map,type=NULL) {
  if (missing(prefix)) stop("Please specify output file")
  ## SITUATION 1: Given readMerlin class data
  if (!missing(x)) {
    if (class(x) == "readMerlin") {
      write.table(x$ped,file=paste(prefix,"ped",sep="."),sep="\t",col.names=F,row.names=F,na="x",quote=F)
      write.table(x$dat,file=paste(prefix,"dat",sep="."),sep="\t",col.names=F,row.names=F,na="x",quote=F)
    }
  }
  ## SITUATION 2: given pedigree, map, and dat data frame  
  if (!missing(ped) & !missing(dat) ) {
    write.table(ped,file=paste(prefix,"ped",sep="."),sep="\t",col.names=F,row.names=F,na="x",quote=F)
    write.table(dat,file=paste(prefix,"dat",sep="."),sep="\t",col.names=F,row.names=F,na="x",quote=F)
    if ( !missing(map)) write.table(map,file=paste(prefix,"map",sep="."),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
  } else if (missing(ped) &  missing(dat)) {
    stop("Please specify Pedigree and dat data")
  } else  
  ## SITUATION 3: pedigree data frame available but no dat --- specify type to create dat automatically
  if (!missing(ped) & missing(dat) & length(type) > 0) {
    cat("generating dat file from pedigree data frame\r")
    ## calculate expected length of data name
    length.type <- length(type[type=="M"])*2+length(type[!type=="M"]) 
    if (length.type != ncol(ped)-5) stop("type has length shorter than pedigree data")
    type.name <- character(length(type))
    colnum <- 6
    allname <- names(ped)
    for ( i in 1:length(type) ) {
      if (type[i] == "M" ) {
        marker <- allname[colname]
        marker <- strsplit(marker,"\\_")
        type.name[i] <- marker
        colnum <- colnum+2
      } else {
        type.name[i] <- allname[colnum]
        colnum <- colnum+1
      }
    }
    dat <- data.frame("type"=type,"name"=type.name)
    ord <- order(ped[,1],ped[,2])
    ped <- ped[ord,]
    write.table(ped,file=paste(prefix,"ped",sep="."),sep="\t",col.names=F,row.names=F,na="x",quote=F)
    write.table(dat,file=paste(prefix,"dat",sep="."),sep="\t",col.names=F,row.names=F,na="x",quote=F)
    if ( !missing(map)) write.table(map,file=paste(prefix,"map",sep="."),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
  } else stop("Please specify type for dat file to write out")
}  
