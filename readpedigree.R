## ------------------------------------------------------------- ##
## 12/20/2011 ##
## obsolete function replalced with individual read.datatype.r
## ------------------------------------------------------------- ##
## A bundle of code to read pedigree from various format and return a data frame
## 1. read.merlin (03/01/10) : READ MERLIN PED FILE WITH DAT FILE
#-------read.merlin(pedfile="",datfile="",...) # "sep" can be specify to tab if needed.
#-------assume two allele columns per one marker
#-------fixed (04/15/10) :read error when map file contain E END-OF-DATA
## 2. write.merlin (03/02/10) : Write PED and DAT file for MERLIN
## 3. freq (03/01/10): calculate effective minor allele frequency from the gwasGEEAdd()
## 4. splitwhite (add 03/03/10): split large white pedigree
## 5. splitblack (add 03/03/10): split large black pedigree
## 6. plink2merlin (add 03/03/10): convert plink export pedigree and map file to merlin format
## 7. read.tbl (add 03/05/10): read merlin tbl format result
## 8. simSNP (add 03/07/10): simulate SNP at specific frequency
## 9. summary table function groups (add 05/31/10) ## -- see more in ~/Thesis/script/results_summary_WM_050510.R
## - annotate.snp(x,qc,col=c("HUGO","ROLE"),newcol=c("HUGO","ROLE"))
##             - merge column "HUGO" "ROLE" from qc to x using "SNP" and change the name to "HUGO" "ROLE"
## - map.snp(x,qc) : add MapInfo or POS to CHR POS
## - load.annot(x) : load my default QC file
## - p.format(x) ## format a p-value vector x for pretty display in LaTex
myformat <- function(x,verbose=FALSE,cut=1e-4) {
    a <- character(length(x))
    index <- 1:length(x)
    index <- index[!is.na(x)]
    for (i in index) {
            a[i] <-    ifelse(abs(x[i])<=cut,paste(gsub("e-0{0,2}","x10$^{-",format(x[i],nsmall=3)),"}$",sep=""),format(x[i],nsmall=3))
            if (verbose) print(a[i])
    }
    return(a)
}
## x <- c(0.01,0.001,0.0001,0.00001,1e-17)
## myformat(x)


load.annot<- function(load=FALSE,annot.path="~/Thesis/results/annotation",qcwf="QC-white-012710.csv",qcbf="QC-black-012710.csv") {
if (load) {
    annot.path <- annot.path
    qcw <- read.csv(file.path(annot.path,"QC-white-012710.csv"),header=T,as.is=T)
            cat("loaded qcw\n")
    qcb <- read.csv(file.path(annot.path,"QC-black-012710.csv"),header=T,as.is=T)
            cat("loaded qcb\n")
}
}
## load.annot(load=TRUE)
## load("~/Thesis/results/annotation/qc-012710.RData")
annotate.snp <- function(x,qc,col=c("HUGO","ROLE"),newcol=c("HUGO","ROLE")) {
  cat("Adding",newcol,"of SNP\n")
  names(qc) <- toupper(names(qc))
  col <- toupper(col)
  x[,c(newcol)] <- qc[match(x$SNP,qc$SNP),c(col)]
  return(x)
}


map.snp <- function(x,qc) {
  cat("Adding CHR & POS of SNP\n")
  x[,c("CHR","POS")] <- qc[match(x$SNP,qc$SNP),c(grep("CHR",names(qc)),grep("MapInfo|POS",names(qc)))]
  return(x)
}


read.merlin <- function(pedfile,datfile,force=1,...) {
dat=read.table(file=datfile,header=F,col.names=c("TYPE","NAME"),as.is=T,...)
expectedcol<-sum(dat[,1]=="M")*2+sum(dat[,1]=="S2")*2+sum(dat[,1]=="S")+sum(dat[,1]=="T")+sum(dat[,1]=="C")+5
pedigree=read.table(file=pedfile,header=F,as.is=T,na.strings="x",...)
# Warn if the number of column in pedigree file is different from expected
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


write.merlin <- function(.x,prefix,.ped,.dat,type=NULL) {
  if (missing(prefix)) stop("Please specify output file")
### SITUATION 1: Given readMerlin class data
  if (!missing(.x)) {
        if (class(.x) == "readMerlin") {
          write.table(.x$ped,file=paste(prefix,"ped",sep="."),sep="\t",col.names=F,row.names=F,na="x",quote=F)
          write.table(.x$dat,file=paste(prefix,"dat",sep="."),sep="\t",col.names=F,row.names=F,na="x",quote=F)
        }
  }
### SITUATION 2: given pedigree and dat data frame  
  if (!missing(.ped) & !missing(.dat)) {
          write.table(.ped,file=paste(prefix,"ped",sep="."),sep="\t",col.names=F,row.names=F,na="x",quote=F)
          write.table(.dat,file=paste(prefix,"dat",sep="."),sep="\t",col.names=F,row.names=F,na="x",quote=F)
        } else if (missing(.ped) &  missing(.dat)) {
          stop("Please specify Pedigree and dat data")
        } else  
### SITUATION 3: pedigree data frame available but no dat --- specify type to create dat automatically
  if (!missing(.ped) & missing(.dat) & length(type) > 0 ) {
        cat("generating dat file from pedigree data frame\r")
        length.type <- length(type[type=="M"])*2+length(type[!type=="M"]) # calculate expected length of data name
        if (length.type != ncol(.ped)-5) stop("type has length shorter than pedigree data")
        type.name <- character(length(type))
        colnum <- 6
        allname <- names(.ped)
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
        .dat <- data.frame("type"=type,"name"=type.name)
    ord <- order(.ped[,1],.ped[,2])
    .ped <- .ped[ord,]
        write.table(.ped,file=paste(prefix,"ped",sep="."),sep="\t",col.names=F,row.names=F,na="x",quote=F)
        write.table(.dat,file=paste(prefix,"dat",sep="."),sep="\t",col.names=F,row.names=F,na="x",quote=F)
  } else stop("Please specify type for dat file to write out")
}  
 
freq <- function(a) {
  if (length(a)==3 ) {
        b <- (2*a[3]+a[2])/(2*sum(a))
  }
  if (is.na(a[3])) a[3] <- 0
  if (is.na(a[2])) a[2] <- 0
  b <- (2*a[3]+a[2])/(2*sum(a))
  return(b)
}
snp.freq <- freq <- function(x) {
    a<-table(x)
  if (length(a)==3 ) {
        b <- (2*a[3]+a[2])/(2*sum(a))
  }
  if (is.na(a[3])) a[3] <- 0
  if (is.na(a[2])) a[2] <- 0
  b <- (2*a[3]+a[2])/(2*sum(a))
  return(b)
}


splitwhite <- function(pedigree,newpedfile,newdatfile,verbose=options()$verbose) {
newFam<-function(.x,oldped,newped,id) {
  if (missing(id)) stop("Please enter at least one ID to change family ID")
  if (missing(oldped)) stop("Please enter PED ID to change")
  if (missing(newped)) stop("Please enter new PED ID")
  if (newped %in% .x[,1]) stop("Please check if new pedigree ID is already exist")
  if (verbose) print("Old Pedigree")
  if (verbose)   print(.x[which(.x[,1]==oldped),1:5])
  for (cid in id) {
        if (verbose) print(paste("Switching pedigree for ID",cid,"to",newped))
        if (!cid %in% .x[,2]) warning(paste(cid," not in current pedigree"))
        .x[which(.x[,1]==oldped & .x[,2]==cid),1]<-newped
  }
  if (verbose)   print("New pedigree:")
  if (verbose)   print(.x[which(.x[,1]==newped),1:5])
  return(.x)
}
makeFounder<-function(.x,pedid,id) {
  if (!id %in% .x[which(.x[,1]==pedid),2]) stop(paste(id," not in current pedigree"))
  if (verbose)        print(paste("Make",id,"founder"))
  .x[which(.x[,1]==pedid & .x[,2]==id),3:4] <- 0
  if (verbose)   print(.x[which(.x[,1]==pedid),1:5])
  return(.x)
}
newParents<-function(.x,pedid,id,parid) {
  if (!length(parid)==2) stop("Mising a pair of new parents ID")
  if (!id %in% .x[which(.x[,1]==pedid),2]) stop(paste(id," not in current pedigree"))
  if (verbose)        print(paste("Make Pedigree:",pedid,"ID:",id,"founder"))
  .x[which(.x[,1]==pedid & .x[,2]==id),3:4] <- parid
  if (verbose)   print(.x[which(.x[,1]==pedid),1:5])
  return(.x)
}
dupID<-function(.x,oldped,newped,id) {
  if (missing(oldped)) stop("Please enter current PedID")
  if (missing(newped)) stop("Please enter the new PedID")
  if (missing(id)) stop("Please enter ID to be duplicated")
  if (verbose)        print(paste("Duplicaing",id,"of pedigree",oldped,"to new pedigree:",newped))
  newline<-.x[which(.x[,1]==oldped & .x[,2]==id),]
  newline[,1]<-newped
  .x<-rbind(.x,newline)
  .x<-.x[order(.x[,1],.x[,2]),]
  if (verbose)   print(.x[which(.x[,1]==newped),1:5]) ## replace 1:6 with 1:5 -- pedigree alone no affection
  return(.x)
}
##HALF-SIBS
#if (verbose)  print("Switching Halfsibs family for ID 6260501")
#pedigree[pedigree[,2]==6260501,3]<-62699
#if (verbose)  print("Removed previous father 62697 of 6260501")
#pedigree<-pedigree[-grep(62697,pedigree[,2]),]


#pedigree[pedigree[,2]==25214001,3:4] <- c(25114001,25114012)
# All split pedigree starts with 98*
# Split PEDIGREE 482 as follow to 98482
#**************************** [STR: LINKAGE ANALYSIS DATA] Comment OUT *************************************#
#pedigree<-newFam(pedigree,482,98482,c(4820301,4820312,4820321,4820322,4820323,4820324,4820325,4820326))
#make 4820301 a founder
#pedigree <- makeFounder(pedigree,98482,4820301)
# Split PEDIGREE 1452
#pedigree<-newFam(pedigree,1452,981452,c(14520725,14520712,14520701,82510101,82510201,82510301,82510401))
#pedigree<-makeFounder(pedigree,981452,14520701)                      
# Split PEDIGREE 11521232
#pedigree<- newFam(pedigree,11521232,981232,c(123299,123298,21011501,22412301,23312301,23312312,22412323,22412324,23312321,23312322))
#pedigree<-makeFounder(pedigree,981232,21011501)
## Split PEDIGREE 415
#pedigree <- newFam(pedigree,415,98415,c(4150021,4150221,4150222,4150223,4150001,4150012,4150101,4150201,4150212,41598,41599))
#pedigree <- dupID(pedigree,98415,415,41599)
#**************************** [STR: LINKAGE ANALYSIS DATA] Comment IN *************************************#
if ( dim(pedigree[pedigree[,1]==482,])[[1]] > 0 ){
  pedigree<-newFam(pedigree,482,98482,c(4820301,4820312,4820321,4820322,4820323,4820324,4820325,4820326))
  pedigree <- makeFounder(pedigree,98482,4820301)
}
pedigree<-newFam(pedigree,1452,981452,c(14520725,14520712,14520701))
pedigree<-makeFounder(pedigree,981452,14520701)                      
pedigree<- newFam(pedigree,11521232,981232,c(123299,123298,21011501,22412301,23312301,23312312,22412323,22412324,23312321,23312322))
pedigree<-makeFounder(pedigree,981232,21011501)
pedigree <- newFam(pedigree,415,98415,c(4150021,4150221,4150222,4150223,4150001,4150012,4150101,4150201,4150212,41598,41599))
pedigree <- dupID(pedigree,98415,415,41599)




return(pedigree)
#write.table(pedigree2,file=newpedfile,col.names=F,row.names=F,na="x",quote=F)
#write.table(dat,file=newdatfile,col.names=F,row.names=F,na="x",quote=F)
}


# ----------------------- AA -----------------------#


splitblack <- function(pedigree,newpedfile,newdatfile,verbose=FALSE) {
newFam<-function(.x,oldped,newped,id) {
  if (missing(id)) stop("Please enter at least one ID to change family ID")
  if (missing(oldped)) stop("Please enter PED ID to change")
  if (missing(newped)) stop("Please enter new PED ID")
  if (newped %in% .x[,1]) stop("Please check if new pedigree ID is already exist")
  if (verbose) print("Old Pedigree")
  if (verbose)   print(.x[which(.x[,1]==oldped),1:5])
  for (cid in id) {
        if (verbose)          print(paste("Switching pedigree for ID",cid,"to",newped))
        if (!cid %in% .x[,2]) warning(paste(cid," not in current pedigree"))
        .x[which(.x[,1]==oldped & .x[,2]==cid),1]<-newped
  }
  if (verbose)   print("New pedigree:")
  if (verbose)   print(.x[which(.x[,1]==newped),1:5])
  return(.x)
}
makeFounder<-function(.x,pedid,id) {
  if (!id %in% .x[which(.x[,1]==pedid),2]) stop(paste(id," not in current pedigree"))
  if (verbose)        print(paste("Make",id,"founder"))
  .x[which(.x[,1]==pedid & .x[,2]==id),3:4] <- 0
  if (verbose)   print(.x[which(.x[,1]==pedid),1:5])
  return(.x)
}
newParents<-function(.x,pedid,id,parid) {
  if (!length(parid)==2) stop("Mising a pair of new parents ID")
  if (!id %in% .x[which(.x[,1]==pedid),2]) stop(paste(id," not in current pedigree"))
  if (verbose) print(paste("Make Pedigree:",pedid,"ID:",id,"founder"))
  .x[which(.x[,1]==pedid & .x[,2]==id),3:4] <- parid
  if (verbose)   print(.x[which(.x[,1]==pedid),1:5])
  return(.x)
}
dupID<-function(.x,oldped,newped,id) {
  if (missing(oldped)) stop("Please enter current PedID")
  if (missing(newped)) stop("Please enter the new PedID")
  if (missing(id)) stop("Please enter ID to be duplicated")
  if (verbose)        print(paste("Duplicaing",id,"of pedigree",oldped,"to new pedigree:",newped))
  newline<-.x[which(.x[,1]==oldped & .x[,2]==id),]
  newline[,1]<-newped
  .x<-rbind(.x,newline)
  .x<-.x[order(.x[,1],.x[,2]),]
  if (verbose)   print(.x[which(.x[,1]==newped),1:5])
  return(.x)
}
# TWINS --> removed phenotype of one of the MZ twins
# -- removed contaminated samples 70110201
# -- removed contaminated samples 70420101
# - TBD -
if (verbose)  print("# HALF-SIBS")
if (nrow(pedigree[pedigree[,2]==71120301,]) > 0) {
  cat("Switching halfsibs 71120301 family\n")
  pedigree[pedigree[,2]==71120301,3:4]<-c(711295,711298)
}
if (nrow(pedigree[pedigree[,2]==71120401,]) > 0) {
  cat("Switching halfsibs 71120401 family\n")
  pedigree[pedigree[,2]==72130401,3:4]<-c(721399,721398)
}  
if (nrow(pedigree[pedigree[,2]==82150101,]) > 0) {
  cat("Switching halfsibs 82159191 family\n")
pedigree[pedigree[,2]==82150101,3:4]<-c(821599,821598)
}
if (nrow(pedigree[pedigree[,2]==71830301,]) > 0) {
  cat("Switching halfsibs 71830301 family\n")
pedigree[pedigree[,2]==71830301,3:4]<-c(718399,718398)
}
if (nrow(pedigree[pedigree[,2]==71850621,]) > 0) {
  cat("Switching halfsibs 71850621 family\n")
  pedigree[pedigree[,2]==71850621,3:4]<-c(71850612,71850601)
}
if (nrow(pedigree[pedigree[,2]==70680212,]) > 0) {
  cat("Fix incorrect nuclear family 70680212\n")
  pedigree[pedigree[,2]==70680121,3:4]<-c(70680812,70680801)
}
if (nrow(pedigree[pedigree[,1]==201,]) > 0 ) {
  cat("SPLIT PEDIGREE 201\n")
  pedigree<-newFam(pedigree,201,98201,c(20199,2010001,2010212,2010201,2010214,2010401,2010413,2010221,2010222,2010241,2010242,2010432,2010433))
  # # Probably have to add one recode to duplicate 20198
  pedigree<-dupID(pedigree,201,98201,20198)
}
  if (nrow(pedigree[pedigree[,1]==255,]) > 0 ) {
    # Split PEDIGREE 255
  pedigree<-newFam(pedigree,255,98255,c(2550112,2550101,2550121,2550122,2550123,2550124))
  pedigree<-makeFounder(pedigree,98255,2550101)
  }
## ----------------------- COMMENT OUT FOR LINKAGE -----------------##
## Split PEDIGREE 70827184
##pedigree <- newFam(pedigree,70827184,7184,c(718499,718498,71840101,71840112,71840201,70820522))
## Duplicate ID
##pedigree <- dupID(pedigree,7184,70827184,71840101)
## ----------------------- COMMENT OUT FOR LINKAGE -----------------##
return(pedigree)
#write.table(ped4,file=newpedfile,col.names=F,row.names=F,na="x",quote=F)
#write.table(dat,file=newdatfile,col.names=F,row.names=F,na="x",quote=F)
}


# Convert PLINK Export in pedigree format to Melrin format with or without phenotype file to merge with
# Phenotype file hsould have the first two column as PEDID and PID following by PHENOTYPE DATA
# Bhoom Suktitipat
# Last update: 03/03/10
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


## FUNCTION TO READ TBL FILE OUTPUT FROM MERLIN
## COMBINE THE RESULTS FROM ALL CHROMOSOME TO ONE FILE
## SPLIT RESULTS BY PHENOTYPE or
## COMBINE THE RESULTS SIDE BY SIDE with only
## P-VALUE or/and H2
read.tbl <- function(prefix="vrfsplwhitefixed", # prefix of the input file with -assoc-chr
                         filter, # only output results this p-value
                         all=missing(filter),
                         map="white-verify-merlin.map",
                         type="merlin",
                         wide=FALSE,
                         chr=c(paste("0",1:9,sep=""),10:22),
                         ...
                         ) {
  if (type == "merlin") {
        map <- read.table(map,header=F,col.names=c("CHR","SNP","CM"),as.is=T)
        result <- list()
        for (i in 1:length(chr)) {
          result[[i]] <- read.table(file=paste(prefix,"-assoc-chr",chr[i],".tbl",sep=""),header=T,as.is=T,na.strings="-")
        }
        resultTable <- do.call(rbind,result)
        resultTable[,"POS"] <- map[match(resultTable$SNP,map$SNP),"CM"]*1000000
        rm(result)
        if (all) {
          print(paste("Writing Results for All phenotype",paste(prefix,"-ALL.csv",sep="")))
          write.csv(resultTable,file=paste(prefix,"-ALL.csv",sep=""),row.names=F,quote=F,na="")
          phenotypes <- unique(resultTable$TRAIT)
          for (phenotype in phenotypes) {
            print(paste("Writing Results for",phenotype,"in",paste(prefix,"-",phenotype,".csv",sep="")))
            write.csv(resultTable[which(resultTable$TRAIT==phenotype),],file=paste(prefix,"-",phenotype,".csv",sep=""),
                      row.names=F,quote=F,na="")
                    }
        } else {
          outfile <- paste(prefix,"-",filter,"-ALL.csv",sep="")
          print(paste("Writing filtered Results filtered at least",
                      filter,"for All phenotype",outfile))
          filter.result <- resultTable[which(resultTable$PVALUE <= filter),]
          filter.result <- filter.result[order(paste(filter.result$TRAIT,filter.result$PVALUE,sep="")),]
          write.csv(filter.result,file=paste(prefix,"-",filter,".csv",sep=""),row.names=F,quote=F,na="")
          phenotypes <- unique(resultTable$TRAIT)
          for (phenotype in phenotypes) {
            outfile <- paste(prefix,"-",filter,"-",phenotype,".csv",sep="")
            print(paste("Writing Results filtered at least",filter,
                        "for",phenotype,"in",outfile))
            write.csv(filter.result[which(filter.result$TRAIT==phenotype),],
                      file=outfile,row.names=F,quote=F,na="")
          }
        }
    
  } else {
        cat("To be implemented")
  }
  return(resultTable)
}


## -------------------------Example of read.tbl ----------------------------------#
if (0) {
  source("~/Thesis/script/readpedigree.R")
  setwd("~/Thesis/results/gee-results-eig2/verifyGEE-cluster/merlin64/splitPed021910/rst/cluster/fixed")
  white.rst <- read.tbl(prefix="vrfsplwhitefixed",map="white-verify-merlin.map")
  setwd("..")
  black.rst <- read.tbl(prefix="vrfsplblack",map="black-verify-merlin.map")
  setwd("~/Thesis/results/candidateGene")
  white.cand <- read.tbl("cand-white",chr=c("01","02","03","04","06","09","11","12","17","19","22"),
                             map="cand-white.map")
}
## -------------------------End Example read.tbl----------------------------------#
patID <- function(ped,id) {
  patid <- ped[which(ped[,2]==id),3]
  return(patid)
}
matID <- function(ped,id) {
  matid <- ped[which(ped[,2]==id),4]
  return(matid)
}
allele.drop <- function(ped,patid,matid,prop.pat,prop.mat,allele1.name,allele2.name) {
  if (prop.pat <0.5) {
        a1 <- ped[which(ped[,2]==patid),allele1.name]
  } else {
        a1 <- ped[which(ped[,2]==patid),allele2.name]
  }
  if (prop.mat <0.5) {
        a2 <- ped[which(ped[,2]==matid),allele1.name]
  } else {
        a2 <- ped[which(ped[,2]==matid),allele2.name]
  }
  alleles <- c(a1,a2)
  return(alleles)
}


simSNP <- function(
                       ped=ped,
                       MAF=0.01,
                       markerName="SNP1",
                       allele=c(1,3),
                       recodeA=FALSE,
                       verbose=FALSE
                       ) {
  .ped <- ped[,1:5]
  if (missing(allele)) {
        if (runif(1)<0.5) {
          if (runif(1)<0.5) {
            allele=c("A","G")
          } else
          {
            allele=c("G","A")
          }
        }
        else {
          if (runif(1)<0.5) {
            allele=c("C","T")
          } else {
            allele=c("T","C")
          }
        }
  }
  markers <- list()
  for (j in 1:length(MAF)) {
        maf <- MAF[j]
        markername <- markerName[j]
        cat("Simulating",markername,"\n")
  allele1 <- paste(markername,"_1",sep="")
  allele2 <- paste(markername,"_2",sep="")
  row.names(.ped) <- .ped[,2]
  founder <- .ped[which(.ped[,3]==0 & .ped[,4]==0),2]
  .ped[which(.ped[,3]==0 & .ped[,4]==0),allele1] <- (runif(length(founder))<maf)+0
  .ped[which(.ped[,3]==0 & .ped[,4]==0),allele2] <- (runif(length(founder))<maf)+0
  others <- .ped[which(!.ped[,3]==0 & !.ped[,4]==0),2]
  family <- unique(.ped[,1])
  pat <- .ped[which(!.ped[,3]==0 & !.ped[,4]==0),3]
  mat <- .ped[which(!.ped[,3]==0 & !.ped[,4]==0),4]
  geta1 <- runif(length(others)) # choice for father allele
  geta2 <- runif(length(others)) # choice for mother allele
  for (i in 1:length(others)) {
        alleles <- allele.drop(.ped,pat[i],mat[i],geta1[i],geta2[i],allele1,allele2)
        a1 <- alleles[1]
        a2 <- alleles[2]
        .ped[which(.ped[,2]==others[i]),allele1] <- a1
        .ped[which(.ped[,2]==others[i]),allele2] <- a2
if (verbose) {
        if (!i%%100==0) {
          cat(".")
        } else {
          cat(i,"\n")
        }   }
  }
        if (verbose) cat(i,"\n")
        ## PICK UP Unassign genotype
        k=1
        maxit=3 ## TO PREVENT UNLIMITED LOOP
        others <- unique(c(.ped[,2][is.na(.ped[,allele1])],.ped[,2][is.na(.ped[,allele2])]))
        while (length(others)>0 | k >maxit) {
          geta1 <- runif(length(others)) # choice for father allele
          geta2 <- runif(length(others)) # choice for mother allele
          if (verbose) print(paste(length(others),"people still have missing genotypes: iteration",k))
          for (i in 1:length(others)) {
            id <- others[i]
#            debug(allele.drop)
            alleles <- allele.drop(.ped,patID(.ped,id),matID(.ped,id),geta1[i],geta2[i],allele1,allele2)
            a1 <- alleles[1]
            a2 <- alleles[2]
        .ped[which(.ped[,2]==others[i]),allele1] <- a1
        .ped[which(.ped[,2]==others[i]),allele2] <- a2
          }
        others <- unique(c(.ped[,2][is.na(.ped[,allele1])],.ped[,2][is.na(.ped[,allele2])]))
        k <- k+1
          if (verbose) {
            if (!i%%100==0) {
              cat(".")
            } else {
              cat(i,"\n")
            }   }
        }
        if (length(others) >0 & k > maxit) stop("Genotype still missing;need to increase number of iteration")
        if (verbose) cat(i,"\n")
        if (recodeA) {
          markers[[markername]] <- .ped[,allele1]+.ped[,allele2]
          .ped[,allele1] <- NULL
          .ped[,allele2] <- NULL
        } else {
          .ped[which(.ped[,allele1]==0),allele1] <- allele[1]
          .ped[which(.ped[,allele1]==1),allele1] <- allele[2]
          .ped[which(.ped[,allele2]==0),allele2] <- allele[1]
          .ped[which(.ped[,allele2]==1),allele2] <- allele[2]
          index <- grep(markerName,names(.ped))
          markers[[markername]] <- .ped[,index]
          .ped[,allele1] <- NULL
          .ped[,allele2] <- NULL
        }
  } # loop for generating several SNPs
#  if (recodeA) {
  allnames <- paste(rep(names(markers),each=2),c("_1","_2"),sep="")
  oldname.ped <- names(ped)
  ped <- cbind(ped,do.call(cbind,markers))
  names(ped) <- c(oldname.ped,allnames)
#  } else {
#        ped <- cbind(.ped,do.call(cbind,marker))
#  }
  return(ped)
}


if (0) {
  setwd("~/Thesis/script")
  source("~/Thesis/script/readpedigree.R")
  .ped <- read.table("~/cluster/whiteR21.fam",col.names=c("PID","IID","PAT","MAT","SEX","AFF"),as.is=T)
  pedlist <- list()
  freq <- rep(c(0.01,0.05,0.1,0.2,0.3,0.4),each=10)
        system.time({
  for (i in 1:length(freq)) {
        newped <- simSNP(.ped,MAF=freq[i],markerName="SNP",allele=c(1,3),verbose=TRUE)
        names(newped)[7:8] <- paste("SNP",i,c("_1","_2"),sep="")
        pedlist[[i]] <- newped[,7:8]
        rm(newped)
  }
  })
  new.ped <- cbind(.ped,do.call(cbind,pedlist))
 
  write.table(new.ped,"simwhite.ped",quote=F,row.names=F)




  recodeA <- function(x,y) {
  geno <- x+y
  return(geno)
}


}




genPheno <- function(
                         ped=sim.additive,
                         maf=0.05,
                         geno=7, # if blank please specify the genotype name
                         genotype="SNP1_A", # if black please specify the genotype column
                         pheno="T1", #
                         traitmean=0, # trait mean
                         totalVar=NA,
                         betaG=NA, # need if totalVariance missing
                         locus.h2=0.1,
                         total.h2=0.4,
                         absence=rep(FALSE,nrow(ped)),
                         returnped=FALSE
                         ) {
  # For family j, individual i
  # Y_ij = Beta0 + alpha_j + Beta1_i*Geno_i + e_ij
  # alpha_j is the random family effect ~ N(fam_mean,total.h2-locus.h2)
  # Error ~ N(0,rand_var)
  # BetaG = Sqrt(loc.h2*tot.h2)/2p
  # CHECK IF ONLY betaG or totalVariance was specified
  if (is.na(betaG) & is.na(totalVar) ) stop("Please specify betaG or totalVar")
  if (missing(locus.h2)) stop("Please specify locus specific heritability")
  fam.mean <- function(ped,famMean=0,famVar=famVar) {
        counts <- tapply(ped[,2],ped[,1],length)
        alpha.j <- list()
        sd <- sqrt(famVar)
        for ( i in 1:length(counts) ) {
          alpha.j[[i]] <- rnorm(counts[i],famMean,sd)
          loc <- sum(counts[1:i])
        }
   alpha.j <- do.call(c,alpha.j)
        return(alpha.j)
  }
  if (missing(geno) ) {
        if (missing(genotype)) {
          stop("Please specify either the column of linked genotype or the genotype name")
          }
        else {
          geno <- grep(genotype,names(ped))
        }
  }
  a <- table(ped[which(ped[,3]==0&ped[,4]==0),geno])
  if (missing(maf)) maf <- freq(a)
#  if (p==1) p=0 # Take care of the problem if p == 1
  if (pheno %in% names(ped) ) warn("Phenotype already exist. Phenotype in pedigree will be replaced")
  ped[,pheno] <- NA
  ### CALCULATE variance(locus) = 2pq(beta^2)
  if (is.na(betaG) & !is.na(totalVar) & locus.h2>0) betaG <- sqrt((locus.h2*totalVar)/(2*maf*(1-maf))) # calculate betaG (see book by David Siegmund & Benjamin Yakir: The Statistics of Gene Mapping page 40)
  if (locus.h2==0) betaG=0
  genetic.effect <- betaG*ped[,geno]
  if (is.na(totalVar) ) {
        G <- var(genetic.effect)
        totalVar <- G
  }
  famVar=(total.h2-locus.h2)*totalVar
#  cat("Total Family Variance=",famVar,"Total Heritability=",total.h2,"\n")
  alpha <- fam.mean(ped,famMean=0,famVar)
  error <- rnorm(nrow(ped),0,sqrt(totalVar*(1-total.h2)))
  trait <- traitmean + alpha + genetic.effect + error
  trait[absence] <- NA
  return(trait)
  if (returnped) {
        ped[,pheno] <- trait
        return(ped)
  }
}
annotate.snp <- function(x,qc) {
  cat("Adding Gene name and possible function of SNP\n")
  x[,c("HUGO","ROLE")] <- qc[match(x$SNP,qc$SNP),c("HUGO","ROLE")]
  return(x)
}
map.snp <- function(x,qc) {
  cat("Adding CHR & POS of SNP\n")
  x[,c("CHR","POS")] <- qc[match(x$SNP,qc$SNP),c(grep("MapInfo|POS",names(qc)),grep("CHR",names(qc)))]
}
if (0) {
source("~/Thesis/script/readpedigree.R")


## Import files from UCSC genome browser database
## header based on this file http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/knownGene.sql
hg18 <- read.table("~/Thesis/results/annotation/galaxy/hg18/knownGene.txt",sep="\t",header=F,col.names=c("name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","proteinID","alignID"),as.is=T) ## from http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/knownGene.txt.gz
## combined with kgAlias
## alias <- read.table("~/Thesis/results/annotation/galaxy/hg18/kgAlias.txt",sep="\t",col.names=c("kgID","alias"),as.is=T)
##known.gene <- read.table("~/Thesis/results/annotation/glist-hg18",header=F,col.names=c("CHR","START","STOP","GENE"),as.is=T)
## hugo <- read.table("~/Thesis/results/annotation/galaxy/hg18/HUGO_update060610.csv",sep=",",header=T,as.is=T)
#library(XML)
#hugo <- xpathApply( htmlTreeParse("~/Thesis/results/annotation/galaxy/hg18/HUGO_download_update060610.html", useInt=T), "//td", function(x) xmlValue(x))


library(foreign)
kgxref <- read.table("~/Thesis/results/annotation/galaxy/hg18/kgXref.txt",header=F,col.names=c("kgID","mRNA","spID","spDisplayID","geneSymbol","refseq","protAcc","description"),as.is=T,sep="\t")
kgxref <- read.dta("~/Thesis/results/annotation/galaxy/hg18/kgxref.dta")#,col.names=c("kgID","mRNA","spID","spDisplayID","geneSymbol","refseq","protAcc","description"))
load("~/Thesis/results/annotation/qc-012710.RData")
names(kgxref) <- c("kgID","mRNA","spID","spDisplayID","geneSymbol","refseq","protAcc","description")
kgXref <- kgxref[grep("^N",kgxref$mRNA),] ## KEEP ONLY GENE WITH REFSEQ
hg18[,"HUGO"] <- kgXref[match(hg18$name,kgXref$kgID),"geneSymbol"]
hg18.unique <- hg18[order(hg18$HUGO,hg18$exonCount,decreasing=T),]
dupID <- duplicated(hg18.unique$HUGO)
hg18.unique <- hg18.unique[!dupID,]
hg18.unique <- hg18.unique[order(hg18.unique$chrom,hg18.unique$txStart),]
gene <- unique(qc$HUGO)
sum(gene %in% hg18.unique$HUGO)
hg18.unique$chr <- gsub("chr|_random|_h2_hap1|_cox_hap1|_qbl_hap2","",hg18.unique$chrom)
hg18.unique$chr <- gsub("X",23,hg18.unique$chr)
hg18.unique$chr <- gsub("Y",24,hg18.unique$chr)
if (0) write.table(hg18.unique,file="hg18-glist",sep=":",row.names=F,col.names=T,na="")
genelist <- hg18.unique
row.names(genelist) <- 1:nrow(genelist)
#qc.genes <- unique(qc$HUGO)
#known.gene <- known.gene[order(known.gene$CHR,known.gene$START),]
if (0) system.time(save.image("~/Thesis/results/annotation/closestgene/genelist_QC-environment.RData",compress=T))
}


closest.gene <- function(snp,qc=qc,genelist=genelist,verbose=FALSE) {
  ## snp = query SNP
  ## qc = database = data frame with at least these columns
  ## CHR, SNP, POS|BP, HUGO|GENE, ROLE
  ## genelist = data frame with :
  qcnames <- toupper(names(qc))
  role.col <- grep("ROLE",qcnames)
  hugo.col <- grep("HUGO|GENE",qcnames)
  pos.col <- grep("POS|BP",qcnames)
  snp.col <- grep("SNP",qcnames)
  chr.col <- grep("CHR",qcnames)
  ##    snp = "rs12063900"
  ## Get SNP info
  snp.pos <- grep(snp,qc[,snp.col])
  chr <- qc[snp.pos,chr.col]
  pos <- qc[snp.pos,pos.col]
  gene.loc <- genelist[which(genelist$chr==chr & genelist$txStart<pos & genelist$txEnd>pos),]
  if ( nrow(gene.loc) ==1) {
        gene <- gene.loc[,"HUGO"]
        exStart <- unlist(strsplit(gene.loc[,"exonStarts"],","))
        exEnd <- unlist(strsplit(gene.loc[,"exonEnds"],","))
        exon <- 1:length(exStart)
        e.start = max(exon[exStart <= pos])
        e.end = min(exon[exEnd >= pos])
        if (e.start == e.end) {
          print(paste(snp,"is in exon",e.start,"of",gene))
          location <- paste("exon",e.start,sep="")
        } else if (e.start < e.end) {
          print(paste(snp,"is in intron",e.start,"of",gene))
          location <- paste("intron",e.start,sep="")
        }
  } else {
        location <- "near gene"
        gene <- ""
  } ## END IF WHEN NO GENE WAS FOUND
  ## genelist[which(genelist$chr == chr & (genelist$txStart <= pos & pos <= genelist$txEnd)),])
  disStart <- pos-genelist$txStart[genelist$chr == chr & genelist$txStart<=pos]
  ## disEnd <- pos-genelist$txEnd[genelist$chr == chr & genelist$txStart <= pos]
  downstream <- (genelist[genelist$chr==chr,][which.min(disStart),c(2,3,4,5,13)])
  upstream <- (genelist[genelist$chr==chr,][which.min(disStart)+1,c(2,3,4,5,13)])
  dd <- pos-downstream[4]
  ud <- upstream[3]-pos
  if (dd < ud) {
        if (options()$verbose) {
          print(downstream)
          print(paste("Upstream from ",upstream[5],"=",upstream[3]-pos,"bp"))
        }
        print(paste("Downstream from ",downstream[5],"=",pos-downstream[4],"bp"))
  } else {
        if (options()$verbose) {
          print(upstream)
          print(paste("Downstream from ",downstream[5],"=",pos-downstream[4],"bp"))
        }
        print(paste("Upstream from ",upstream[5],"=",upstream[3]-pos,"bp"))
        ##genelist[genelist$chr==1,][which.min(disEnd),c(2,3,4,5,13)]
  }
  output <- list(chr=chr,pos=pos,snp=snp,gene=gene,location=location,upstream=upstream,downstream=downstream)
  return(output)
}


if (0) {
  closest.gene("rs954213",verbose=T,qc=qc,genelist=genelist)
        closest.gene("rs1383089",verbose=T,qc=qc,genelist=genelist)
  genelist[12675,]
## CRP signals
closest.gene("rs12063900",qc=qc,genelist=genelist)
closest.gene("rs2167529",qc=qc,genelist=genelist)
closest.gene("rs13098723",qc=qc,genelist=genelist)
closest.gene("rs7695413",qc=qc,genelist=genelist)
closest.gene("rs3864148",qc=qc,genelist=genelist)
closest.gene("rs9444762",qc=qc,genelist=genelist)
closest.gene("rs9506242",qc=qc,genelist=genelist)
closest.gene("rs7165083",qc=qc,genelist=genelist)
## IL-6 signals
ril6 <- sapply(c("rs954213","rs4428953","rs1383089","rs957215",
                     "rs7719405","rs1984703","rs1584535","rs11563579","rs10110604"),
                   closest.gene,qc=qc,genelist=genelist)
## MCP-1 signals
rmcp1 <- sapply(c("rs6689939","rs1885011","rs6428744","rs2884802",
                      "rs9999992","rs2263190","rs9519286","rs6068934"),
                    closest.gene,qc=qc,genelist=genelist)
}


read.mach <- function(genfile=mldose,mapfile=mlinfo, pedfile=pedfile, sep.gen=" ",sep.ped=",") {
  ## obsolete updated with read.mldose
  print("Reading in Data")
  ped.dat <- read.table(genfile, header = FALSE, sep = sep.gen, as.is=T)
  ped.dat[,1] <- matrix(unlist(strsplit(ped.dat[,1],"->")),byrow=T,ncol=2)[,2]
  ped.dat[,2] <- NULL
  mlinfo <- read.table(mapfile,header=T,as.is=T)
  snp.names <- paste(gsub("X","",mlinfo[,1]),mlinfo[,2],sep="_")
  names(ped.dat) <- c("id",snp.names)
  pedigree <- read.table(pedfile, header = TRUE, sep = sep.ped)
  system.time(pedigree[,names(ped.dat)[-1]] <- ped.dat[match(pedigree[,2], ped.dat[,1]),-1])
  print("Done reading in data")
  data = list(snps = snp.names, geno = pedigree)
  class(data)="mach.data"
  return(data)
}


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


read.washu.mldose <- function(mldosefile,mlinfofile,pedfile,sep.ped=",",
                                sep.mldose=" ",sep.mlinfo=" ") {
  message("Read WashU-modified mldose file for AA imputed Data\nUpdated 11/29/10")
  ped.dat <- read.table(mldosefile, header = FALSE, sep = sep.mldose, as.is=T)
  message("Reading mlinfo files")
  mlinfo <- read.table(mlinfofile,header=FALSE,col.names=c("SNP","AL1","AL2","FREQ1","MAF","QUALITY","RSQ"),as.is=T,sep=sep.mlinfo)
  snp.names <- paste(gsub("X","",mlinfo[,1]),mlinfo[,2],sep="_")
  names(ped.dat) <- c("id",snp.names)
  pedigree <- read.table(pedfile, header = TRUE, sep = sep.ped,col.names=c("fid","id","pat","mat","sex"))
  genotype <- merge(pedigree, ped.dat, by = "id")
  names(genotype)[1:5] <- c("IID","FID","PAT","MAT","SEX")
  print("Done reading in imputed genotype data")
##  return(list(genotype = genotype, snps = snp.names, ped = ped))
  return(genotype)
}  


## http://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session
.ls.objects <- function (pos = 1, pattern, order.by,
                            decreasing=FALSE, head=FALSE, n=5) {
        napply <- function(names, fn) sapply(names, function(x)
                                             fn(get(x, pos = pos)))
        names <- ls(pos = pos, pattern = pattern)
        obj.class <- napply(names, function(x) as.character(class(x))[1])
        obj.mode <- napply(names, mode)
        obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
        obj.size <- napply(names, object.size)
        obj.dim <- t(napply(names, function(x)
                            as.numeric(dim(x))[1:2]))
        vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
        obj.dim[vec, 1] <- napply(names, length)[vec]
        out <- data.frame(obj.type, obj.size, obj.dim)
        names(out) <- c("Type", "Size", "Rows", "Columns")
        if (!missing(order.by))
            out <- out[order(out[[order.by]], decreasing=decreasing), ]
        if (head)
            out <- head(out, n)
        out
}
lsos <- function(..., n=10) {
        .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}
