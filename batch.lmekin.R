message("batch.lmekin() wrapper of lmekin in kinship library")
message("AUTHOR: BHOOM SUKTITIPAT")
message("Initial Version Date: 12/06/2010")
message("Update: 09/27/2011")
#12/28/2010 Fix: makekinship function to be more reliable on FID,IID,PAT,MAT structure of PLINK
#09/27/2011 Fix: assuming that the missing chr = 1 [ line 40-42 ]
debug=0
if (debug) {
  setwd("~/Thesis/script/debug_data")
  ## IMPUTE DATA
  imputedgeno <- read.washu.mldose("mach.chr19.7.out.mldose.gz","mach.chr19.7.out.mlinfo.gz","blackped.csv")
  ## GENOTYPE DATA
  typedgeno <- read.table("black.geno.raw",sep=" ",header=TRUE,as.is=TRUE)
  ## SURVIVAL DATA
  time2event <- read.table("chd_event.csv",sep=",",header=TRUE,as.is=TRUE)
  genodata=genotype;phendata=phenotype;k=10;genostart=6;genoend=100;outcomes="IncidentAllCHD" ## genostart=6 for imputed data
  covar="+AGE+SEX+V1+V2"
  suffix="cox.chd1p"
  outdir="COGENT/COXCHD"
  chr=19
  time="FUTIME"
  outcome="IncidentAllCHD"

  ## CONTINUOUS PHENOTYPE
  cont <- read.table("caccogent.csv",header=TRUE,as.is=TRUE,sep=",")
  genodata=typedgeno;phendata=cont;k=10;genostart=7;genoend=100;outcome="LOGCAC" ## genostart=7 for genotype data
  outdir="test_lmekin"
  suffix="lmekin"
  chr=5
  covar="+AGE + SEX + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10"
  famid="FID"
}
## EXAMPLE
##
##batch.lmekin(genodata=typedgeno,phendata=cont,k=100,genostart=7,genoend=100,outcomes="LOGCAC1",covar="+AGE+SEX+V1+V2",suffix="lmekin",outdir="result")
batch.lmekin <- function(genodata=genotype,phendata=phenotype,k=1000,datatype,
                       genoend=ncol(genodata),outcomes,famid="FID",
                       covar="",suffix="",outdir="result",chr=chr,genostart
                       ,...) {
  if (missing(chr) ) {
    stop("Missing chr option")
    #chr=1
  }
  ##combine phenotype with genotype
  if ( missing(genostart)) {
    if ( datatype=="genotyped") {
      message("Analyzing Genotype Data")
      genostart=7
    } else if( datatype=="imputed") {
      message("Analyzing Imputed Data")
      genostart=6
    }
  }
  system(paste("mkdir -p",outdir))
  require(kinship)
  totalsnp <- genoend-genostart+1
  pedhead <- names(genodata)[1:5]
  message("Calculate KINSHIP coefficient matrix using kinship()")
  kinmat <- makekinship(famid=genodata[,"FID"],id=genodata[,"IID"],father.id=genodata[,"PAT"],mother.id=genodata[,"MAT"])
  ## example of covar "+AGE + SEX + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10"
  cv <- gsub("\\s+","",covar,perl=T) ## remove extra white spaces
  cv <- gsub("^\\+","",cv,perl=T) ## remove the beginning + sign
  covarList <- unlist(strsplit(cv,"\\+"))
  for (outcome in outcomes) {
    ## Creating Output file
    outfilename <- paste(outdir,paste(paste(outcome,suffix,sep="."),chr,"txt",sep="."),sep="/")
    ##             coef exp(coef)  se(coef) robust se      z Pr(>|z|)  
    cat("SNP,coef,se,t,Pr(>|t|),NFAM,FREQ,Memo\n",sep="",file=outfilename)
    p.lme <- data.frame(SNP=character(0),coef=double(0),se=double(0),
                        t=double(0),p=double(0),NCluster=double(0),FREQ=character(0),Memo=double(0),stringsAsFactors=F)
    ## Split the analysis into pieces
    begin <- genostart
    last <- max(min(k+genostart-1,genoend),genostart)
    allpcs <- ceiling(totalsnp/k)
    pcs <- 1
    while (begin <= genoend) {
      snp <- names(genodata)[begin:last]
      ## Check if Phendata contains "ID/id" variable
      phenid.col <- grep("^id$",names(phendata),ignore.case=TRUE)[1]
      if (is.na(phenid.col) ) {
        cat("Variable name 'ID or id' not found in phenotype data. Please check\n")
        stop
      } else {
        names(phendata)[phenid.col] <- "ID"
      }
      ## Check to see if phenotype file contain family ID "family" "FID" or "famid"
      fid.col <- grep("^fid$|^family$|^famid$",names(phendata),ignore.case=TRUE)[1]
      if ( is.na(fid.col) ) {
        rem.col <- phenid.col
      } else {
        rem.col <- c(phenid.col,fid.col)*(-1)
        names(phendata)[fid.col] <- "FID"
      }
      data <- genodata[,c(pedhead,snp)]
      data[,names(phendata)[rem.col]] <- phendata[match(data$IID,phendata[,phenid.col]),rem.col]
      allsnp=length(snp)
      for (j in 1:allsnp) {
        refpoint <- begin-genostart
        ## Identify SNP for testing
        mysnp <- snp[j]
        message(paste("Analyzing",mysnp))
        if (cv == "") {
          datfm <- as.formula(paste(outcome,"~",mysnp,"+ IID +",famid))
        }
        else {
          datfm <- as.formula(paste(outcome,"~",mysnp,"+",paste(covarList,collapse="+"),"+ IID +",famid))
        }
        gwdat <- model.frame(datfm,data=data,na.action=na.omit)
        gwdat <- gwdat[order(gwdat[,famid]),]
        af <- table(round(gwdat[,2]))
        NClust <- length(unique(gwdat[,famid]))
        snpcall <- round(gwdat[,2])
        freq <-  paste(paste(names(af),af,sep=":"),collapse="-")
        ##        cellcount <- table(gwdat[,1],snpcall)
        ## Output of mean allele frequency in case/control
        CASE=paste(dimnames(af)[[1]],af,sep=":",collapse="|")
        ## if both cases and controls exist, do the analysis.
        ##if (sum(cellcount[1,]!=0)>1 & sum(cellcount[2,]!=0)>1) {
        if (! is.na(var(af))) {
          message(paste("Analyzing Piece",pcs,"snp",((pcs-1)*k)+j,":",mysnp,"of",genoend-genostart))
          if (! cv == "") {
            mf <- as.formula(paste(outcome,"~",mysnp,"+",paste(covarList,collapse="+")))
          }
          else {
            mf <- as.formula(paste(outcome,"~",mysnp))
          }
          assign("id",gwdat$IID,inherit=TRUE,envir=.GlobalEnv)
          sLME <- try(lmekin(mf ,random= ~ 1|id,data=gwdat,varlist=list(kinmat)))
          if ("try-error" %in% class(sLME)) {
            ## LMEKIN failed as well.
            message(paste(mysnp,"LMEME failed.\n"))
            rst <- c(mysnp,NA,NA,NA,NA,NClust,CASE,1)
            rst <- data.frame(t(rst))
            names(rst) <- names(p.lme)
            p.lme <- rbind(p.lme,rst)
          } else {
            rLME <- sLME$ctable[mysnp,]
            if (options()$verbose) print(sLME)
            rst <- c(mysnp,signif(rLME[1],3),signif(rLME[2],3),signif(rLME[3],4),signif(rLME[4],4),NClust,CASE,0)
            rst <- data.frame(t(rst))
            names(rst) <- names(p.lme)
            p.lme <- rbind(p.lme,rst)
          }
        } 
        ## when outcome only exist in one category for the marker
        else {
          message(paste("No variation of outcome for SNP",mysnp,"for",outcome))
          rst <- c(mysnp,NA,NA,NA,NA,NClust,CASE,2)
          rst <- data.frame(t(rst))
          names(rst) <- names(p.lme)
          p.lme <- rbind(p.lme,rst)
        }
      } ## Loop through all SNPs in piece k_th
      message(paste(outcome,":Chr",chr,":Pc",pcs,"of",allpcs,"SNP",mysnp,":",refpoint+j,"of",genoend-genoend+1,"analyzed!\n",sep=" "))
      write.table(p.lme,sep=",",col.names=F,file=outfilename,row.names=F,quote=F,append=T)
      rm(data,fid.col,allsnp,refpoint,mysnp,datfm,gwdat,af,NClust,snpcall,CASE,sLME,rLME,p.lme,rst)
      p.lme <- data.frame(SNP=character(0),coef=double(0),se=double(0),
                          t=double(0),p=double(0),NCluster=double(0),FREQ=character(0),Memo=double(0),stringsAsFactors=F)
      gc()
      begin <- last+1
      last <- min(begin+k-1,genoend)
      pcs <- pcs+1
      if (options()$verbose) cat("Repeat while loop\n")
    } # while loop speed up data access
    outfilename <- paste(outdir,paste(paste(outcome,suffix,sep="."),chr,"txt",sep="."),sep="/")
  }# for each outcome
}

