summary.machcoxme <- function(x) {
  if (!class(x)=="coxme") stop("X is not a coxme object")
  beta <- x$coefficients
  nvar <- length(beta$fixed)
  nfrail <- nrow(x$var) - nvar
  se <- sqrt(diag(x$var)[nfrail + 1:nvar])
  xout <- list(coef=beta$fixed[1],se=se[1])
  return(xout)
}

## Memo code for mach.GeeLogit
## 1: COXME Error
## 0: COXME Success

debug=0
if (debug) {
  setwd("debug_data")
  genotype <- read.washu.mldose("mach.chr19.7.out.mldose.gz","mach.chr19.7.out.mlinfo.gz","blackped.csv")
  phenotype <- read.table("chd_event.csv",sep=",",header=TRUE,as.is=TRUE)
  genodata=genotype;phendata=phenotype;k=10;genostart=6;genoend=100;outcomes="IncidentAllCHD"
  covar="+AGE + SEX + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10"
  covar="+AGE+SEX+V1+V2"
  suffix="cox.chd1p"
  outdir="COGENT/COXCHD"
  chr=19
  time="FUTIME"
  outcome="IncidentAllCHD"
  famid="FID"
}

mach.coxme <- function(genodata=genotype,phendata=phenotype,k=1000,genostart=6,
                       genoend=ncol(genodata),outcomes=names(phenotype)[-2:-1],famid="FID",time,
                       covar="",suffix="",outdir="./result/white",mincount=5,chr=1
                       ,...) {
  ##combine phenotype with genotype
  message("mach.coxme hard code to adjusted for AGE SEX and V1-V10, 11/29/10")
  system(paste("mkdir -p",outdir))
  require(kinship)
  totalsnp <- genoend-genostart+1
  pedhead <- names(genodata)[1:5]
  ## example of covar "+AGE + SEX + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10"
  cv <- gsub("\\s+","",covar,perl=T) ## remove extra white spaces
  cv <- gsub("^\\+","",cv,perl=T) ## remove the beginning + sign
  covarList <- unlist(strsplit(cv,"\\+"))
  covar.string <- ifelse(length(covarList)==0,"",paste("+",paste(covarList,collapse="+"),sep=""))
  ##    cv <- gsub("\\+","-",cv,perl=T) ## example of cv out ""AGE-SEX-V1-V2-V3-V4-V5-V6-V7-V8-V9-V10"
  for (outcome in outcomes) {
    ## Creating Output file
    outfilename <- paste(outdir,paste(paste(outcome,suffix,sep="."),chr,"txt",sep="."),sep="/")
    cat("SNP,coef,ExpCoef,SE,Z,Pr(>|z|),Cluster,FREQ,Memo\n",sep="",file=outfilename)
    p.cox <- data.frame(SNP=character(0),coef=double(0),expCoef=double(0),se=double(0),
                        z=double(0),p=double(0),NCluster=double(0),FREQ=character(0),Memo=double(0),stringsAsFactors=F)
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
        rem.col <- (-1)*phenid.col
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
        datfm <- as.formula(paste(outcome,"~",mysnp,covar.string,"+",time,"+",famid))
        gwdat <- model.frame(datfm,data=data,na.action=na.omit)
        gwdat[,"censoredOutcome"] <- Surv(gwdat[,time],gwdat[,outcome])
        gwdat <- gwdat[order(gwdat[,famid]),]
        af <- round(tapply(round(gwdat[,2]),gwdat[,1],mean,na.rm=TRUE),2)
        NClust <- length(unique(gwdat[,famid]))
        snpcall <- round(gwdat[,2])
        ## freq <-  paste(paste(names(af),af,sep=":"),collapse="-")
        cellcount <- table(gwdat[,1],snpcall)
        ## Output of mean allele frequency in case/control
        CASE=paste(dimnames(af)[[1]],af,sep=":",collapse="|")
        ## if both cases and controls exist, do the analysis.
        if (sum(cellcount[1,]!=0)>1 & sum(cellcount[2,]!=0)>1) {
          message(paste("Analyzing Piece",pcs,"snp",((pcs-1)*k)+j,":",mysnp,"of",genoend-genostart))
          mf <- as.formula(paste("censoredOutcome ~ ",mysnp,covar.string))
          assign("groupID",gwdat[,famid],inherit=TRUE,envir=.GlobalEnv)
          sCOX <- try(coxme(mf ,random= ~ 1|groupID,data=gwdat))
          if ("try-error" %in% class(sCOX)) {
            ## COX failed as well.
            message(paste(mysnp,"COXME failed.\n"))
            rst <- c(mysnp,NA,NA,NA,NA,NA,NClust,CASE,1)
            rst <- data.frame(t(rst))
            names(rst) <- names(p.cox)
            p.cox <- rbind(p.cox,rst)
          } else {
            rCOX <- summary.machcoxme(sCOX)
            if (options()$verbose) print(sCOX)
            rst <- c(mysnp,signif(rCOX$coef,3),signif(exp(rCOX$coef),3),signif(rCOX$se,4),
                     round(rCOX$coef/rCOX$se,2),signif(1-pchisq((rCOX$coef/rCOX$se)^2,1),2),NClust,CASE,0)
            rst <- data.frame(t(rst))
            names(rst) <- names(p.cox)
            p.cox <- rbind(p.cox,rst)
          }
        }
        ## when outcome only exist in one category for the marker
        else {
          message(paste("No variation of outcome for SNP",mysnp,"for",outcome))
          rst <- c(mysnp,NA,NA,NA,NA,NA,NClust,CASE,2)
          rst <- data.frame(t(rst))
          names(rst) <- names(p.cox)
          p.cox <- rbind(p.cox,rst)
        }
      } ## Loop through all SNPs in piece k_th
      message(outcome,":Chr",chr,":Pc",pcs,"of",allpcs,"SNP",mysnp,":",refpoint+j,"of",totalsnp,"analyzed!\n",sep=" ")
      write.table(p.cox,sep=",",col.names=F,file=outfilename,row.names=F,quote=F,append=T)
      rm(data,fid.col,allsnp,refpoint,mysnp,datfm,gwdat,af,NClust,snpcall,cellcount,CASE,sCOX,rCOX,p.cox,rst)
      p.cox <- data.frame(SNP=character(0),coef=double(0),expCoef=double(0),se=double(0),
                          z=double(0),p=double(0),NCluster=double(0),FREQ=character(0),Memo=double(0),stringsAsFactors=F)
      gc()
      begin <- last+1
      last <- min(begin+k-1,genoend)
      pcs <- pcs+1
      if (options()$verbose) cat("Repeat while loop\n")
    } # while loop speed up data access
    outfilename <- paste(outdir,paste(paste(outcome,suffix,sep="."),chr,"txt",sep="."),sep="/")
  }# for each outcome
}
