gwasGEE <- function(genodata=genotype,phendata=phenotype,k=1000,genostart=7,
                       genoend=ncol(genodata),outcomes=names(phenotype)[-2:-1],
                       famid="FID",covar="",suffix="",outdir="./results",mincount=0 ,...) {
  ## -------------------------------------------------------------------------------------------------------------------- ##
  ## Analyze family data for association using GEE with exchangeable correlation structure and robust variance estimators
  ## By: Bhoom Suktitipat, MD, PhD
  ## Date: 07/25/2012
  ## -------------------------------------------------------------------------------------------------------------------- ##
  ## genodata : a data frame containing genotype data from PLINK obtained 
  ##            using plink --bfile [plink_file] --recodeA --out [output_file]
  ## phendata : a data frame containing phenotype data with "ID" column for combining with genodata
  ## k        : an integer used to opimize the size of data to be handled at a time. Default = 1000
  ## genostart: an integer denoting the first column of the genotype in genodata
  ## genoend  : an integer denoting the column number of the last marker to be analyzed, 
  ##            by default this is the last column of genodata
  ## outcomes : a character string or array containing the name of phenotypes in "phendata" to be analyzed
  ## famid    : a character string denoting the variable name in genodata that is corresponding to the cluster ID 
  ##            to be used with GEE model
  ## covar    : a character string containing covariates in the form
  ##            "+X1+X2+X3"
  ## suffix   : a character string to be added to the result filename in addition to the phenotype name
  ## outdir   : a character string denoting the output directory 
  ## mincount : an integer denoting the minimum genotype count of SNPs to be analyzed e.g. 
  ##            specifying mincount=5 will only analyze SNP with at least 5 counts for any genotype
  ## ...      : additional options to be passed onto gee() function
  ## -------------------------------------------------------------------------------------------------------------------- ##
  ## Example (NOT RUN)
  ## phenotype <- read.table(")
  ## DEBUGGING PARAMTERS 
  ## genodata=ped2;phendata=pheno;k=1000;genostart=7;genoend=ncol(ped2);outcomes="phen1";famid="FID";covar="";suffix="";outdir=".";mincount=0
  #combine phenotype with genotype
  system(paste("mkdir -p",outdir))
  require(gee)
  totalsnp <- genoend-genostart+1
  pedhead <- names(genodata)[1:5]
  cv <- gsub("\\s+","",covar,perl=T)
  cv <- gsub("^\\+","",cv,perl=T)
  cv <- gsub("\\+","-",cv,perl=T)
  for (outcome in outcomes) {
    ### VECTORS CONTAINING THE RESULTS
    NSNPs=genoend-genostart+1
    SNP=character(NSNPs);Est=numeric(NSNPs);NaiveSE=numeric(NSNPs);NaiveZ=numeric(NSNPs);
    RobustSE=numeric(NSNPs);RobustZ=numeric(NSNPs);CASE=character(NSNPs)
    begin <- genostart
    last <- max(min(k+genostart-1,genoend),genostart)
    allpcs <- ceiling(totalsnp/k)
    pcs <- 1
    while (begin <= genoend) {
      snp <- names(genodata)[begin:last]
      ## find ID column in phenotype file
      phen.id <- grep("^I?ID",names(phendata),ignore.case=TRUE)[1]
      if (length(phen.id)==0) {
        stop("Variable name 'ID or id' not found in phenotype data. Please check\n")
      }
      data <- genodata[,c(pedhead,snp)]
      data[,outcome] <- phendata[match(data$IID,phendata[,phen.id]),outcome]
      allsnp=length(snp)
      #            cat("before snp loop\n")
      for (j in 1:allsnp) {
        index=(pcs-1)*k+j ## index of the results
        N=allsnp ## total SNPs in the current piece of genotype data being analyzed
        refpoint <- begin-genostart
        mysnp <- snp[j]
        mf <- as.formula(paste(outcome,"~",mysnp,covar))
        datfm <- as.formula(paste(outcome, "~",mysnp,"+",famid,covar))
        gwdat <- model.frame(datfm,data=data,na.action=na.omit)
        af <- table(gwdat[,2])
        NClust <- length(unique(gwdat[,3]))
        ##freq <-  paste(paste(names(af),af,sep=":"),collapse="-")
        cellcount <- table(round(gwdat[,mysnp]))
        af <- table(round(gwdat[,2]))
        allelecount=paste(dimnames(af)[[1]],af,sep=":",collapse="|")
        if (!length(af)==1) {
          if (any(cellcount<mincount)) {
            cat(paste(mysnp,"has low genotype counts\n"))
            SNP[index]=mysnp;Est[index]=NA;NaiveSE[index]=NA;NaiveZ[index]=NA;
            RobustSE[index]=NA;RobustZ[index]=NA;CASE[index]=allelecount
          } else {
            # detect whether any cells is less than 5 #
            sGee <- try(gee(mf,family=gaussian,data=data,id=FID,corstr="exchangeable",na.action=na.omit,...))
            if ( ! "try-error" %in% class(sGee)) { 
              rGee <- summary(sGee)
              coeff <- rGee$coefficients[2,]
              SNP[index]=mysnp;Est[index]=coeff[1];NaiveSE[index]=coeff[2];NaiveZ[index]=coeff[3];
              RobustSE[index]=coeff[4];RobustZ[index]=coeff[5];CASE[index]=allelecount
            } else {
              cat(paste(mysnp,"has GEE error\n",sGee,"\n------------------------------------------------------\n"))
              SNP[index]=mysnp;Est[index]=NA;NaiveSE[index]=NA;NaiveZ[index]=NA;
              RobustSE[index]=NA;RobustZ[index]=NA;CASE[index]=allelecount
            }
          }
        } else {
          cat(paste(mysnp,"is a monomorphic SNP\n"))
          SNP[index]=mysnp;Est[index]=NA;NaiveSE[index]=NA;NaiveZ[index]=NA;
          RobustSE[index]=NA;RobustZ[index]=NA;CASE[index]=allelecount
        }
        cat(outcome,suffix," Pc",pcs,"of",allpcs,"SNP",mysnp,":",refpoint+j,"of",totalsnp,"analyzed!\n",sep=" ")
      } # for loop over k SNPs in piece pcs-th
      
      rm(rGee,data,snp,mf,datfm,gwdat,coeff,mysnp)
      gc()
      begin <- last+1
      last <- min(begin+k-1,genoend)
      pcs <- pcs+1
    } # while loop speed up data access
    ## WRITE OUT RESULTS 
    outfilename <- paste(outdir,paste(paste(outcome,suffix,sep=""),"csv",sep="."),sep="/")
    p.Gee <- data.frame("SNP"=SNP,"EST"=signif(Est,3),"NaiveSE"=signif(NaiveSE,3),"NaiveZ"=signif(NaiveZ,3),
                        "RobustSE"=signif(RobustSE,3),"RobustZ"=signif(RobustZ,3),
                        "RobustP"=signif(pnorm(abs(RobustZ),lower.tail=FALSE),3),
                        "N"=CASE,stringsAsFactors = FALSE)
    write.table(p.Gee,sep=",",col.names=TRUE,file=outfilename,row.names=FALSE,quote=FALSE)
    message("Finished analysis of ",outcome)
  }# for each outcome
}