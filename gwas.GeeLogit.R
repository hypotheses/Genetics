## 02/08/11: adding threshold for singularity check threshold 
cor.snp <- function(y, x) {
  if (!is.numeric(y)) 
    y <- as.numeric(as.factor(y))
  return(sd(y) == 0 || abs(cor(y, x, use = "complete")) > 0.99999999)
}
gwasGeeLogit <- function(covar="",famid="FID",genodata=genotype,phendata=phenotype,outcomes,chr,
                         suffix="",outdir="./result/white",
                         k=1000,genostart=7,genoend=ncol(genodata),minprop=0.001,singularity.threshold=20,...) {
#combine phenotype with genotype
  ## Check if phenotype contain any characters
  pheno.type <- grep("character",sapply(phendata,class))
  if ( length(pheno.type) > 0) {
    message(paste(names(phendata)[pheno.type],"is not numeric"))
    stop
  }
  message("gwasGeeLogit version 02/08/11")
  system(paste("mkdir -p",outdir))
  require(Zelig)
  totalsnp <- genoend-genostart+1
  pedhead <- names(genodata)[1:5]
  cv <- gsub("\\s+","",covar,perl=T)
  cv <- gsub("^\\+","",cv,perl=T)
  cv <- gsub("\\+","-",cv,perl=T)
  for (outcome in outcomes) {
    outfilename <- paste(outdir,paste(paste(outcome,"-",suffix,sep=""),chr,"csv",sep="."),sep="/")
    cat("SNP,Estimate,Naive SE,Naive.z,Robust SE,Robust.z,Cluster,FREQ,Memo\n",sep="",file=outfilename)
    p.Gee <- data.frame(SNP=character(0),Estimate=double(0),NaiveS.E.=double(0),NaiveZ=double(0),
                        RobustS.E.=double(0),RobustZ=double(0),NCluster=double(0),CASE=character(0),Memo=double(0),stringsAsFactors=F)
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
        names(phendata)[fid.col] <- famid
      }
      data <- genodata[,c(pedhead,snp)]
      data[,names(phendata)[rem.col]] <- phendata[match(data$IID,phendata[,phenid.col]),rem.col]
      data <- data[order(data$FID),]
      allsnp=length(snp)
      singular <- list()
      for (j in 1:allsnp) {
        refpoint <- begin-genostart
        mysnp <- snp[j]
        mf <- as.formula(paste(outcome,"~",mysnp,covar))
        datfm <- as.formula(paste(outcome, "~",mysnp,"+",famid,covar))
        gwdat <- model.frame(datfm,data=data,na.action=na.omit)
        snpcount <- round(gwdat[,2])
        af <- table(round(gwdat[,2]))
        freq <-  paste(paste(names(af),af,sep=":"),collapse="-")
        cellcount <- table(round(gwdat[,2]),gwdat[,1])
        NClust <- length(unique(gwdat[,3]))
        ## Check to see if all independent variables are colinear -- here removing outcome and famid variables
        singular[[j]] <- determinant(as.matrix(t(gwdat[,c(-1,-3)]))%*%as.matrix(gwdat[,c(-1,-3)]))$modulus[[1]]
        ## Check if there is variation in the case status ## 
        if (dim(cellcount)[2]>1 & dim(cellcount)[1]>1) {
          prop <- min(cellcount)/sum(cellcount)
          casecount <- paste("Case",paste(dimnames(cellcount)[[1]],cellcount[,2],sep=":",collapse="|"),collapse="",sep="-")
          controlcount <- paste("Control",paste(dimnames(cellcount)[[1]],cellcount[,1],sep=":",collapse="|"),collapse="",sep="-")
          CASE <- paste(casecount,controlcount,sep="-")
          ## check singularity of the covariance matrix
          if (options()$verbose) message(paste("Determinant of covariance matrix = ",singular[[j]]))
          if (sum(cellcount[,1]!=0)>1 & sum(cellcount[,2]!=0)>1) {
            if (singular[[j]] < singularity.threshold | prop < minprop) {
              ## 
              sGLM <- try(glm(mf,family=binomial,data=data,na.action=na.omit,...))
              if ("try-error" %in% class(sGLM)) {
                if (options()$verbose == TRUE) message("Fail logistic regression. Result missing")
                rst <- c(mysnp,NA,NA,NA,NA,NA,NClust,CASE,1)
                ## Memo 1: Singularity Covariance Matrix : try(GLM) failed
                rst <- data.frame(t(rst))
                names(rst) <- names(p.Gee)
                p.Gee <- rbind(p.Gee,rst)
              } else {
                if (options()$verbose == TRUE) message("Return results from Logistic regression")
                rGLM <- summary(sGLM)
                coeff <- rGLM$coefficients[2,]
                coeff[4] <- NA
                coeff[5] <- NA
                coeff[1:3] <- format(coeff[1:3],digits=5)
                rst <- data.frame(t(c(mysnp,coeff,NClust,CASE,2))) # for case/control should be more useful to tabulate case-control
                ## Memo 2: Singularity Covariance Matrix: try(GLM) worked
                ## rst <- data.frame(t(rst))
                names(rst) <- names(p.Gee)
                p.Gee <- rbind(p.Gee,rst)
              }
            } else {
              ## non-singular covariance matrix : GEE should work
              if (options()$verbose == TRUE) message("Fit Logistic regression with GEE")
              sGee <- try(zelig(mf,model="logit.gee",data=gwdat,id=famid,corstr="exchangeable",...))
              if ("try-error" %in% class(sGee)) {
                if (options()$verbose == TRUE) message("GEE failed. Results = missing")
                sGLM <- try(glm(mf,family=binomial,data=data,na.action=na.omit,...))
                if ("try-error" %in% class(sGLM)) {
                  if (options()$verbose == TRUE) message("Fail logistic regression. Result missing")
                  rst <- c(mysnp,NA,NA,NA,NA,NA,NClust,CASE,3)
                  ## Memo 3: Non-singularity Covariance Matrix : try(GEE) failed, and try(GLM) failed
                  rst <- data.frame(t(rst))
                  names(rst) <- names(p.Gee)
                  p.Gee <- rbind(p.Gee,rst)
                } else {
                  if (options()$verbose == TRUE) message("Return results from Logistic regression")
                  rGLM <- summary(sGLM)
                  coeff <- rGLM$coefficients[2,]
                  coeff[4] <- NA
                  coeff[5] <- NA
                  coeff[1:3] <- format(coeff[1:3],digits=5)
                  rst <- data.frame(t(c(mysnp,coeff,NClust,CASE,4))) 
                  ## Memo 4: Non-singularity covariance matrix: try(GEE) failed, try(GLM) worked
                  names(rst) <- names(p.Gee)
                  p.Gee <- rbind(p.Gee,rst)
                }
                ## rst <- c(mysnp,"Estimate"=NA,"Naive S.E."=NA,"Naive z"=NA,"Robust S.E."=NA,
                ##          "Robust z"=NA,NClust,CASE,3)
                ## ## Non-singular, failed GEE
                ## rst <- data.frame(t(rst))
                ## names(rst) <- names(p.Gee)
                ## p.Gee <- rbind(p.Gee,rst)
              } else if ( sGee$error == 0) {
                if (options()$verbose == TRUE) message("Results based on GEE with logit link")
                rGee <- summary(sGee)
                coeff <- rGee$coefficients[2,]
                rst <- c(mysnp,format(coeff,digits=5),NClust,CASE,0) # for case/control should be more useful to tabulate case-control
                rst <- data.frame(t(rst))
                names(rst) <- names(p.Gee)
                p.Gee <- rbind(p.Gee,rst)
              } else {
                cat(paste(mysnp,"has GEE error\n",sGee,"\n------------------------------------------------------\n"))
                rst <- c(mysnp,NA,NA,NA,NA,NA,NClust,CASE,5)
                ## Memo 5: Non-singularity covariance matrix: try(GEE) worked, GEE has error?
                print(paste(mysnp," caused GEE error"))
                rst <- data.frame(t(rst))
                names(rst) <- names(p.Gee)
                p.Gee <- rbind(p.Gee,rst)
              }
            }
          } else {
            message("Low variation among cases/controls")
            casecount <- paste("Case",paste(dimnames(cellcount)[[1]],cellcount[,2],sep=":",collapse="|"),collapse="",sep="-")
            controlcount <- paste("Control",paste(dimnames(cellcount)[[1]],cellcount[,1],sep=":",collapse="|"),collapse="",sep="-")
            CASE <- paste(casecount,controlcount,sep="-")
            rst <- c(mysnp,NA,NA,NA,NA,NA,NClust,CASE,6)
            ## Memo 6: Low variation among cases/controls
            rst <- data.frame(t(rst))
            names(rst) <- names(p.Gee)
            p.Gee <- rbind(p.Gee,rst)
          }
        } ## IF CHECK IF BOTH CASE AND CONTROL EXISTS FOR A GIVEN SNP
        else {
          message("Low variation of case status")
          casecount <- paste("Case",paste(dimnames(cellcount)[[1]],cellcount[,2],sep=":",collapse="|"),collapse="",sep="-")
          controlcount <- paste("Control",paste(dimnames(cellcount)[[1]],cellcount[,1],sep=":",collapse="|"),collapse="",sep="-")
          CASE <- paste(casecount,controlcount,sep="-")
          rst <- c(mysnp,NA,NA,NA,NA,NA,NClust,CASE,7)
          ## Memo 7: No variation between cases/controls
          rst <- data.frame(t(rst))
          names(rst) <- names(p.Gee)
          p.Gee <- rbind(p.Gee,rst)
        }
        cat(outcome,":Chr",chr,":Pc",pcs,"of",allpcs,"SNP",mysnp,":",refpoint+j,"of",totalsnp,"analyzed!\n",sep=" ")
        summary(singular)
      } # for loop over k SNPs in piece pcs-th
      tmpfile <- paste(outdir,paste(paste(outcome,suffix,sep=""),chr,pcs,"RData",sep="."),sep="/")
      write.table(p.Gee,sep=",",col.names=F,file=outfilename,row.names=F,quote=F,append=T)
      rGee <- "0"
      coeff <- "0"                        
      rm(rGee,data,snp,mf,datfm,gwdat,coeff,NClust,freq,mysnp,rst,p.Gee)
      p.Gee <- data.frame(SNP=character(0),Estimate=double(0),NaiveS.E.=double(0),NaiveZ=double(0),
                          RobustS.E.=double(0),RobustZ=double(0),NCluster=double(0),CASE=character(0),stringsAsFactors=F)
      gc()
      begin <- last+1
      last <- min(begin+k-1,genoend)
      pcs <- pcs+1
      cat("Repeat while loop\n")
    } # while loop speed up data access
  }# for each outcome
  if (options()$verbose) return(singular)
  
}
