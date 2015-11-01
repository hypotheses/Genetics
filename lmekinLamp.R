lmekin.Lamp <- function(genodata=genotype,phendata=phenotype,lampdata,k=1000,genostart=7,
                        genoend=ncol(genodata),outcomes,
                        famid="FID",covar="",suffix="",outdir="./lampResult",chr,kinmat,...) {
  ##combine phenotype with genotype
  system(paste("mkdir -p",outdir))
  require(kinship)
  totalsnp <- genoend-genostart+1
  pedhead <- names(genodata)[1:5]
  cv <- gsub("\\s+","",covar,perl=T)
  cv <- gsub("^\\+","",cv,perl=T)
  cv <- gsub("\\+","-",cv,perl=T)
  for (outcome in outcomes) {
    outfilename <- paste(outdir,paste(paste(outcome,suffix,sep=""),chr,"csv",sep="."),sep="/")
    ## FIX
    cat("SNP,EFFSNP,SESNP,TSNP,PSNP,EFFANC,SEANC,TANC,PANC\n",sep="",file=outfilename)
    p.LME <- data.frame(SNP=character(0),EFFSNP=double(0),SESNP=double(0),TSNP=double(0),PSNP=double(0),
                        EFFANC=double(0),SEANC=double(0),TANC=double(0),PANC=double(0),
                        stringsAsFactors=F)
    begin <- genostart
    last <- max(min(k+genostart-1,genoend),genostart)
    allpcs <- ceiling(totalsnp/k)
    pcs <- 1
    while (begin <= genoend) {
      snp <- names(genodata)[begin:last]
      names(phendata)[1:2] <- toupper(names(phendata)[1:2])
      if (! "ID" %in% names(phendata)) {
        cat("Variable name 'ID or id' not found in phenotype data. Please check\n")
        stop
      }
      data <- genodata[,c(pedhead,snp)]
      snpname <- gsub("_.$","",snp)
      col.snp <- match(snpname,names(lampdata))
      lamp <- lampdata[,c(1,col.snp)]
      names(lamp)[-1] <- paste(snpname,"a",sep="")
      data[,names(phendata)[-1]] <- phendata[match(data$IID,phendata[,1]),-1]
      data <- cbind(data,lamp[,-1])
      allsnp=length(snp)
      ##cat("before snp loop\n")
      for (j in 1:allsnp) {
        refpoint <- begin-genostart
        mysnp <- snp[j]
        anc <- paste(gsub("_.$","",mysnp),"a",sep="")
        mf <- as.formula(paste(outcome,"~",mysnp,"+",anc,covar))
        datfm <- as.formula(paste(outcome, "~",mysnp,"+IID+",anc,covar))
        gwdat <- model.frame(datfm,data=data)
#        af <- table(gwdat[,2])
#        NClust <- length(unique(gwdat[,3]))
#        freq <-  paste(paste(names(af),af,sep=":"),collapse="-")
#        cellcount <- table(gwdat[,mysnp])
#        allelecount <- paste(paste(cellcount,sep=":",collapse="|"),collapse="",sep="-")
#        gid <- gwdat[,"IID"]
#        assign("gid",gid,env=.GlobalEnv,inherits=TRUE)
        if (!length(unique(gwdat[,2]))==1) {
          ## lmekin(as_log_dgla ~ rs3802985_T + rs3802985a,random=~1|IID,varlist=kinmat,data=gwdat)
          ## lmekin(mf,random=~1|IID,varlist=kinmat,data=gwdat)
          SLK <- try(lmekin(mf,varlist=kinmat,random=~1|IID,data=gwdat))
          if ( !"try-error" %in% class(SLK) ) { 
            rSNP <- SLK$ctable[2,]
            rANC <- SLK$ctable[3,]
            rst <- c(mysnp,rSNP,rANC)
            names(rst) <- names(p.LME)
            p.LME <- rbind(p.LME,data.frame(t(rst)))
            if (options()$verbose) print(SLK)
          }
          else {
            cat(paste(mysnp,"has error\n------------------------------------------------------\n"))
            rst <- c(mysnp,NA,NA,NA,NA,NA,NA,NA,NA)
            names(rst) <- names(p.LME)
            names(rst) <- names(p.LME)
            p.LME <- rbind(p.LME,data.frame(t(rst)))
          }
        } ## IF try-error
        else {
          cat(paste(mysnp,"is a monomorphic SNP\n"))
          rst <- c(mysnp,NA,NA,NA,NA,NA,NA,NA,NA)
          names(rst) <- names(p.LME)
          p.LME <- rbind(p.LME,data.frame(t(rst)))
        }
        cat(outcome,":Chr",chr,":Pc",pcs,"of",allpcs,"SNP",mysnp,":",refpoint+j,"of",totalsnp,"analyzed!\n",sep=" ")
      } # for loop over k SNPs in piece pcs-th
      write.table(p.LME,sep=",",col.names=F,file=outfilename,row.names=F,quote=F,append=T)
      rGee <- "0"
      coeff <- "0"                        
      rm(SLK,data,snp,mysnp,mf,datfm,gwdat,rst,p.LME)
      p.LME <- data.frame(SNP=character(0),EFFSNP=double(0),SESNP=double(0),TSNP=double(0),PSNP=double(0),
                          EFFANC=double(0),SEANC=double(0),TANC=double(0),PANC=double(0),
                          stringsAsFactors=F)
      
      gc()
      begin <- last+1
      last <- min(begin+k-1,genoend)
      pcs <- pcs+1
      cat("Repeat while loop\n")
    } # while loop speed up data access
    ##  tmpfile <- paste(outdir,paste(paste(outcome,suffix,sep=""),chr,"RData",sep="."),sep="/")
    ##  outfilename <- paste(outdir,paste(paste(outcome,suffix,sep=""),chr,"txt",sep="."),sep="/")
    ##  system(paste("gzip -f",outfilename))
    ##  system(paste("rm -f",tmpfile))
  }# for each outcome
}
