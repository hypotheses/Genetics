makeFamily <- function(nfam,nsibs){
  ## create 3 generations "nfam" extended families from 2 founds with "nsibs" siblings
  famsize <- 2+ (2* nsibs) + (nsibs^2)
  peds <- list() 
  founders <- c(1,2)
  siblings <- c(3:(2+nsibs))
  spouses  <- c((3+nsibs):(2+(2*nsibs)))
  (offspring <- c((3+2*nsibs):famsize))
  for (i in 1:nfam) {
    ## random sex
    iid <- c(1:famsize)
    sib.sex <- sample(1:2,size=nsibs,replace=TRUE,prob=c(.5,.5))
    spouse.sex <- ifelse(sib.sex==1,2,1)
    sex <- c(1,2,sib.sex,spouse.sex,sample(1:2,size=nsibs^2,replace=TRUE,prob=c(.5,.5)))
    fa <- c(0,0,rep(1,nsibs),rep(0,nsibs),rep(ifelse(sib.sex==1,siblings,spouses),each=nsibs))
    mo <- c(0,0,rep(2,nsibs),rep(0,nsibs),rep(ifelse(sib.sex==1,spouses,siblings),each=nsibs))
    peds[[i]] <- as.data.frame(cbind(rep(i,each=famsize),iid,fa,mo,sex),stringsAsFactors=FALSE)
  }
  ped <- do.call(rbind,peds)
  names(ped) <- c("FID","IID","PAT","MAT","SEX")
  ped$IID <- paste(ped$FID,"_",ped$IID,sep="")
  ped$PAT <- ifelse(ped$PAT==0,0,paste(ped$FID,ped$PAT,sep="_"))
  ped$MAT <- ifelse(ped$MAT==0,0,paste(ped$FID,ped$MAT,sep="_"))
  ped
}