################################################################################
## PROGRAM: splitPedigree.r
## BY: Bhoom Suktitipat, MD,PhD
## DATE: Tue Dec 20 15:38:18 EST 2011
################################################################################
## GOAL: Functions for splitting the specified pedigree to a smaller pedigree
##       - newFam() recreate a new pedigree id from a set of ID's
##       - makeFounder() make the specified 'id' in 'pedid' a founder
##       - newParents() assign new parent ID c(PAT,MAT) to 'id' in 'pedid'
##       - dupID(): duplicate individual to make a complete new pedigree
##       - .x in all functions is the pedigree data in LINKAGE FORMAT
##       - splitwhite() : split EA data from GeneSTAR for PhD Dissertation
##       - splitblack() : split AA data from GeneSTAR for PhD Dissertation
## OUTPUT: new pedigree data frame with split family
################################################################################
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
splitwhite <- function(pedigree,newpedfile,newdatfile,verbose=options()$verbose) {
  ##HALF-SIBS
  ##if (verbose)  print("Switching Halfsibs family for ID 6260501")
  ##pedigree[pedigree[,2]==6260501,3]<-62699
  ##if (verbose)  print("Removed previous father 62697 of 6260501")
  ##pedigree<-pedigree[-grep(62697,pedigree[,2]),]
  ##pedigree[pedigree[,2]==25214001,3:4] <- c(25114001,25114012)
  ## All split pedigree starts with 98*
  ## Split PEDIGREE 482 as follow to 98482
  ##***************************[STR: LINKAGE ANALYSIS DATA] Comment OUT *************************************#
  ##pedigree<-newFam(pedigree,482,98482,c(4820301,4820312,4820321,4820322,4820323,4820324,4820325,4820326))
  ##make 4820301 a founder
##pedigree <- makeFounder(pedigree,98482,4820301)
  ## Split PEDIGREE 1452
  ##pedigree<-newFam(pedigree,1452,981452,c(14520725,14520712,14520701,82510101,82510201,82510301,82510401))
  ##pedigree<-makeFounder(pedigree,981452,14520701)                      
  ## Split PEDIGREE 11521232
  ##pedigree<- newFam(pedigree,11521232,981232,c(123299,123298,21011501,22412301,23312301,23312312,22412323,22412324,23312321,23312322))
  ##pedigree<-makeFounder(pedigree,981232,21011501)
### Split PEDIGREE 415
  ##pedigree <- newFam(pedigree,415,98415,c(4150021,4150221,4150222,4150223,4150001,4150012,4150101,4150201,4150212,41598,41599))
  ##pedigree <- dupID(pedigree,98415,415,41599)
  ##**************************** [STR: LINKAGE ANALYSIS DATA] Comment IN ************************************#
  pedigree<-newFam(pedigree,482,98482,c(4820301,4820312,4820321,4820322,4820323,4820324,4820325,4820326))
  pedigree <- makeFounder(pedigree,98482,4820301)
  pedigree<-newFam(pedigree,1452,981452,c(14520725,14520712,14520701))
  pedigree<-makeFounder(pedigree,981452,14520701)                      
  pedigree<- newFam(pedigree,11521232,981232,c(123299,123298,21011501,22412301,23312301,23312312,22412323,22412324,23312321,23312322))
  pedigree<-makeFounder(pedigree,981232,21011501)
  pedigree <- newFam(pedigree,415,98415,c(4150021,4150221,4150222,4150223,4150001,4150012,4150101,4150201,4150212,41598,41599))
  pedigree <- dupID(pedigree,98415,415,41599)
  return(pedigree)
  ##write.table(pedigree2,file=newpedfile,col.names=F,row.names=F,na="x",quote=F)
  ##write.table(dat,file=newdatfile,col.names=F,row.names=F,na="x",quote=F)
}

splitblack <- function(pedigree,newpedfile,newdatfile,verbose=FALSE) {
  ## TWINS --> removed phenotype of one of the MZ twins
  ## -- removed contaminated samples 70110201
  ## -- removed contaminated samples 70420101
  ## - TBD -
  if (verbose)  print("# HALF-SIBS")
  cat("Switching halfsibs 71120301 family\n")
  pedigree[pedigree[,2]==71120301,3:4]<-c(711295,711298)
  cat("Switching halfsibs 71120401 family\n")
  pedigree[pedigree[,2]==72130401,3:4]<-c(721399,721398)
  cat("Switching halfsibs 82159191 family\n")
  pedigree[pedigree[,2]==82150101,3:4]<-c(821599,821598)
  cat("Switching halfsibs 71830301 family\n")
  pedigree[pedigree[,2]==71830301,3:4]<-c(718399,718398)
  cat("Switching halfsibs 71850621 family\n")
  pedigree[pedigree[,2]==71850621,3:4]<-c(71850612,71850601)
  cat("Fix incorrect nuclear family 70680212\n")
  pedigree[pedigree[,2]==70680121,3:4]<-c(70680812,70680801)
  cat("SPLIT PEDIGREE 201\n")
  pedigree<-newFam(pedigree,201,98201,c(20199,2010001,2010212,2010201,2010214,2010401,2010413,2010221,2010222,2010241,2010242,2010432,2010433))
  ## Probably have to add one recode to duplicate 20198
  pedigree<-dupID(pedigree,201,98201,20198)
  ## Split PEDIGREE 255
  pedigree<-newFam(pedigree,255,98255,c(2550112,2550101,2550121,2550122,2550123,2550124))
  pedigree<-makeFounder(pedigree,98255,2550101)
  ## ----------------------- COMMENT OUT FOR LINKAGE -----------------##
  ## Split PEDIGREE 70827184
  ##pedigree <- newFam(pedigree,70827184,7184,c(718499,718498,71840101,71840112,71840201,70820522))
  ## Duplicate ID
  ##pedigree <- dupID(pedigree,7184,70827184,71840101)
  ## ----------------------- COMMENT OUT FOR LINKAGE -----------------##
  return(pedigree)
  ##write.table(ped4,file=newpedfile,col.names=F,row.names=F,na="x",quote=F)
  ##write.table(dat,file=newdatfile,col.names=F,row.names=F,na="x",quote=F)
}
