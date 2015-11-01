## -------------------------------------------------------------------- ##
## Use Gene-Dropping Algortihm to Simulate SNPs for any given pedigree
## Bhoom Suktitipat, MD, PhD
## 10/02/2012
## -------------------------------------------------------------------- ##
## INPUT
## ped: data frame with the first 5 columns 1) family ID, 2) unique personal ID, 3) Father ID, 4) Mother ID, 5) Sex
##      Founders have to have both parents coded as 0
## MAF: MAF of the SNPs to be simulated, for multiple SNPs use c( 0.1, 0.2, 0.3, 0.4, 0.5)
## Prefix: Prefix of the markername to be simulated, e.g. "SNP"
## recodeA: Recode the simulated SNPs according to the number of minor allele
simSNP <- function(ped=ped,  MAF=0.01,  markerName="SNP",  recodeA=FALSE,  verbose=FALSE) {
  if ( length(unique(ped[,2])) < nrow(ped) ) stop("Individual ID is not unique")
  allele.drop <- function(ped,patid,matid,prop.pat,prop.mat,allele1.name,allele2.name) {
    ## Internally called function to do gene-dropping algorithm
    ## OPTIONS: 
    ##  ped = pedigree data frame containing unique individual ID
    ##  patid = father ID
    ##  matid = mother ID
    ##  prop.pat = values to determine which paternal allele will be inherited
    ##  prop.mat = values to determine which maternal allele will be inherited
    ##  allele1.name = column name of the inherited paternal allele
    ##  allele2.name = column name of the inherited maternal allele
    ## OUTPUT: 
    ##  Return two alleles for the offspring of "patid" and "matid"
    if (prop.pat <0.5) {
      ## if prop.pat < .5 inherit paternal's left allele (allele1)
      a1 <- ped[which(ped[,2]==patid),allele1.name]
    } else {
      a1 <- ped[which(ped[,2]==patid),allele2.name]
    }
    if (prop.mat <0.5) {
      ## if prop.mat < .5 inherit paternal's left allele (allele1)
      a2 <- ped[which(ped[,2]==matid),allele1.name]
    } else {
      a2 <- ped[which(ped[,2]==matid),allele2.name]
    }
    alleles <- c(a1,a2)
    return(alleles) ## return two alleles for the offspring of "patid" and "matid"
  }
  
  ## Keep the first five columns of the pedigree 
  .ped <- ped[,1:5]
  if (length(markerName)>1) markerName<- markerName[1]
  markers <- list()
  for (j in 1:length(MAF)) {
    ## Randomly choose allele 1 and allele 2
    if (runif(1)<0.5) {
      if (runif(1)<0.5) {
        allele=c("A","G")
      } else
      {
        allele=c("G","A")
      }
    } else {
      if (runif(1)<0.5) {
        allele=c("C","T")
      } else {
        allele=c("T","C")
      }
    }
    maf <- MAF[j]
    markername <- paste(markerName[1],j,sep="")
    cat("Simulating",markername,"\n")
    allele1 <- paste(markername,"_1",sep="") ## column in .ped to store allele1 
    allele2 <- paste(markername,"_2",sep="") ## column in .ped to store allele2
    row.names(.ped) <- .ped[,2]
    founder <- .ped[which(.ped[,3]==0 & .ped[,4]==0),2] ## ID of founders
    ## assign 0 or 1 allele: if runif value is less than MAF assign 1 (minor allele)
    ## + 0 will convert the logical test to number
    .ped[which(.ped[,3]==0 & .ped[,4]==0),allele1] <- (runif(length(founder))<maf)+0 
    .ped[which(.ped[,3]==0 & .ped[,4]==0),allele2] <- (runif(length(founder))<maf)+0 
    others <- .ped[which(!.ped[,3]==0 & !.ped[,4]==0),2] ## ID of non-founders
    family <- unique(.ped[,1]) ## Family ID
    pat <- .ped[which(!.ped[,3]==0 & !.ped[,4]==0),3]
    mat <- .ped[which(!.ped[,3]==0 & !.ped[,4]==0),4]
    geta1 <- runif(length(others)) # choice for father allele (whether left or right allele from father will be inherited)
    geta2 <- runif(length(others)) # choice for mother allele
    for (i in 1:length(others)) {
      ## for any given pairs of father and mother, find a pair of alleles to be inherited to the offspring
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
    while (length(others)>0 & k <=maxit) {
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
    } ## end while loop
    if (length(others) >0 & k > maxit) stop("Genotype still missing;need to increase number of iteration")
    if (verbose) cat(i,"\n")
    if (recodeA) {
      markers[[markername]] <- .ped[,allele1]+.ped[,allele2]
      .ped[,allele1] <- NULL
      .ped[,allele2] <- NULL
    } else {
      ## assign allele values 0 = common allele (allele[1], 1=rare allele (allele[2]))
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
  if (recodeA) {
    newname.ped <- paste(markerName,1:length(MAF),sep="")
    oldname.ped <- names(ped)
    ped <- cbind(ped,do.call(cbind,markers))
    names(ped) <- c(oldname.ped,newname.ped)
  } else {
    newname.ped <- paste(rep(paste(markerName,1:length(MAF),sep=""),each=2),1:2,sep="_")
    oldname.ped <- names(ped)
    ped <- cbind(ped,do.call(cbind,markers))
    names(ped) <- c(oldname.ped,newname.ped)
  }
  return(ped)
}
recode <- function(ped,model="A") {
  ## From pedigree data frame with first 6 columns as 
  ##    FID, IID, PAT, MAT, SEX, AFF
  ##    follow by two columns per marker
  ## recodeA(ped) will recode the marker according to the speficied model based on minor allele
  
  ## check if the format seems right 
  ped.ncol = (ncol(ped))
  if ( ped.ncol %% 2 ) stop("check the number of column in pedfile")
  ped.structure <- ped[,1:6]
  genotype <- ped[,7:ped.ncol]
  nmarkers = (ped.ncol-6)/2
  if (length(model) > 1) message("More than one model is speficied. Assume Additive Model")
  recoded <- matrix(ncol=nmarkers,nrow=nrow(ped))
  for (i in 0:nmarkers-1) {
    alleles <- as.matrix(genotype[,(2*i)+1:2])
    allele <- names(sta <- sort(table(alleles)))
    alleles.f <- as.factor(alleles)
    alleles.n <- matrix(c(0,1)[alleles.f],ncol=2)
    recoded[,i+1] <- apply(alleles.n,1,sum)
  }
  model=toupper(model)
  if ( model == "A" ) {
    return(list(geno=recoded,ped=ped.structure))
  } else if ( model=="R") {
    recodeR <- ifelse(recoded>0,recoded-1,recoded)
    return(list(geno=recodeR,ped=ped.structure))
  } else if ( model=="D") {
    recodeD <- ifelse(recoded==1,recoded+1,recoded)
    return(list(geno=recodeD,ped=ped.structure))
  }
}

## 01/21/2013
geno012 <- function(genotype,sep.allele="/") {
  ## recode a matrix of character allele of size n x 2 to 
  ## a vector of number of minor alleles vector coded as 0,1, or 2
  ## required input as a matrix of genotypes with 2 columns/marker
  if ( ncol(genotype) %% 2 ) stop("Number of alleles is not a denomination of 2")
  genotype <- as.matrix(genotype)
  genotype <- sub(sep.allele,"",genotype)
  nsnp=ncol(genotype)/2
  recoded <- matrix(ncol=nsnp,nrow=nrow(genotype))
  for (i in 1:nsnp) {
    alleles <- as.matrix(genotype[,(2*(i-1)+1:2)])
    allele <- names(sta <- sort(table(alleles))) ## find minor allele by sorting by allele frequency
    alleles.f <- as.factor(alleles) ## convert to factor for conversion
    alleles.n <- matrix(c(0,1)[alleles.f],ncol=2) ## convert alllele to number of minor allele
    recoded[,i] <- apply(alleles.n,1,sum) ## sum number of minor allele in the genotype
  }
  return(recoded)
}