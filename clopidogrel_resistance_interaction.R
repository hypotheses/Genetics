## clopidogrel_resistance_interaction.R
## 05/02/2013

## Dataset containing ADP agonist activated MEA measure
## Clopidogrel resistant is defined as MEA > 50

library(Rcmdr)
library(aod)

## 1) Setting Up Analysis Environment
setwd("C:/Users/suktitipatb/Dropbox/Research/clopidogrel/") ## set working directory
qtl <-   read.table("clopidogrel.csv", 
                    header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE) ## Read in data set
win.graph() ## open a plot window

## 2) Coding for dominant model, recessive, and create a factor variable
## Trancate extreme ADP value
qtl$ADP[qtl$ADP>150] <- 119
qtl$NR <- 1*(qtl$ADP50U==0)

## coding for dominant model
kqtl$CYP2d <- ifelse(qtl$CYP2==2,1,qtl$CYP2)
qtl$CYP3d <- ifelse(qtl$CYP3==2,1,qtl$CYP3) ## same as CYP3
qtl$PONd <- ifelse(qtl$PON==2,1,qtl$PON)

## coding for recessive model
qtl$CYP2r <- ifelse(qtl$CYP2==1,0,ifelse(qtl$CYP2==2,1,qtl$CYP2))
qtl$CYP3r <- ifelse(qtl$CYP3==1,0,qtl$CYP3) 
qtl$PONr <- ifelse(qtl$PON==1,0,ifelse(qtl$PON==2,1,qtl$PON))

## coding as factor for plotting
qtl$PON.f <- factor(qtl$PON,label=c("QQ)","QR","RR"))


## 3) Interaction Test using Cordell's Likelihood Ratio Test, colEpistatic() function in scrime library
library(scrime)
y = qtl$NR

## 3.1) CYP2 (recessive model) vs PON
set.seed(10592)
(cyp2rxpon <- colEpistatic(as.matrix(qtl[,c("CYP2r","PON")]),qtl$NR,
                          genes=c("*2","PON"))) ## p-value=0.352
perm_cyp2rxpon <- (rep(NA,10000))
for (i in 1:10000) {
  cat("permuting",i,"\r")
  x <- runif(n=211)
  perm_cyp2rxpon[[i]] <- colEpistatic(as.matrix(qtl[,c("CYP2r","PON")]),y[order(x)],
                                      genes=c("*2","PON"))[[4]]
}
sum(perm_cyp2rxpon < cyp2rxpon[[4]])/10000 ## 0.5412

## 3.2) CYP3 vs PON
perm_cyp3xpon <- (rep(NA,10000))
set.seed(54862)
(cyp3xpon <- colEpistatic(as.matrix(qtl[,c("CYP3","PON")]),qtl$NR,
                         genes=c("*3","PON"))) ## 0.930
## 3.2.1) permutation test
for (i in 1:10000) {
  cat("permuting",i,"\r")
  x <- runif(n=211)
  perm_cyp3xpon[[i]] <- colEpistatic(as.matrix(qtl[,c("CYP3","PON")]),y[order(x)],
                                     genes=c("*3","PON"))[[4]]
}
sum(perm_cyp3xpon < cyp3xpon[[4]])/10000 ## 0.4013

## 3.3) CYP2 vs CYP3
perm_cyp23 <- (rep(NA,10000))
set.seed(510932)
## Unadjusted
(cyp2xcyp3 <- colEpistatic(as.matrix(qtl[,c("CYP2","CYP3")]),qtl$NR,
                           genes=c("*2","*3"))) ## 0.865
## Adjusted for age and sex
(cyp2xcyp3 <- cordellsTest(mat.snp=as.matrix(qtl[,c("CYP2","CYP3")]),
                           cl=qtl$NR,cov=as.matrix(qtl[,c("AGE","SEX")]),
                           genes=c("*2","*3"))) ## 0.865
## 3.3.1) permutation test
for (i in 1:10000) {
  cat("permuting",i,"\r")
  x <- runif(n=211)
  perm_cyp23[[i]] <- colEpistatic(as.matrix(qtl[,c("CYP2","CYP3")]),y[order(x)],
                                  genes=c("*2","*3"))[[4]]
}
sum(perm_cyp23 < cyp2xcyp3[[4]])/10000 ## 0.303

