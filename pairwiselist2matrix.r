## Pairwise-Table to Matrix
setwd("C:/Users/suktitipatb/Dropbox/ClinSeq/20120306datadump/Summary/")
library(reshape)
data <- read.table('C:/Users/suktitipatb/Dropbox/ClinSeq/20120306datadump/Summary/MIC.csv',sep=",",
                   header=TRUE,as.is=TRUE)
data <- data[order(data$r,decreasing=TRUE),];

## Find all variable
V1 <- unique(data$V1)
V2 <- unique(data$V2)
V <- unique(c(V1,V2))
N <- length(V)

pairMat <- matrix(nrow=N,ncol=N)
dimnames(pairMat)[[1]] <- V
dimnames(pairMat)[[2]] <- V


## Fill Matrix
for (i in 1:N) {
   pairwise <- data[which( data$V1 %in% V[i] | data$V2 %in% V[i]),]
   pairwise$V1 <- gsub(paste("^",V[i],"$",sep=""),"",pairwise$V1)
   pairwise$V2 <- gsub(paste("^",V[i],"$",sep=""),"",pairwise$V2)
   pairwise$V1V2 <- with(pairwise,paste(V1,V2,sep=""))
   pairMat[i,] <- pairwise[match(V,pairwise$V1V2),"value"]
}
(diag(pairMat) <- 1)

write.csv(pairMat,"C:/Users/suktitipatb/Dropbox/ClinSeq/20120306datadump/Summary/MIC_Matrix.csv",row.names=TRUE)
data<-read.csv("C:/Users/suktitipatb/Dropbox/ClinSeq/20120306datadump/Summary/MIC_Matrix.csv",as.is=TRUE)
row.names(data) <- data[,1]
dataMat<- as.matrix(data[,-1])
library(corrgram); install.packages("corrgram")
corrgram(dataMat,type="cor",order=TRUE,lower.panel=panel.shade,upper.panel=NULL,
         text.panel=panel.txt)
cmat <- dataMat
x.eigen <- eigen(cmat)$vectors[, 1:2]
e1 <- x.eigen[, 1]
e2 <- x.eigen[, 2]
alpha <- ifelse(e1 > 0, atan(e2/e1), atan(e2/e1) + pi)
ord <- order(alpha)
write.csv(cmat[ord,ord],"C:/Users/suktitipatb/Dropbox/ClinSeq/20120306datadump/Summary/MIC_Matrix.csv",row.names=TRUE)

## LOAD IMPUTED DATA (CHAINED EQUATION IN STATA)
## Drop these variable with > 5% missing rate
## drop echo_ivc echo_la2d diabetes echo_ascendingao echo_laarea 
##      echo_raarea cbc_mch cbc_mchc cbc_mpv thyroxine_ft4
data <- read.csv("../imputedQTL.csv",as.is=TRUE)
row.names(data) <- data[,1]
dataMat<- as.matrix(data[,c(-1,-82:-84)])
corrgram(dataMat,type="data",order=TRUE,lower.panel=panel.shade,upper.panel=NULL,
         text.panel=panel.txt)
cmat <- cmat <- cor(dataMat, use = "pairwise.complete.obs")
x.eigen <- eigen(cmat)$vectors[, 1:2]
e1 <- x.eigen[, 1]
e2 <- x.eigen[, 2]
alpha <- ifelse(e1 > 0, atan(e2/e1), atan(e2/e1) + pi)
ord <- order(alpha)
write.csv(cmat[ord,ord],"C:/Users/suktitipatb/Dropbox/ClinSeq/20120306datadump/Summary/ImputedData_corrMatrix.csv",row.names=TRUE)

