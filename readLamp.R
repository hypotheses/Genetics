## Read-in LAMP data

readLamp <- function(lampData,lampMap) {
  map <- read.table(lampMap,header=F,col.names=c("SNP","Allele"),as.is=TRUE,sep="\t")
  data <- read.table(lampData,sep="\t",header=c("ID",map[,1]),as.is=TRUE)
  data[,1] <- gsub(":","",data[,1]) ## Clean extra colon in ID
  oddRow<- seq(1,nrow(data),by=2)
  evenRow <- seq(2,nrow(data),by=2)
  dataOut=list(ANC1=data[oddRow,],ANC2=data[evenRow,])
  return(dataOut)
}
                     
