## -----------------------------------------------
## Transpose data: transposeR.r
## -----------------------------------------------
## Bhoom Suktitipat, MD,PhD
## Date: 2015-11-01
## USAGE Rscript transposeR.r *.txt
## Accept multiple filenames or wildcard arguments
## -----------------------------------------------
options <- commandArgs(trailingOnly = TRUE)
require(data.table)
files=system(paste("ls",paste(options,collapse = " ")),intern=TRUE)
for (file in files) {
  print(paste("Processing",file))
  data<-fread(file)
  tdata <- data[, data.table(t(.SD), keep.rownames=TRUE)]
  write.table(tdata,file=paste("trans_",file,sep=""),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,na = "")
}
