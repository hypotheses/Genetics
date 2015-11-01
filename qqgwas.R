# OLD VERSION DOESN'T WORK WITH R2.9
# 10/07/09 : Edit to contain option to output as png; removed the printing of observed > expected in the graph
# 10/08/09 : Fix the limit of the plot from 7 to maximum value
# 10/13/09 : Additional Fix for limit
# 10/24/09 : Fixed the expected line location, using library("sma")
# 11/04/09 v.06: Remove embeded figure generating code
# v.07 used library("sma") to plot the qqline
# 12/23/09: v.08: Upgrade to R version 2.10 which removed sma library
# 02/14/10: v.11 -- add confidence band in the plots
# 03/08/11: v.12 -- add function for lambda, fixed extra dots from plot.qqline()

## Calculate labmda
lambda <- function(p) {
  p <- na.omit(p)
  z <- abs(qnorm(p/2))
  l_gif <- median((z)^2,na.rm=TRUE)/0.456
  message("Genomic Inflation Factor (lambda) = ",signif(l_gif,5))
  return(l_gif)
}

#************* Required: library "sma"
# install.packages("sma")
#-- plot.qqline from sma library
plot.qqline <- function (x, y, a = 0.25, ...)
{
    y <- quantile(y[!is.na(y)], c(a, 1 - a))
    x <- quantile(x[!is.na(x)], c(a, 1 - a))
#    points(x, y, ...)
    slope <- diff(y)/diff(x)
    int <- y[1] - slope * x[1]
    abline(int, slope, ...)
}
#-- end plot.qqline
qq.gwas <- function(x,graphTitle,hpoint,limit=7,...) { 
  observed <- sort(na.omit(x))
  lobs <- -(log10(observed))
  expected <- c(1:length(observed))
  lexp <- -(log10(expected/(length(expected))))
  lim.obs <- ceiling(max(lobs))
  lim.exp <- ceiling(max(lexp))
  rm(observed);gc(reset=TRUE)  
  ##  diff <- lobs-lexp
  explim<-ceiling(max(lim.exp,lim.obs))
  plot(c(0,explim), c(0,explim), col="red", lwd=3, type="n",
       xlab="Expected (-logP)", ylab="Observed (-logP)",
       xlim=c(0,lim.exp), ylim=c(0,lim.obs), las=1, xaxs="i", yaxs="i", bty="l",...)
       if (!missing(graphTitle)) title(main=graphTitle)
##   if (0) {
##     if (missing(hpoint)) hpoint<-min(lobs[diff>0])
##     abline(h=hpoint,lty="dashed",col="red")
##     text(4,hpoint+.2,paste("Observed > Expected:",length(lobs[lobs>hpoint])),pos=4,cex=0.8)
##     print(paste("Observed > Expected:",length(lobs[lobs>hpoint])))
##   }
  plot.qqline(lexp,lobs,col="red")
  points(lexp, lobs, pch=23, cex=.4, bg="black",...)
#---- REMOVE EMBEDED FIGURE GENERATING CODES
#  if ( pdfout || pngout) dev.off()
#---- REMOVE EMBEDED FIGURE GENERATING CODES	
}

## RUN ##
PREFIX=NULL
INFILE="PHE1.txt"
OUTFILE="QQPlot_PHE1.png"
P_ONLY=FALSE
batch::parseCommandArgs()
if ( P_ONLY ) {
  if ( ! is.null(PREFIX) ) {
    files <- dir()[grep(paste("^",PREFIX,sep=""),dir())]
    print(paste("Result files included in the plot",files))
    results <- list()
    for (i in 1:length(files)) {
      message("Loading chromosome ",i," in ",files[i]," p-values")
      results[[i]] <- read.table(files[i],header=TRUE,as.is=TRUE)
    }
    result <- do.call(rbind,results)
    rm(results);gc(reset=TRUE)
  } else {
    ## for P_ONLY no header required!
    result <- read.table(INFILE,header=TRUE,as.is=TRUE,col.names="P")
  }
} else {
  result <- read.table(INFILE,header=TRUE,as.is=TRUE)
  names(result) <- c("CHR","SNP","POSITION","P")
}
png(file=OUTFILE,width=6,height=6,units="in",res=150)
print(paste("Genomic control factor =",lambda(result$P)))
gc(reset=TRUE)
qq.gwas(result$P) 
dev.off()
