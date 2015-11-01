# Making a Manhattan Plot of a data file
# Need columns with these names "CHR","BP","p","MAF","HWEP"
# Will plot to the standard output. Need to embed this function in pdf(), png(), etc.
# The plot does not have title, but can be added by running title("YOUR TITLE") after running manhattan.plot()

mht.plot <- function(P,CHR,POSITION,gwas.cutoff=5e-8,signal.from=5,plotcolor=c("blue","cyan"),signalcolor="red",limy,chr,...) {
# Recommended plot size: width=1280,height=.5*768
is.missing <- (is.na(CHR) | is.na(P) | is.na(POSITION))
CHR <- CHR[!is.missing]
P <- P[!is.missing]
POSITION <- POSITION[!is.missing]
print(paste(sum(is.missing), "SNPs with missing chromosomal position or P-values will not be plotted"))
if (gwas.cutoff < 1) {
 	GWAS.cutoff=-log(gwas.cutoff,10)
} 
else {
	GWAS.cutoff=gwas.cutoff
}
chrlist <- sort(unique(CHR))
chr.table <- data.frame("CHR"=chrlist,"chrMax"=as.numeric(tapply(POSITION,CHR,max)),
				"chrMin"=as.numeric(tapply(POSITION,CHR,min)))
chr.table$range <- chr.table$chrMax - chr.table$chrMin
chr.table$start <- 1
offscale <- 20000000
for ( i in 2:nrow(chr.table)) {
  chr.table[i,"start"] <- chr.table[i-1,"start"]+chr.table[i-1,"range"]+offscale
}
newscale <- POSITION+chr.table[match(CHR,chr.table$CHR),"start"]-chr.table[match(CHR,chr.table$CHR),"chrMin"]
chr.table$tick <- as.numeric(tapply(newscale,CHR,mean))
nlogp <- -log10(P)
##cat(dim(filtered.result)[1],"\n")
### MAKING PLOT FRAME
if (missing(limy)) limy=ceiling(max(nlogp))
plot(newscale,nlogp,type="n",xlab="Chromosome",ylab="-log10(p-value)",xaxt="n",xlim=c(0,max(newscale)),ylim=c(0,limy)) # ,log="y"
abline(h=GWAS.cutoff,lty="dashed",col="red")

### MAKING PLOT
mycolor <- plotcolor
odd<-chr.table$CHR[as.logical(chr.table$CHR %%2)]
even<-chr.table$CHR[!as.logical(chr.table$CHR %%2)]
cat("Plotting chromosome",odd,"\n")
points(newscale[CHR%in%odd],nlogp[CHR%in%odd],pch=1,cex=.5,lwd=.5,col=mycolor[1])
cat("Plotting chromosome",even,"\n")
points(newscale[CHR%in%even],nlogp[CHR%in%even],pch=1,cex=.5,lwd=.5,col=mycolor[2])
chrlab <- c(1:22,"X","XY","Y","MT")
chr <- length(unique(CHR))
axis(1,at=chr.table$tick[1:chr],label=chrlab[1:chr],las=3)
abline(h=signal.from,lty="dashed",col="darkblue")
myline<-signal.from
}

### Running Interactive Manhattan Plots
INFILE="PHE1.txt"
OUTFILE="MHTPlot_PHE1.pdf"
batch::parseCommandArgs()
if (INFILE == "") {
  message("Missing INFILE")
  stop('Example: R --vanilla --args INFILE "PHE1.txt" OUTFILE "MHTPlot_PHE1.pdf" < ~/script/R/mhtplot.v072011.R ')
}
result <- read.table(INFILE,header=TRUE,as.is=TRUE)
names(result) <- c("CHR","SNP","POSITION","P")
if ( OUTFILE == "") {
  message("Saving plot to MHTPlot.pdf")
  OUTFILE="MHTPlot.pdf"
}
png(file=OUTFILE,width=11,height=4,units="in",res=150)
mht.plot(result$P,result$CHR,result$POSITION)
dev.off()
