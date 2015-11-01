######## Yoonhee Kim (April. 2011)
#######################################################################
### P.plot function for plotting colorful manhattan plot
### input files   :  plink.qassoc (plink association result) file  
###                        or                  
###                  Not limited number of columns, but required below 3 columns
###
###                  CHR :: chromosome
###                  BP :: positions
###                  P :: P value of GWAS results
###                  hwe :: HWE pvalue (required for flag=T option)
###
#######################################################################
P.plot <- function(d, flag=F, hwe.cut=1e-6, suggestiveline=0, title=NULL,genomewideline=-log10(5e-8)) {  
  
  d$logp = -log10(d$P)
  ch.t<-data.frame(table(sort(d$CHR)))    ### Generating the table of number of markers in each chromosome
  d$BP<-d$BP/10^6 					              ### BP positions to Mb s
  
  s<-c();s[1]<-0
  for(i in 1:nrow(ch.t)){
    s[i+1]<-s[i]+ch.t$Freq[i]}       ## accumlated number of markers in previous chromosomes
  
  first.ch<-c();first.ch[1]<-1
  for(i in 1:(nrow(ch.t)-1)){
    first.ch[i+1]<-first.ch[i]+ch.t$Freq[i]}       ## accumlated number of markers in previous chromosomes (idicating the row)
  
  first.bp.ch<-d$BP[first.ch]
  
  last.bp.ch<-d$BP[s[-1]]
  
  last.each<-c();last.each[1]<-last.bp.ch[1]
  for(i in 1:(nrow(ch.t)-1)){
    last.each[i+1]<-last.each[i]+as.numeric(last.bp.ch[i+1])-as.numeric(first.bp.ch[i+1])}       ## accumlated position of markers in previous chromosomes
  
  
  region.length<-d$BP[s[-1]]- d$BP[first.ch]
  total.length<-sum(as.numeric(region.length))
  multiplier<-0.01
  
  spacer.length<-round((total.length * multiplier),digits=0)
  spacer.length.list<-spacer.length*(1:length(unique(sort(d$CHR))))-(1:length(unique(sort(d$CHR))))
  
  d$newBP<-rep(c(0,last.each[-length(last.each)]+spacer.length.list[-length(spacer.length.list)]),ch.t$Freq)
  d$pos<-d$newBP + (d$BP-rep(first.bp.ch,ch.t$Freq))
  
  
  #color<-(d$CHR%%2)+8
  color<-d$CHR
  if (flag) {     							      ## Do you want to flag by HWE violated markers?
    color[d$hwe < hwe.cut]<-"yellow"
  }
  
  ticks<-NULL
  
  for ( i in sort(unique(d$CHR))) {
    ticks<-c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])			}			
  ticklim<-c(min(d$pos),max(d$pos))				
  
  maxy=ceiling(max(d$logp,na.rm=TRUE))
  #print(maxy)
  #if (maxy<8){ maxy=8 }
  
  plot(d$pos,d$logp,ylab=expression(-log[10](italic(p))), xlab="Chromosome",ylim=c(0,maxy),col=c("black","grey"),pch=20,cex=0.8,xaxt="n",main=title)
  #plot(d$pos,d$logp,ylab=expression(-log[10](italic(p))), xlab="Chromosome",ylim=c(0,maxy),col=c(color),pch=20,cex=0.8,xaxt="n",main=title)	## dot plot	
  #plot(d$pos,d$logp,ylab=expression(-log[10](italic(p))), xlab="Chromosome",ylim=c(0,maxy),col=c(color),pch=20,xaxt="n",main=title,type="h")	## line plot	
  box()
  axis(1,at=c(ticks),sort(unique(d$CHR)))
  abline(h=genomewideline,lty=2,col="red")
  abline(h=-log10(0.05/nrow(d)),lty=2,col="black")   ### divided the number of markers used
  abline(h=suggestiveline,lty=2,col="blue")							
  
  
  
}
if (0) {
  femaledata <- read.table("Z:/Research/GeneSTAR/CHARGE_CRP_Inflammation/Result021811/Summary/Upload/GeneSTAR_cohort_female_031111.txt.gz",header=TRUE,as.is=TRUE)
  femaledata2 <- femaledata[,c("chr","SNPID","position","p","HWE_pval")]; names(femaledata2) <- c("CHR","SNP","BP","P","hwe")
  maledata <- read.table("Z:/Research/GeneSTAR/CHARGE_CRP_Inflammation/Result021811/Summary/Upload/GeneSTAR_cohort_male_031111.txt.gz",header=TRUE,as.is=TRUE)
  maledata2 <- maledata[,c("chr","SNPID","position","p","HWE_pval")]; names(maledata2) <- c("CHR","SNP","BP","P","hwe")
  png(filename="maleCRP.png",width=11,height=8.5,res=300,units="in")
  P.plot(d=maledata2,flag=TRUE,hwe.cut=1e-6,suggestiveline=5,title="Male CRP")
  dev.off()
  png(filename="femaleCRP.png",width=11,height=8.5,res=300,units="in")
  P.plot(d=femaledata2,flag=TRUE,hwe.cut=1e-6,suggestiveline=5,title="Female CRP")
  dev.off()
}

