#-----------------------------------------------------------------------------####
### This script is for checking phenotypes and generating the summary tables. ####
### SEP. 8.2010 ~ NOVEMBER 30, 2011
### Yoonhee Kim, BHOOM SUKTITIPAT
#-----------------------------------------------------------------------------####
# Summarize Continuous Variables
# REQUIRE INPUT: phenofile = phenotype file in CSV format 
#                with first column name "ID" follow by phenotype data
# OUTPUT: For each phenotype
# 1) PHE.pdf: Histograms of distribution + Quantile-Quantile Normal Plot
# 2) phenotypes_correlation.png: 
#    All pair-wise correlation correlation plot with loess smoothing plots
# 3) cor_table.csv: summary table of pair-wise correlation 
# 4) summary_untran_pheno.csv: summary table of summary statistics,
#                              N missing, kurtosis, skewness, N 
#-----------------------------------------------------------------------------####
summary.qtl <- function(phenofile) {
  require(moments)
  ##cor<-read.csv("cor_table.csv")
  pheno<-read.csv(file=phenofile,sep=",",header=TRUE,as.is=TRUE)
  ## remove ID
  ph<-pheno[,-1*grep("id",names(pheno))]
  ## Summary of Quantitative Phenotypes
  for ( i in c(1:ncol(ph))) {
  win.graph(width=8,height=11) ## open x11 graphics width=8, height=11
  par(mfrow=c(2,1),new=FALSE)
  h<-hist(ph[,i],plot=TRUE,label=TRUE,srt=90,xlab="",main=colnames(ph)[i])
  #h<-hist(ph[,i],breaks=50,main=colnames(ph)[i],plot=TRUE,label=TRUE,srt=90)
  #legend(max(ph[,i],na.rm=TRUE)*5/8, max(h$counts),names(summary(ph[,i])),box.col="white")
  #legend(max(ph[,i],na.rm=TRUE)*3/4, max(h$counts),summary(ph[,i]),box.col="white")
  qqnorm(ph[,i])
  qqline(ph[,i])
  #skew=agostino.test(ph[,i])
  #mtext(paste("skeness = ",round(skew[[1]][1],3),",p = ",
  #            round(skew[[2]],3),sep=""),line=-1)
  mtext(paste("skeness",round(skewness(ph[,i],na.rm=TRUE),3),sep=" = "),line=-1)
  mtext(paste("kurtosis",round(kurtosis(ph[,i],na.rm=TRUE),3),sep=" = "),line=-2)
  pname<-paste("plot_",names(ph)[i],sep="")
  savePlot(pname,type="pdf")
  dev.off()
  }
  
  ### Qualitative phenotype
#  par(mfrow=c(2,2))
#  for ( j in c(5,6,8,9)){
#  barplot(table(ph[,j]),legend.text=table(ph[,j]),main=colnames(ph)[j])
#  }
#  savePlot("quali.pdf",type="pdf")
  ### Corr matrix of phenotype variables
  
  panel.hist <- function(x, ...)
  {
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(usr[1:2], 0, 1.5) )
      h <- hist(x, plot = FALSE)
      breaks <- h$breaks; nB <- length(breaks)
      y <- h$counts; y <- y/max(y)
      rect(breaks[-nB], 0, breaks[-1], y,  ...)
  }
  
  panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
  {
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r <- abs(cor(x, y,"complete.obs"))
      p <- round(cor.test(x,y)$p.value,5)
      sig<-symnum(p,c(0,0.001,0.01,0.05,1),symbols=c("***","**","*"," "))
      txt <- format(c(r, 0.123456789), digits=digits)[1]
      txt <- paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex = 1 )
      text(0.7,0.7,sig, cex=2)
  
  }
  #win.graph(10,10)
  png(file="phenotypes_correlation.png",width=10,height=10,units="in",res=300)
  pairs(ph[,c(1:10)],diag.panel=panel.hist,upper.panel=panel.cor,lower.panel=panel.smooth)
  dev.off()
  #savePlot("phenotypes_correlation",type="pdf")
  pair.cor <- function(x, y)
  {
      r <- abs(cor(x, y,"complete.obs"))
      p <- round(cor.test(x,y)$p.value,5)
      sig<-symnum(p,c(0,0.001,0.01,0.05,1),symbols=c("***","**","*"," "))
    return(r,sig,p)
    cat(r,p,sig)  
  }

  #### Pairwise correlation matrix
  x <- ph
  n <- ncol(x) 
  cor.mat <- matrix(0, n, n) 
  for (i in seq(1, n)) { 
  cor.mat[i,i] <- sum( complete.cases(x[, i]) )   
  for (j in seq(i+1, length=n-i)) { 
  cor.mat[j,i]  <- abs(cor(x[, i], x[, j],"complete.obs" )) ### correlation
  cor.mat[i,j] <- round(cor.test(x[, i], x[, j])$p.value,5)  ## pvalue
  } 
  } 
  
  colnames(cor.mat)<-rownames(cor.mat)<-colnames(x)
  
  acc<-0
  table <- data.frame()
  
  for (i in seq(1, n)) { 
         for (j in seq(i+1, length=n-i)) { 
          row<-acc+j        
          table[row,1]<-colnames(x)[i]
          table[row,2]<-colnames(x)[j]
          table[row,3]  <- abs(cor(x[, i], x[, j],"complete.obs" )) ### correlation
          table[row,4] <- round(cor.test(x[, i], x[, j])$p.value,5)  ## pvalue
  } 
  acc<-acc+n-i
  } 
  
  table<-na.omit(table)
  T<-subset(table,table$V3>0.5 & table$V4<0.05)
  write.csv(table,"cor_table.csv")
  
  
  #### Kurtosis, Skewness test + summary table
  
  for ( i in 1:ncol(ph)) {
  kurtose<-anscombe.test(ph[,i])
  round(kurtose$statistic[1],3)
  round(kurtose$p.value,3)
  skew<-agostino.test(ph[,i])
  round(skew$statistic[1],3)
  round(skew$p.value,3)
  }
  
  
  totname<-vector(length=ncol(ph))
  tot<-matrix(nrow=ncol(ph),ncol=11)
  for ( i in 1:ncol(ph)) { 
     totname[i]<-colnames(ph)[i]
    stot<-summary(ph[,i])
     tot[i,1:6] <- stot
     tot[i,7] <- stot[7]
     kurtose<-anscombe.test(ph[,i])
     skew<-agostino.test(ph[,i])
     tot[i,8:11]<-c(round(kurtose$statistic[1],3),
                    round(kurtose$p.value,3),
                    round(skew$statistic[1],3),
                    round(skew$p.value,3))
  }
  dimnames(tot)[[1]]<-totname 
  dimnames(tot)[[2]]<-c("Min","1Q","Med","Mean","3Q","Max","Na's","Kurtosis","Kurtosis.P","Skewness","Skewness.P")
  write.csv(tot,"summary_untrans_pheno.csv")
}
#-----------------------------------------------------------------------------####
## EXAMPLE NOT RUN
#-----------------------------------------------------------------------------####
## > setwd("C:\Documents and Settings\suktitipatb.WISTERIA08\My Documents\Dropbox\Functional_Capacity\sas\Cognitive_vs_Symptom_Score")
## To generate the PDF file including histogram and qqplots with summary stats
## > phenofile="sz_Symptoms_vs_Cognition.csv"
## Phenofile is in CSV format with First Column == ID
## > summary.qtl(phenofile=phenofile)
