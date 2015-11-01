#!/bin/sh

##########################
#--- Extract + strand SNPs with HWEP > 1e-6
##########################
#bsub plink --noweb --bfile /research/2nd_imp_data/shared/washu/genotype/plink/whiteR21 --extract aasnps_1000g.txt --make-bed --out whiteImp
#bsub plink --noweb --bfile /research/2nd_imp_data/shared/washu/genotype/plink/blackR21 --extract aasnps_1000g.txt --make-bed --out blackImp

##########################
## Download liftover files
##########################
#if [[ ! -e liftOver ]]; then wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver; fi
#if [[ ! -e hg18ToHg19.over.chain.gz ]]; then wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz ;fi

##########################
## Convert bim file to bed format
##########################
#awk '{print "chr"$1"\t"$4"\t"$4+1"\t"$2}' blackImp.bim > blackImp.bim.bed
#awk '{print "chr"$1"\t"$4"\t"$4+1"\t"$2}' whiteImp.bim > whiteImp.bim.bed

##########################
## Run liftover to get hg19 map position
##########################
## gzip whiteImp.bim blackImp.bim 
## liftOver whiteImp.bim.bed hg18ToHg19.over.chain.gz whiteImp.bim
## liftOver need a library
##./liftOver: /lib64/tls/libc.so.6: version `GLIBC_2.4' not found (required by ./liftOver)
## Use web interface instead genome.ucsc.edu/cgi-bin/hgLiftOver using default version
## minimum ratio of bases that must remap = 0.95
## download output to hglft_genome_1aaf_4653d0.bed -- 818154 records
## failed SNPs reported in hglft_genome_1aaf_4653d0.err -- 37415 records

##########################
## RE-extract hg19 SNP
##########################
#awk '{print $4}' hglft_genome_1aaf_4653d0.bed > hg19.snp
#bsub plink --noweb --bfile /research/2nd_imp_data/shared/washu/genotype/plink/whiteR21 --extract hg19.snp --make-bed --out whiteImp
#bsub plink --noweb --bfile /research/2nd_imp_data/shared/washu/genotype/plink/blackR21 --extract hg19.snp --make-bed --out blackImp

##########################
## Update map position in PLINK map file
##########################
## white
cat <<EOF > updateWhiteImp.R
## read old map
data <- read.table("whiteImp.bim",header=FALSE,sep="\t",as.is=TRUE)
## read new position
new_pos <- read.table("hglft_genome_1aaf_4653d0.bed",header=FALSE,as.is=TRUE)
data[,"newpos"] <- new_pos[match(data[,2],new_pos[,4]),2] ## add new_pos to the map file (col 2= start position of each SNP)
## check if all SNPs has a new position
if ( nrow(data) == nrow(na.omit(data)) ) {
 ## output, CHR, SNP (newname as chr_pos), 0, newposition, a1, a2
 data[,2] <- paste(data[,1],data[,"newpos"],sep="_") 
 write.table(data[,c(1,2,3,7,5,6)],file="whiteImp.newbim",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
} else {
 message("Please check for unupdated SNPs")
} 
EOF
## R --no-save --no-restore < updateWhiteImp.R | tee updateWhiteImp.Rout
## black
cat <<EOF > updateBlackImp.R
## read old map
data <- read.table("blackImp.bim",header=FALSE,sep="\t",as.is=TRUE)
## read new position
new_pos <- read.table("hglft_genome_1aaf_4653d0.bed",header=FALSE,as.is=TRUE)
data[,"newpos"] <- new_pos[match(data[,2],new_pos[,4]),2] ## add new_pos to the map file
## check if all SNPs has a new position
if ( nrow(data) == nrow(na.omit(data)) ) {
 ## output, CHR, SNP (newname as chr_pos), 0, newposition, a1, a2
 data[,2] <- paste(data[,1],data[,"newpos"],sep="_") 
 write.table(data[,c(1,2,3,7,5,6)],file="blackImp.newbim",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
} else {
 message("Please check for unupdated SNPs")
} 
EOF
## R --no-save --no-restore < updateBlackImp.R | tee updateBlackImp.Rout

##########################
## create plink file by chromosome
## -- also remove individual with missing data more than 10% (founders and other individuals without genotype data)
##########################
echo -ne "conversion started: " >> conversion_date.txt
date >> conversion_date.txt
mkdir -p ea_chr
mkdir -p aa_chr
for i in {22..1}; do
  plink --noweb --bed whiteImp.bed --bim whiteImp.newbim --fam whiteImp.fam --chr $i --mind .1 --recode --out ea_chr/ea_b37_fw${i} 
  plink --noweb --bed blackImp.bed --bim blackImp.newbim --fam blackImp.fam --chr $i --mind .1 --recode --out aa_chr/aa_b37_fw${i} 
done
echo -ne "conversion finished: " >> conversion_date.txt
date >> conversion_date.txt
##########################
## Prephase with shapeit ?
##########################

##########################
## Convert to impute format using gtool ?
##########################
#wget http://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool_v0.7.5_x86_64.tgz
#tar xvzf gtool_v0.7.5_x86_64.tgz 
