## ----------------------------------------------------------------------------------
## SUMMARY OF MAP AND INTER MARKER DISTANCE OF THE SIMULATED GENOTYPES USED FOR TRAP
## TO TEST THE EFFECT OF DENSER VS LOOSER MARKER MAPS ON TYPE I ERROR RATE
## 03/16/2012
## BHOOM SUKTITIPAT, MD, PHD
## ----------------------------------------------------------------------------------
## LOAD TRAP MAP ASSIGNED TO TILE 
setwd("../Desktop")
hapmap300 <- read.table("hapmap300.CHR2.map",sep=" ",header=TRUE,as.is=TRUE) ## DENSE MAP
hapmap300t <- read.table("hapmap300t.CHR2.map",sep=" ",header=TRUE,as.is=TRUE) ## LOOSE MAP

## FUNCTION TO CALCULATE INTER-MARKER DISTANCE
winsorize <- function(x, q=0.05) {
 extrema <- quantile(x, c(q, 1-q))
 x[x<extrema[1]] <- extrema[1]
 x[x>extrema[2]] <- extrema[2]
 x
} 
intermarker.distance <- function(x,winsorized=TRUE) {
  x2 <- c(x[-1],x[length(x)])
  dist <- x2-x
  if (winsorized) dist <- winsorize(dist)
  return(dist)
}
## INTER-MARKER DISTANCE FOR CHROMSOME 1 or 2
dense.chr1 <- intermarker.distance(hapmap300$POSITION)
loose.chr1 <- intermarker.distance(hapmap300t$POSITION)
dense.chr1.bytile <- sapply(tapply(hapmap300$POSITION,hapmap300$TILENO,intermarker.distance),median)
loose.chr1.bytile <- sapply(tapply(hapmap300t$POSITION,hapmap300t$TILENO,intermarker.distance),median)
## SUMMARY STATISTICS
# 1) Number of Total SNPs
dim(hapmap300);dim(hapmap300t)
# 2) Mean Global Intermarker Distance
mean(dense.chr1); mean(loose.chr1)
# 3) Mean Intermarker distance per tile
mean(dense.chr1.bytile); mean(loose.chr1.bytile)
# 4) Number of Tiles
length(unique(hapmap300$TILENO)) ; length(unique(hapmap300t$TILENO))
# 5) Number of SNPs per tile
## Histogram distribution of number of SNPs per tile
hist(dense.pertile <- tapply(hapmap300$SNP,hapmap300$TILENO,length),xlab="SNPs per Tile",main="More Dense")
hist(loose.pertile <- tapply(hapmap300t$SNP,hapmap300t$TILENO,length),xlab="SNPs per Tile",main="Less Dense")
summary(dense.pertile); summary(loose.pertile)


