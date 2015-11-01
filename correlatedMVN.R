## SIMULATING CORRELATED DATA 
## JD LONG's Blog
## http://www.cerebralmastication.com/2010/08/even-simpler-multivariate-correlated-simulations/

require(MASS)

#make this reproducible
set.seed(2)

#how many draws in our starting data set?
n <- 1e4

# how many draws do we want from this distribution?
drawCount <- 1e4

myData   <- rnorm(n, 5, .6)
yourData <- myData  + rnorm(n, 8, .25)
hisData  <- myData  + rnorm(n, 6, .4)
herData  <- hisData + rnorm(n, 8, .35)

ourData <- data.frame(myData, yourData, hisData, herData)

# now we have raw correlations in the 70s, 80s, and 90s. Funny how that
# works out
cor(ourData)

#build the mvrnorm and take draws
#this replaces all my normalizing and
#copula building in the previous example
myDraws <- mvrnorm(1e5, mu=mean(ourData),
                   Sigma=cov(ourData)
                   )
myDraws <- data.frame(myDraws)

#check the mean and sd
apply(myDraws, 2, mean)
apply(myDraws, 2, sd)

#let's look at the mean and standard dev of the starting data
apply(ourData, 2, mean)
apply(ourData, 2, sd)
# so myDraws contains the final draws
# let's check Kolmogorov-Smirnov between the starting data
# and the final draws

for (i in 1:ncol(ourData)){
  print(ks.test(myDraws[[i]], ourData[[i]]))
}

#look at the correlation matrices
cor(myDraws)
cor(ourData)

#it's fun to plot the variables and see if the PDFs line up
#It's a good sanity check. Using ggplot2 to plot

require(ggplot2)

# rearrange the data to be "tall" not "wide"
meltDraws <-melt(myDraws)
meltDraws$dataSource <- "simulated"
meltData <- melt(ourData)
meltData$dataSource <- "original"
plotData <- rbind(meltData, meltDraws)

qplot(value, colour=dataSource, data=plotData, geom="density")+ facet_wrap(~variable)

## Another example from other websites with 
## Linear Algebra explanation
## http://comisef.wikidot.com/tutorial:correlation
nA         <- 5      # number of assets
nT         <- 100    # number of obs
rho        <- 0.8    # correlation between two assets

## create uncorrelated observations
X         <- rnorm(nA * nT) * 0.05
dim(X)    <- c(nT, nA)

## check
pairs(X, col = grey(0.4)); cor(X)

## set correlation matrix
M         <- array(rho, dim = c(nA, nA))
diag(M)    <- 1
## compute cholesky factor
cF        <- chol(M)

## induce correlation, check
Y        <- X %*% cF
pairs(Y, col = grey(0.4)); cor(Y)
