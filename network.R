## network example
## http://cran.r-project.org/web/packages/network/network.pdf
## From http://www.r-bloggers.com/how-to-plot-a-network-package-network-tip-2/

install.packages('network')
library(network)

## generate example data
generateBA = function(n = 100, n0 = 2){
  mat = matrix(0, nrow= n, ncol = n)
  for(i in 1:n0){
    for(j in 1:n0){
      if(i != j){
        mat[i,j] = 1
        mat[j,i] = 1
      }
    }
  }
  for(i in n0:n){
    list = c()
    for(k in 1:(i-1)){
      list = c(list, sum(mat[,k]))
    }
    link = sample(c(1:(i-1)), size = 1, prob = list)
    mat[link,i] = 1
    mat[i,link] = 1
  }
  return(mat)
}

artificialNet= generateBA(200)
a = network(artificialNet, directed = FALSE)
plot.network(a)

## network summary statistics
network.density(a)
network.size(a)
network.edgecount(a)

## create blank template
net = network.initialize(10)
plot(net)

## connect node 1 and 2
net[1,2] = 1
plot(net)