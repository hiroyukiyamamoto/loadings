library(geigen)

source("C:/Users/yamamoto.HMT/Documents/R/multiomics/multipls_rog.R")

### Multiset PLS
multipls_geigen <- function(X,Y,tau){
  p <- ncol(Y)
  D <- diag(1, p)
  multipls <- multipls_rog(X, Y, tau, D, kappa = 0)
  return(multipls)
}
