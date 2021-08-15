pls_svd <- function(X,Y){
  
  source("R/loadings/plsrog/pls_rog.R")
  pls_svd <- pls_rog(X,Y,kappa=0)
  
  return(pls_svd)
  
}