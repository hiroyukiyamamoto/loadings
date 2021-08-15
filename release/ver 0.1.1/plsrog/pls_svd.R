pls_svd <- function(X,Y){
  
  pls_svd <- pls_rog(X,Y,kappa=0)
  
  return(pls_svd)
  
}