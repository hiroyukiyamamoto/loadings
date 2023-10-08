pls_svd <- function(X,Y){
@p <- ncol(Y)
  D <- diag(1, p)
  pls_svd <- pls_rog(X, Y, D, kappa = 0)
  return(pls_svd)
}