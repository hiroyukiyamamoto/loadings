multipls_geigen <- function(X,Y,tau){
  p <- ncol(Y)
  D <- diag(1, p)
  multipls <- multipls_rog(X, Y, tau, D, kappa = 0)
  return(multipls)
}
