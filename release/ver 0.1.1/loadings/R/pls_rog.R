pls_rog <- function(X, Y, D, kappa=0.999){
  
  # penalized matrix
  P <- NULL
  p <- colSums(Y)
  for(i in 1:ncol(Y)){
    P <- cbind(P,Y[,i]/p[i])
  }
  P <- t(P)
  
  # autoscaling
  X <- scale(X)
  Y <- scale(Y,scale=FALSE)
  
  # sample size-1
  N <- nrow(X)-1
  
  # smoothing parameter
  C <- kappa*t(Y)%*%t(P)%*%t(D)%*%D%*%P%*%Y+(1-kappa)*diag(1,g)
  
  # cholesky decomposition
  Rx <- chol(solve(C))
  Ry <- chol(C)
  
  # singular value decomposition
  USVx <- svd(Rx%*%t(Y)%*%X/N)
  USVy <- svd(t(X)%*%Y%*%solve(Ry)/N)
  
  # weght vector
  Wx <- USVx$v
  Wy <- solve(Ry)%*%USVy$v
  
  # score
  T <- X%*%Wx
  S <- Y%*%Wy
  
  list(P=Wx, T=T, Q=Wy, U=S)
}