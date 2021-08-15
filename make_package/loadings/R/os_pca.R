# ------------
#   OS-PCA
# ------------
os_pca <- function(X,Y,D,kappa=0.999, M=diag(1,nrow(X))){
  
  X <- scale(X)
  MX <- scale(M%*%X)
  
  E <- (1-kappa)*diag(1,ncol(MX))+kappa*t(MX)%*%t(D)%*%D%*%MX
  
  G <- chol(solve(E))
  W0 <- svd(G%*%t(MX)%*%MX)$v
  
  # score
  t <- X%*%W0
  
  R <- chol(E)
  z <- svd(t(MX)%*%MX%*%solve(R))$v
  W2 <- solve(R)%*%z
  
  S <- MX%*%W2
  
  list(P=W0, T=t, Q=W2, U=S)
  
}
