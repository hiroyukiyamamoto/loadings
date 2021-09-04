# ------------
#   OS-PCA
# ------------
os_pca <- function(X,D,kappa=0.999, M=diag(1,nrow(X))){

  MX <- scale(M%*%X)
  X <- scale(X)

  E <- (1-kappa)*diag(1,ncol(MX))+kappa*t(MX)%*%t(D)%*%D%*%MX

  G <- chol(solve(E))
  W0 <- svd(G%*%t(MX)%*%MX)$v

  # score
  t <- X%*%W0
  Mt <- MX%*%W0

  R <- chol(E)
  z <- svd(t(MX)%*%MX%*%solve(R))$v
  W2 <- solve(R)%*%z

  Ms <- MX%*%W2

  list(P=W0, T=t, MT=Mt, Q=W2, U=Ms)

}
