rcca_da <- function(X, Y, tau, k){

  # autoscaling
  X <- scale(X)
  Y <- scale(Y,scale=FALSE)

  # sample size-1
  N <- nrow(X)-1

  Wx <- NULL; Wy <- NULL; P <- NULL; T <- NULL; U <- NULL
  for(i in 1:k){

    # Cholesky decomposition
    R <- chol((1-tau)*(1/N)*t(X)%*%X+tau*diag(1,ncol(X)))

    # singular value decomposition
    USVx <- svd(t(Y)%*%X%*%solve(R)*(1/N))
    #USVy <- svd(solve(R)%*%t(X)%*%Y*(1/N))

    # weight vector
    wx <- solve(R)%*%USVx$v[,1]
    #wy <- solve(R)%*%USVy$v[,1]

    Wx <- cbind(Wx,wx)
    #Wy <- cbind(Wy,wy)

    # score
    T0 <- X%*%wx
    #U0 <- Y%*%wy

    T <- cbind(T,T0)
    #U <- cbind(U,U0)

    # loading
    p <- t(X)%*%T0[,1]/as.numeric(t(T0[,1])%*%T0[,1])
    P <- cbind(P,p)

    # deflation
    X <- X-T0[,1]%*%t(p)

  }

  rcca_da <- list()
  rcca_da$Wx <- Wx
  rcca_da$P <- P
  rcca_da$T <- T
  #rcca_da$Wy <- Wy
  #rcca_da$S <- U

  return(rcca_da)

}
