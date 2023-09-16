pls_da <- function(X, Y, k){

  # autoscaling
  X <- scale(X)
  Y <- scale(Y,scale=FALSE)

  # sample size-1
  N <- nrow(X)-1

  Wx <- NULL; Wy <- NULL; P <- NULL; T <- NULL; U <- NULL
  for(i in 1:k){

    # singular value decomposition
    USVx <- svd(t(Y)%*%X/N)
    USVy <- svd(t(X)%*%Y/N)

    # weight vector
    wx <- USVx$v[,1]
    wy <- USVy$v[,1]

    Wx <- cbind(Wx,wx)
    Wy <- cbind(Wy,wy)

    # score
    T0 <- X%*%wx
    U0 <- Y%*%wy

    T <- cbind(T,T0)
    U <- cbind(U,U0)

    # loading
    p <- t(X)%*%T0[,1]/as.numeric(t(T0[,1])%*%T0[,1])
    P <- cbind(P,p)

    # deflation
    X <- X-T0[,1]%*%t(p)

  }

  plsda <- list()
  plsda$Wx <- Wx
  plsda$P <- P
  plsda$T <- T
  plsda$Wy <- Wy
  plsda$S <- U

  return(plsda)

}
