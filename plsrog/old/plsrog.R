plsrog <- function(X,class, kappa){
    
  # response variable
  Y0 <- factor(class)
  Y <- model.matrix(~ Y0 + 0)
  
  # penalized matrix
  P <- NULL
  p <- colSums(Y)
  for(i in 1:ncol(Y)){
    P <- cbind(P,Y[,i]/p[i])
  }
  P <- t(P)
  
  # differential matrix
  g <- ncol(Y)
  D <- diff(diag(1,g))
    
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
    
  # PLS loading
  R <- NULL
  for(i in 1:(ncol(Y)-1)){
      lambdax <- cov(T[,i],S[,i])
      r <- (sqrt(N)*lambdax*Wx[,i])/as.numeric(sqrt(t(Wy[,i])%*%t(Y)%*%Y%*%Wy[,i]))
      R <- cbind(R,r)
  }
  
  # statistical test of PLS loading
  P <- NULL
  for(i in 1:(ncol(Y)-1)){
    p <- 2*pt(abs(R[,i])*sqrt(nrow(X)-2)/sqrt(1-R[,i]^2), nrow(X)-2, lower.tail=FALSE)
    P <- cbind(P,p)
  }
  
  all <- list(T,S,Wx,Wy,R,P)
}
