library(geigen)

### Multiset PLS-ROG
multipls_rog <- function (X,Y,tau,D,kappa = 0.999)
{
  P <- NULL
  p <- colSums(Y)
  for(i in 1:ncol(Y)){
    P <- cbind(P,Y[,i]/p[i])
  }
  P <- t(P)
  
  XX <- NULL
  for(i in 1:length(X)){
    xx <- as.matrix(X[[i]])
    XX[i] <- list(scale(xx,scale=TRUE))
  }
  
  Y <- scale(Y,scale=FALSE)
  
  Z <- c(X=XX,Y=list(Y))
  
  N <- nrow(Z[[1]])-1;ZZ <- NULL
  for(i in 1:length(Z)){
    zz <- NULL
    for(j in 1:length(Z)){
      z <- (tau[i,j]/N)*t(Z[[i]])%*%Z[[j]]
      zz <- cbind(zz,z)
    }
    ZZ <- rbind(ZZ,zz)
  }
  
  pnum <- 0
  for(i in 1:length(Z)){
    pnum <- pnum+ncol(Z[[i]])
  }
  
  #B <- diag(2,pnum)
  #B[(nrow(B)-ncol(Y)+1):nrow(B),(nrow(B)-ncol(Y)+1):nrow(B)] <- 2*((1-kappa)*diag(1,ncol(Y))+kappa*t(Y)%*%t(P)%*%t(D)%*%D%*%P%*%Y)
  
  B <- diag(1,pnum)
  if (kappa!=0){
    B[(nrow(B)-ncol(Y)+1):nrow(B),(nrow(B)-ncol(Y)+1):nrow(B)] <- ((1-kappa)*diag(1,ncol(Y))+kappa*t(Y)%*%t(P)%*%t(D)%*%D%*%P%*%Y) # Multiset PLS-ROG
  }else{
    B[(nrow(B)-ncol(Y)+1):nrow(B),(nrow(B)-ncol(Y)+1):nrow(B)] <- (1-kappa)*diag(1,ncol(Y)) # Multiset PLS
  }

  eig_plsrog <- geigen(ZZ,B,symmetric=TRUE)
  w_multiplsrog <- eig_plsrog$vector
  
  lambda <- eig_plsrog$values
  lambda_index <- order(lambda,decreasing=TRUE)
  
  W <- NULL
  index1 <- 1
  for(i in 1:length(Z)){
    W[[i]] <- w_multiplsrog[index1:(index1+ncol(Z[[i]])-1),lambda_index]
    index1 <- index1+ncol(Z[[i]])
  }
  Wx <- W[-length(W)]
  Wy <- W[[length(W)]]
  
  T <- NULL
  for(i in 1:length(Z)){
    T[[i]] <- Z[[i]]%*%W[[i]]
  }
  
  S <- T[[length(T)]]
  T <- T[-length(T)]
  
  #list(P = Wx, T = T, Q = Wy, U = S)
  
  multipls_rog <- NULL
  
  multipls_rog$P <- Wx
  multipls_rog$T <- T
  multipls_rog$Q <- Wy
  multipls_rog$U <- S
  
  return(multipls_rog)
  
}