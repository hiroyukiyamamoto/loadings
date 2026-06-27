one_kpca <- function(X,K){
  m <- nrow(X)
  X <- scale(X)
  K <- (K + t(K)) / 2
  eig <- eigen(K, symmetric = TRUE)
  keep <- eig$values > sqrt(.Machine$double.eps) * max(abs(eig$values))
  if (!any(keep)) {
    stop("K has no positive eigenvalues.")
  }
  Z <- sweep(eig$vectors[, keep, drop = FALSE], 2, sqrt(eig$values[keep]), "*")
  K <- Z %*% t(Z)
  C <- t(X) %*% K %*% X/m
  C <- (C + t(C)) / 2
  eig <- eigen(C, symmetric = TRUE)
  index <- order(eig$values, decreasing = TRUE)
  ev1 <- eig$vectors[, index]
  T <- X %*% ev1
  C <- K %*% X %*% t(X) %*% K/m
  C <- (C + t(C)) / 2
  eig <- eigen(C, symmetric = TRUE)
  index <- order(eig$values, decreasing = TRUE)
  ev <- eig$vectors[, index]
  TK <- K %*% ev
  one_kpca <- NULL
  one_kpca$T <- T
  one_kpca$U <- TK
  one_kpca$P <- ev1
  return(one_kpca)
}
