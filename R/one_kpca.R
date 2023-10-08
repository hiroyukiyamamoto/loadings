one_kpca <- function(X,K){
  m <- nrow(X)
  X <- scale(X)
  eig <- eigen(t(X) %*% K %*% X/m)
  index <- order(eig$values, decreasing = TRUE)
  ev1 <- eig$vectors[, index]
  T <- X %*% ev1
  eig <- geigen::geigen(t(K) %*% X %*% t(X) %*% K/m, K, symmetric = FALSE)
  index <- order(eig$values, decreasing = TRUE)
  ev <- eig$vectors[, index]
  TK <- K %*% ev
  one_kpca <- NULL
  one_kpca$T <- T
  one_kpca$U <- TK
  one_kpca$P <- ev1
  return(one_kpca)
}
