pls_rog <- function(X, Y, D, kappa=0.999){
  P <- NULL
  p <- colSums(Y)
  for (i in 1:ncol(Y)) {
    P <- cbind(P, Y[, i]/p[i])
  }
  P <- t(P)
  X <- scale(X)
  Y <- scale(Y, scale = FALSE)
  N <- nrow(X) - 1
  C <- kappa * t(Y) %*% t(P) %*% t(D) %*% D %*% P %*% Y + (1 - kappa) * diag(1, ncol(Y))
  Rx <- chol(solve(C))
  Ry <- chol(C)
  USVx <- svd(Rx %*% t(Y) %*% X/N)
  USVy <- svd(t(X) %*% Y %*% solve(Ry)/N)
  Wx <- USVx$v
  Wy <- solve(Ry) %*% USVy$v
  T <- X %*% Wx
  S <- Y %*% Wy
  for (i in 1:ncol(T)) {
    if (stats::cov(T[, i], S[, i]) < 0) {
      S[, i] <- -S[, i]
      Wy[, i] <- -Wy[, i]
    }
  }
  list(P = Wx, T = T, Q = Wy, U = S)
}
