# Internal helper for symmetric generalized eigenvalue problems used by
# multipls_rog(). This is not a full generalized eigenvalue solver.
# It assumes B is positive definite and solves A v = lambda B v by
# Cholesky transformation to a standard symmetric eigenvalue problem.
.generalized_eigen_symmetric <- function(A, B) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  n <- nrow(A)

  if (n != ncol(A) || nrow(B) != n || ncol(B) != n) {
    stop("A and B must be square matrices with the same dimensions.")
  }

  R <- chol(B)
  R_inv <- backsolve(R, diag(1, n))
  C <- t(R_inv) %*% A %*% R_inv
  C <- (C + t(C)) / 2
  eig <- eigen(C, symmetric = TRUE)
  vectors <- R_inv %*% eig$vectors

  list(values = eig$values, vector = vectors, vectors = vectors)
}
