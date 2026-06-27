pls_loading <- function(pls){
  T <- data.frame(pls$T)
  P <- data.frame(pls$P)
  U <- data.frame(pls$U)
  Q <- data.frame(pls$Q)
  n <- nrow(T)
  PLS_loadings <- NULL
  p_PLS <- NULL
  for (i in 1:ncol(T)) {
    Wx <- P[, i]
    Wy <- Q[, i]
    S <- U[, i]
    lambda <- stats::cov(T[, i], U[, i])
    loading <- (sqrt(n - 1) * lambda * Wx)/as.numeric(sqrt(t(S) %*% S))
    p <- 2 * stats::pt(abs(loading) * sqrt(n - 2)/sqrt(1 - loading^2), n - 2, lower.tail = FALSE)
    PLS_loadings <- cbind(PLS_loadings, loading)
    p_PLS <- cbind(p_PLS, p)
  }
  pls$loading$R <- PLS_loadings
  pls$loading$p.value <- p_PLS
  return(pls)
}
