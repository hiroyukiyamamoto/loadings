onekpca_loading <- function(onekpca){
  k <- min(ncol(onekpca$T), ncol(onekpca$U))
  n <- nrow(onekpca$T)
  ONEKPCA_loadings <- NULL
  p_ONEKPCA <- NULL
  for (i in 1:k) {
    lambdax <- (1/2) * stats::cov(onekpca$T[, i], onekpca$U[,i])
    loading <- (2 * lambdax * onekpca$P[, i]/sqrt(stats::var(onekpca$U[,i])))
    p <- 2 * stats::pt(abs(loading) * sqrt(n - 2)/sqrt(1 - loading^2), n - 2, lower.tail = FALSE)
    ONEKPCA_loadings <- cbind(ONEKPCA_loadings, loading)
    p_ONEKPCA <- cbind(p_ONEKPCA, p)
  }
  onekpca$loading$R <- ONEKPCA_loadings
  onekpca$loading$p.value <- p_ONEKPCA
  return(onekpca)
}
