multipls_loading <- function(multipls) 
{
  T <- multipls$T
  for (j in 1:length(T)) {
    t <- multipls$T[[j]]
    s <- multipls$U
    Wx <- multipls$P[[j]]
    Wy <- multipls$Q
    tau <- multipls$tau
    multiPLS_score <- NULL
    multiPLS_loading <- NULL
    p_multiPLS_loading <- NULL
    for (k in 1:nrow(t)) {
      index <- 1:length(multipls$T)
      lambda0 <- 0
      for (i in index[-j]) {
        u <- multipls$T[[i]]
        lambda0 <- lambda0 + tau[j, i] * stats::cov(t[, k], u[, k])
      }
      lambda0 <- lambda0 + tau[j, ncol(tau)] * stats::cov(s[, k], t[, k])
      lambda0 <- lambda0 / (2 * as.numeric(t(Wx[, k]) %*% Wx[, k]))
      r0 <- 0
      for (i in index[-j]) {
        u <- multipls$T[[i]]
        r0 <- r0 + tau[j, i]*u[, k]
      }
      r0 <- r0 + tau[j, ncol(tau)]*s[, k]
      loading <- (Wx[, k] * 2 * lambda0)/sqrt(stats::var(r0))
      n <- nrow(t)
      p <- 2 * stats::pt(abs(loading) * sqrt(n - 2)/sqrt(1 - loading^2), n - 2, lower.tail = FALSE)
      multiPLS_score <- cbind(multiPLS_score, r0)
      multiPLS_loading <- cbind(multiPLS_loading, loading)
      p_multiPLS_loading <- cbind(p_multiPLS_loading, p)
    }
    multipls$loading$Score[[j]] <- multiPLS_score
    multipls$loading$R[[j]] <- multiPLS_loading
    multipls$loading$p.value[[j]] <- p_multiPLS_loading
  }
  multiPLS_score <- NULL
  multiPLS_loading <- NULL
  p_multiPLS_loading <- NULL
  for (k in 1:ncol(s)) {
    lambda0 <- 0
    for (i in 1:length(T)) {
      u <- multipls$T[[i]]
      lambda0 <- lambda0 + tau[length(T)+1, i] * stats::cov(u[, k], s[, k])
    }
    lambda0 <- lambda0 / (2 * as.numeric(t(Wy[, k]) %*% Wy[, k]))
    r0 <- 0
    for (i in 1:length(T)) {
      u <- multipls$T[[i]]
      r0 <- r0 + tau[length(T)+1, i]*u[, k]
    }
    loading <- (Wy[,k] * 2 * lambda0)/sqrt(stats::var(r0))
    n <- nrow(s)
    p <- 2 * stats::pt(abs(loading) * sqrt(n - 2)/sqrt(1 - loading^2), n - 2, lower.tail = FALSE)
    multiPLS_score <- cbind(multiPLS_score, r0)
    multiPLS_loading <- cbind(multiPLS_loading, loading)
    p_multiPLS_loading <- cbind(p_multiPLS_loading, p)
  }
  multipls$loading$Score[[length(T)+1]] <- multiPLS_score
  multipls$loading$R[[length(T)+1]] <- multiPLS_loading
  multipls$loading$p.value[[length(T)+1]] <- p_multiPLS_loading
  return(multipls)
}