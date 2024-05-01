plsda_loading <- function (plsda){
  n <- nrow(plsda$T)
  PLSDA_loadings <- NULL
  p_PLSDA <- NULL
  for(i in 1:ncol(plsda$T)){
    loading <- plsda$P[,i]*(sqrt(t(plsda$T[,i])%*%plsda$T[,i]/(nrow(plsda$T)-1)))[1,1]
    p <- 2 * stats::pt(abs(loading) * sqrt(n - 2)/sqrt(1 - loading^2), n - 2, lower.tail = FALSE)
    PLSDA_loadings <- cbind(PLSDA_loadings, loading)
    p_PLSDA <- cbind(p_PLSDA, p)
  }
  plsda$loading$R <- PLSDA_loadings
  plsda$loading$p.value <- p_PLSDA
  return(plsda)
}
