ospca_loading <- function(ospca){

  MT <- data.frame(ospca$MT)
  P <- data.frame(ospca$P)
  U <- data.frame(ospca$U)
  Q <- data.frame(ospca$Q)

  n <- nrow(MT)

  OSPC_loadings <- NULL
  p_OSPCA <- NULL
  for(i in 1:ncol(MT)){
    Wx <- P[,i];Wy <- Q[,i]; S <- U[,i]
    lambda <- stats::cov(MT[,i], U[,i])
    #loading <- (sqrt(n-1)*lambda*Wx)/as.numeric(sqrt(t(Wy[,1])%*%t(Y)%*%Y%*%Wy[,1]))
    loading <- (sqrt(n-1)*lambda*Wx)/as.numeric(sqrt(t(S)%*%S))
    p <- 2*stats::pt(abs(loading)*sqrt(n-2)/sqrt(1-loading^2), n-2, lower.tail=FALSE)
    OSPC_loadings <- cbind(OSPC_loadings,loading)
    p_OSPCA <- cbind(p_OSPCA,p)
  }

  ospca$loading$R <- OSPC_loadings
  ospca$loading$p.value <- p_OSPCA

  return(ospca)

}


