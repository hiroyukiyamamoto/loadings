multipls_loading <- function(multipls){
  
  # each data set
  T <- multipls$T
  for(j in 1:length(T)){
    
    t <- multipls$T[[j]]
    s <- multipls$U
    
    Wx <- multipls$P[[j]]
    Wy <- multipls$Q
    
    # each component
    multiPLS_loading <- NULL
    p_multiPLS_loading <- NULL
    for(k in 1:nrow(t)){
      
      index <- 1:length(multipls$T)  
      lambda0 <- NULL
      
      for(i in index[-j]){
        u <- multipls$T[[i]]
        lambda0 <- (tau[j,i]*cov(t[,k], u[,k]))
        lambda0 <- (lambda0 + tau[j,ncol(tau)]*cov(s[,k], t[,k]))/(2*as.numeric(t(Wx[,k])%*%Wx[,k]))
      }
      
      r0 <- NULL
      for(i in index[-j]){
        u <- multipls$T[[i]]
        r0 <- tau[j,i]*sqrt(var(u[,k]))
      }
      
      loading <- (Wx[,k]*2*lambda0)/(r0+tau[j,ncol(tau)]*sqrt(var(s[,k])))
      
      n <- nrow(t)
      z_fisher <- (1/2)*log((1+abs(loading))/(1-abs(loading)))
      p <- 2*pnorm(z_fisher, 0, 1/sqrt(n-3), lower.tail=FALSE) 
      
      multiPLS_loading <- cbind(multiPLS_loading, loading)
      p_multiPLS_loading <- cbind(p_multiPLS_loading, p)
      
    }
    
    multipls$loading$R[[j]] <- multiPLS_loading
    multipls$loading$p.value[[j]] <- p_multiPLS_loading
    
  }
  
  return(multipls)
  
}