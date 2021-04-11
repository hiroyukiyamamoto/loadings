pca_loading <- function(pca){
  # loadings
  cr_sd <- summary(pca)$importance[1,]
  PC_weight <- pca$rotation
  lambda_sqrt <- matrix(1,nrow(PC_weight))%*%cr_sd
  PC_loading <- lambda_sqrt*PC_weight
  # p-value
  n <- nrow(pca$x)
  p_loading <- function(n,PC_loading){
    pv_PC_loading <- 2*pt(abs(PC_loading)*sqrt(n-2)/sqrt(1-PC_loading^2),n-2,lower.tail=FALSE)
    #return(pv_PC_loading)
  }
  pv_PC_loading <- apply(PC_loading,2,p_loading,n=n)
  
  pca$loading$R <- PC_loading
  pca$loading$p.value <- pv_PC_loading
  
  return(pca)
}


