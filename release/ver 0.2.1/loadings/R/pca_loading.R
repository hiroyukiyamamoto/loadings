pca_loading <- function(pca){
  # loadings
  cr_sd <- summary(pca)$importance[1,]
  PC_weight <- pca$rotation
  lambda_sqrt <- matrix(1,nrow(PC_weight))%*%cr_sd
  PC_loading <- lambda_sqrt*PC_weight
  # p-value
  n <- nrow(pca$x)

  pv_PC_loading <- NULL
  for(i in 1:ncol(PC_loading)){
    p_loading <- 2*stats::pt(abs(PC_loading)*sqrt(n-2)/sqrt(1-PC_loading^2),n-2,lower.tail=FALSE)
    pv_PC_loading <- cbind(pv_PC_loading,p_loading)
  }
  pca$loading$R <- PC_loading
  pca$loading$p.value <- pv_PC_loading

  return(pca)
}


