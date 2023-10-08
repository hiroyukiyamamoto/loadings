pca_loading <- function(pca){
  cr_sd <- summary(pca)$importance[1, ]
  PC_weight <- pca$rotation
  lambda_sqrt <- matrix(1, nrow(PC_weight)) %*% cr_sd
  PC_loading <- lambda_sqrt * PC_weight
  n <- nrow(pca$x)
  p_loading <- 2 * stats::pt(abs(PC_loading) * sqrt(n - 2)/sqrt(1 - PC_loading^2), n - 2, lower.tail = FALSE)
  pca$loading$R <- PC_loading
  pca$loading$p.value <- p_loading
  return(pca)
}
