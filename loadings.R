template_function <- function(method_name){
  
  source("R/loadings/pca_loading.R")
  source("R/loadings/pls_loading.R")
  
  template_function <- NULL
  if (method_name=="pca"){
    template_function <- pca_loading
  }
  if (method_name=="pls"){
    template_function <- pls_loading
  }
  if (method_name=="plsrog"){
    template_function <- pls_loading
  }
  if (method_name=="ospca"){
    template_function <- pls_loading
  }

  # error   
  if (is.null(template_function)){
    message("Method name does not match.")
  }
  
  return(template_function)

}