template_function <- function(method_name){
  
  template_function <- NULL
  if (method_name=="pca"){
    method_loadings <- pca_loading
  }
  if (method_name=="pls"){
    method_loadings <- pls_loading
  }
  if (method_name=="plsrog"){
    method_loadings <- pls_loading
  }
  if (method_name=="ospca"){
    method_loadings <- pls_loading
  }

  # error   
  if (is.null(template_function)){
    message("Method name does not match.")
  }
  
  return(template_function)

}