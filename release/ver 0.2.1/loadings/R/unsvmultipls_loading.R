unsvmultipls_loading <- function(unsvmultipls){
  
  unsvmultipls$U <- unsvmultipls$T[[length(unsvmultipls$T)]]
  unsvmultipls$Q <- unsvmultipls$P[[length(unsvmultipls$P)]]
  
  unsvmultipls <- multipls_loading(unsvmultipls)
  
  unsvmultipls$U <- NULL
  unsvmultipls$Q <- NULL
  
  return(unsvmultipls)
}