unsvmultipls_loading <- function (unsvmpls) 
{
  unsvmultipls <- unsvmpls
  unsvmultipls$U <- unsvmpls$T[[length(unsvmpls$T)]]
  unsvmultipls$Q <- unsvmpls$P[[length(unsvmpls$P)]]
  unsvmultipls$T <- unsvmpls$T[-length(unsvmpls$T)]
  unsvmultipls$P <- unsvmpls$P[-length(unsvmpls$P)]
  unsvmultipls <- multipls_loading(unsvmultipls)
  unsvmultipls$T <- unsvmpls$T
  unsvmultipls$P <- unsvmpls$P
  unsvmultipls$U <- NULL
  unsvmultipls$Q <- NULL
  return(unsvmultipls)
}
