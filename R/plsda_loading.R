plsda_loading <- function (plsda){
  R <- NULL
  for(i in 1:ncol(plsda$T)){
    r <- plsda$P[,i]*(sqrt(t(plsda$T[,i])%*%plsda$T[,i]/(nrow(plsda$T)-1)))[1,1]
    R <- cbind(R,r)
  }
  plsda$R <- R
  return(plsda)
}
