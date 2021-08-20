rm(list=ls(all=TRUE))

library(loadings)

### PCA
data(fasting)

pca <- prcomp(X, scale=TRUE)

pca <- pca_loading(pca)
pca$loading$R
pca$loading$p.value

### PLS
data(whhl)

pls <- pls_svd(X,Y)
# library(chemometrics)
# pls <- pls_eigen(X,Y,a=3)

pls <- pls_loading(pls)
pls$loading$R
pls$loading$p.value

### PLS-ROG
data(whhl)

plsrog <- pls_rog(X,Y,0.999)

plsrog <- pls_loading(plsrog)
plsrog$loading$R
plsrog$loading$p.value

### OS-PCA
data(turnover)

plsrog <- os_pca(X,Y,0.999)

plsrog <- pls_loading(plsrog)
plsrog$loading$R
plsrog$loading$p.value


