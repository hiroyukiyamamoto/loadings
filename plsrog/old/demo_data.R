rm(list=ls(all=TRUE))

X <- matrix(rnorm(9*100),ncol=9)
Y <- c(1,1,1,2,2,2,3,3,3)

Z <- list(X,Y)
names(Z) <- c("X","Y")

save(Z,file="C:/R/demo.Rdata")

