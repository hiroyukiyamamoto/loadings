rm(list=ls(all=TRUE))

# -------------
#   Setting
# -------------
load(file="C:/R/demo.Rdata")

X <- Z$X # data
class <- Z$Y # class

# ------------------
#   data preparing
# ------------------
# data matrix
X <- as.matrix(X)
X <- matrix(as.numeric(X),nrow=nrow(X)) # metabolites*samples

# missing value imputation(half of minimum value)
for(i in 1:nrow(X)){
  m <- min(X[i,],na.rm=TRUE)
  na_index <- which(is.na(X[i,]))
  X[i,na_index] <- m/2
}

X1 <- X

# response variable
Y0 <- factor(c(1,1,1,2,2,2,3,3,3))
Y <- model.matrix(~ Y0 + 0)

# penalized matrix
P <- NULL
p <- colSums(Y)
for(i in 1:ncol(Y)){
  P <- cbind(P,Y[,i]/p[i])
}
P <- t(P)

# differential matrix
g <- ncol(Y)
D <- diff(diag(1,g))

# ------------
#   PLS-ROG
# ------------
# autoscaling
X <- scale(t(X))
Y <- scale(Y,scale=FALSE)

# samplesize-1
N <- nrow(X)-1

# smoothing parameter
kappa <- 0.5
C <- kappa*t(Y)%*%t(P)%*%t(D)%*%D%*%P%*%Y+(1-kappa)*t(Y)%*%Y+(10^-10)*diag(1,g)

# cholesky decomposition
Rx <- chol(solve(C))
Ry <- chol(C)

# singular value decomposition
USVx <- svd(Rx%*%t(Y)%*%X/N)
USVy <- svd(t(X)%*%Y%*%solve(Ry)/N)

# weight vector
Wx <- USVx$v
Wy <- solve(Ry)%*%USVy$v

# score
T <- X%*%Wx
S <- Y%*%Wy

# -------------------
#   factor loading
# -------------------
R <- NULL
for(i in 1:(ncol(Y)-1)){
  lambdax <- cov(T[,i],S[,i])
  r <- (sqrt(N)*lambdax*Wx[,i])/(sqrt(t(Wy[,i])%*%t(Y)%*%Y%*%Wy[,i]))
  R <- cbind(R,r)
}

## test
#R <- NULL
#for(i in 1:nrow(X1)){
#  R[i] <- cor.test(S[,1],X1[i,])$estimate
#}