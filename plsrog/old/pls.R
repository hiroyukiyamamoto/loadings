rm(list=ls(all=TRUE))

# test

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

# autoscaling
D <- scale(t(X))

# response variable
Y0 <- factor(class)
Y <- model.matrix(~ Y0 + 0)
Y <- scale(Y,scale=FALSE)

# ----------------
#   ordinary PLS
# ----------------
# (sample size)-1
N <- nrow(D)-1

# singular value decomposition
USVx <- svd(t(Y)%*%D/N)
USVy <- svd(t(D)%*%Y/N)

# weight vector matrix
Wx <- USVx$v
Wy <- USVy$v

# score matrix
T <- D%*%Wx
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
#for(i in 1:ncol(D)){
#  R[i] <- cor.test(S[,1],X1[i,])$estimate
#}
