rm(list=ls(all=TRUE))

# Turnover data by Nakayama

# -------------------------------------
#   OS-PCA : Orthogonal smoothed PCA
# -------------------------------------
# // file
file <- "C:/R/data/Turnover_nakayama.csv"
DATA0 <- read.csv(file)
DATA <- DATA0[,1:(33+7)]

# // metabolome data
X <- as.matrix(DATA[,-c(1:7)])
X <- t(X) # sample*metabolites
X <- scale(X, scale=TRUE)

# // class label
Z <- c(1:11)
class <- c(as.numeric(matrix(1,11)),as.numeric(matrix(2,11)),as.numeric(matrix(3,11)))
Y0 <- factor(Z)
Y <- model.matrix(~ Y0 + 0)
Y <- scale(Y, scale=FALSE)

# ------------------------
#   differential matrix
# ------------------------
D <- diff(diff(Y)) # Second-order differential

DA <- D
D0 <- matrix(0,nrow(D),ncol(D))
D1 <- cbind(D,D0,D0)
D2 <- cbind(D0,D,D0)
D3 <- cbind(D0,D0,D)
DA <- rbind(D1,D2,D3)

kappa <- 0.999 # smoothing parameter
#kappa <- 0 # PCA

E <- (1-kappa)*diag(1,ncol(X))+kappa*t(X)%*%t(DA)%*%DA%*%X

# ------------
#   OS-PCA
# ------------
G <- chol(solve(E))
W0 <- svd(G%*%t(X)%*%X)$v

# score
t <- X%*%W0

R <- chol(E)
z <- svd(t(X)%*%X%*%solve(R))$v
W2 <- solve(R)%*%z
s <- X%*%W2

# factor loading
rs <- NULL;rt <- NULL;pt <- NULL
for(i in 1:ncol(X)){
  rs <- cor.test(s[,2],X[,i])
  rt[i] <- rs$estimate
  pt[i] <- rs$p.value
}
q <- p.adjust(pt,method="BH")

# ---------
#   plot
# ---------
par(mfrow=c(1,2)) 
#plot(t[,1],t[,2],col=class, pch=class)
#plot(s[,1],s[,2],col=class, pch=class)

plot(t[,1],t[,2], pch=class)
plot(s[,1],s[,2], pch=class)

# ---------------
#   save result
# ---------------
ALL <- data.frame(DATA[,1], DATA[,2], w=W0[,1],R=rt,p=pt)
write.csv(ALL, file="C:/R/all.csv")

R1 <- NULL; P1 <- NULL
for(i in 1:ncol(X)){
  r <- cor.test(t[,1], X[,i])
  R1[i] <- r$estimate
  P1[i] <- r$p.value
}

R2 <- NULL; P2 <- NULL
for(i in 1:ncol(X)){
  r <- cor.test(t[,2], X[,i])
  R2[i] <- r$estimate
  P2[i] <- r$p.value
}

R <- cbind(R1,R2)
P <- cbind(P1,P2)

# -----------------------
#   Contribution Ratio
# -----------------------
l <- NULL
for(i in 1:(nrow(X)-1)){
  l[i] <- cov(t[,i],s[,i])
}
print(100*abs(l)/sum(abs(l)))

lambda <- svd(G%*%t(X)%*%X)$d[1:(nrow(X)-1)]
print(100*lambda/sum(lambda))

           