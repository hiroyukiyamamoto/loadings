rm(list=ls(all=TRUE))

# GreenTea metabolome data by Tsugawa

# -------------------------------------------------------
#   OS-PCA : Orthogonal smoothed PCA for repeated data
# -------------------------------------------------------
# // file
file <- "C:/R/data/GreenTeaMetabolome.csv"
DATA0 <- read.csv(file)
DATA <- DATA0[-1,-1]

X0 <- as.matrix(DATA)
X <- matrix(as.numeric(X0),nrow=nrow(DATA))

X <- t(X) # sample*metabolites
X <- scale(X, scale=TRUE)

T <- prcomp(X)

plot(T$x[,1],T$x[,2], pch=c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10), cex=1.5)

# # // class label
Z <- c(1:10)
Y0 <- factor(Z)
Y <- model.matrix(~ Y0 + 0)
Y <- scale(Y, scale=FALSE)

# ------------------------
#   differential matrix
# ------------------------
D <- diff(diff(Y)) # Second-order differential

kappa <- 0.1 # smoothing parameter
#kappa <- 0 # PCA

M <- NULL
for(i in seq(0,27,3)){
  M <- rbind(M,c(rep(0,i), rep(1,3), rep(0,27-i)))
}
M <- M/3

MX <- scale(M%*%X)

E <- (1-kappa)*diag(1,ncol(X))+kappa*t(MX)%*%t(D)%*%D%*%MX

# ------------
#   OS-PCA
# ------------
G <- chol(solve(E))
W0 <- svd(G%*%t(MX)%*%MX)$v

# score
t <- X%*%W0

R <- chol(E)
z <- svd(t(MX)%*%MX%*%solve(R))$v
W2 <- solve(R)%*%z

s <- MX%*%W2

# factor loading
rs <- NULL;rt <- NULL;pt <- NULL
for(i in 1:ncol(X)){
  rs <- cor.test(s[,1],MX[,i])
  rt[i] <- rs$estimate
  pt[i] <- rs$p.value
}

# ---------
#   plot
# ---------
par(mfrow=c(1,2)) 
plot(t[,1],t[,2], pch=c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10), cex=1.5)
plot(s[,1],s[,2], pch=c(1,2,3,4,5,6,7,8,9,10), cex=1.5)

r <- cor.test(t[,1],c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10))
print(r$estimate)

# ---------------
#   save result
# ---------------
ALL <- data.frame(M=as.character(DATA0[-1,1]),w=W0[,1],R=rt,p=pt, q=p.adjust(pt,method="BH"))
write.csv(ALL, file="C:/R/PC1.csv")

#plot(rt, W0[,1])

# -----------------------
#   Contribution Ratio
# -----------------------
l <- NULL
for(i in 1:(nrow(MX)-1)){
  l[i] <- cov(MX%*%W0[,i],MX%*%W2[,i]) # cov(Mt,Ms)
}
print(100*abs(l)/sum(abs(l)))

lambda <- svd(t(MX)%*%MX%*%solve(R))$d[1:(nrow(MX)-1)]
print(100*lambda/sum(lambda))