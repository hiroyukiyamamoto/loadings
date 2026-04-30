library(loadings)
library(geigen)

data(turnover)

X <- scale(turnover$X, scale = TRUE)
D <- turnover$D

A <- t(X) %*% X
n <- ncol(X)

# --- B ---
kappa <- 0.999
B <- (1 - kappa) * diag(1, n) +
  kappa * t(X) %*% t(D) %*% D %*% X

I <- diag(1, n)
O <- matrix(0, n, n)

# --- 目的関数 ---
L <- rbind(
  cbind(O, A),
  cbind(A, O)
)

# --- 制約（正定値にする） ---
eigB <- eigen(B, symmetric = TRUE)$values
eta_max <- sqrt(min(eigB))

eta <- 0.9 * eta_max   # ←ここが重要
tau <- 0            # 数値安定

R <- rbind(
  cbind((1 + tau) * I, -eta * I),
  cbind(-eta * I, B + tau * I)
)

# --- 固有値問題 ---
Z <- geigen::geigen(L, R, symmetric = TRUE)

lambda <- Z$values
lambda_index <- order(lambda, decreasing = TRUE)

V <- Z$vectors

wx <- V[1:n, lambda_index[1:2]]
wy <- V[(n + 1):(2 * n), lambda_index[1:2]]

# --- 合成ベクトル（理論に合わせる） ---
w <- wx - eta * wy

# --- スコア ---
score_y <- X %*% wy[,1]
score_x <- X %*% wx[,1]
score_w <- X %*% w[,1]

# --- プロット ---
group_col <- c(rep(1,11), rep(2,11), rep(3,11))

plot(score_y[,1], X %*% wy[,2],
     col = group_col, pch = 19,
     main = "wy score")

plot(score_w[,1], X %*% w[,2],
     col = group_col, pch = 19,
     main = "w score")

# --- 相関（wyベース） ---
r <- sapply(1:ncol(X), function(i) {
  cor(score_y, X[,i])
})

# --- 理論的比例対象 ---
v1 <- wx[,1] - eta * wy[,1]

# --- 確認 ---
plot(v1, r, pch = 19)
abline(lm(r ~ v1), col = 2)

cor(r, v1)
r / v1

v <- wx - eta * wy
score_v <- X %*% v
plot(score_v[,1], score_v[,2],
     col = group_col, pch = 19,
     main = "v score")
