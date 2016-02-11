library(MASS)
library(ggplot2)
library(dplyr)
set.seed(1)
theme_set(theme_grey(base_size=24))
distance <- function(x, y=x) outer(x, y, FUN="-")
kernel.se <- function(x, y=x) exp(-distance(x, y)**2)
add.posterior.plot <- function(gp, df., colour) {
  gp +
    geom_point(data=df., color=colour) +
    geom_ribbon(data=df., aes(ymin=f-f.var, ymax=f+f.var), alpha=.3, fill=colour)
}
solve.chol <- function(Q, b) backsolve(Q, backsolve(Q, b, transpose=TRUE))
solve.pivoted.chol <- function(Q, b) {
  pivot <- attr(Q, 'pivot')
  # r <- attr(Q, 'rank')
  backsolve(Q, backsolve(Q, b[pivot], transpose=TRUE))[order(pivot)]
}
N <- 35
M <- 20
sigma <- .4
inputs.inducing <- seq(-.5, 4.5, length.out=M)
train <- data.frame(input=rnorm(N, 2, 2)) %>%
  mutate(f=mvrnorm(1, mu=rep(0, n()), kernel.se(input) + sigma**2 * diag(N)), f.var=NA, method=NA)
inputs.test <- seq(min(train$input), max(train$input), length.out=300)
K <- kernel.se(train$input) + 1e-9 * diag(nrow(train))
K.noisy <- K + sigma**2 * diag(nrow(train))
cond.number <- function(X) {
  eig <- eigen(X, only.values=TRUE)
  eig$values[1] / eig$values[nrow(X)]
}
Kfu <- kernel.se(train$input, inputs.inducing)
Kuu <- kernel.se(inputs.inducing)
Kuu.inv <- ginv(Kuu)
Kstaru <- kernel.se(inputs.test, inputs.inducing)
Qstarstar <- Kstaru %*% Kuu.inv %*% t(Kstaru)
Qff <- Kfu %*% Kuu.inv %*% t(Kfu)
# cond.number(K.chol)
K.chol <- chol(K.noisy)
#
# Test Cholesky decomposition of non-negative definite matrix
R <- chol(K, pivot=TRUE)
rank. <- attr(R, 'rank')
pivot <- attr(R, 'pivot')
oo <- order(pivot)
max(abs(t(R[1:rank.,oo])%*%R[1:rank.,oo] - K))
## solve for RPx
a <- Kfu[,1]
b <- K %*% a
RPx <- backsolve(R, b[pivot], k=rank., transpose=TRUE)
x <- backsolve(R, RPx, k=rank.)
x - a[pivot[1:rank.]]
K %*% x[oo] - a
max(abs(crossprod(R[1:rank.,oo]) - K))
rank.

rank.
R[21:24,21:24]

# max(abs(crossprod(K.chol) - K.noisy))
K.inv.f.chol <- solve.chol(K.chol, train$f)
# max(abs(train$f - K.noisy %*% K.inv.f.chol))

Kstarstar <- kernel.se(inputs.test)
Kstartrain <- kernel.se(inputs.test, train$input)
exact <- data.frame(input=inputs.test) %>%
  mutate(f=Kstartrain %*% solve.chol(K.chol, train$f),
         f.var=pmax(0, diag(Kstarstar - Kstartrain %*% solve.chol(K.chol, t(Kstartrain)))),
         method='exact')
#
# Subset of regressors (notation from Quinonero & Rasmussen (2005), eqn 16b)
Sigma <- ginv(1/sigma**2 * crossprod(Kfu) + Kuu)
Sigma.starstar <- Kstaru %*% Sigma %*% t(Kstaru)
sor <- data.frame(input=inputs.test) %>%
  mutate(f=1 / sigma**2 * Kstaru %*% Sigma %*% t(Kfu) %*% train$f,
         f.var=pmax(0, diag(Sigma.starstar)),
         method='SoR')
#
# Deterministic training conditional (eqn 20b)
dtc <- data.frame(input=inputs.test) %>%
  mutate(f=1 / sigma**2 * Kstaru %*% Sigma %*% t(Kfu) %*% train$f,
         f.var=pmax(0, diag(Kstarstar - Qstarstar + Sigma.starstar)),
         method='DTC')
#
# Fully independent training conditional (eqn 24b)
Lambda.inv.fitc <- 1 / (diag(K) - diag(Qff) + sigma**2)
Sigma.fitc <- ginv(Kuu + t(Kfu) %*% diag(Lambda.inv.fitc) %*% Kfu)
fitc <- data.frame(input=inputs.test) %>%
  mutate(f=Kstaru %*% Sigma.fitc %*% t(Kfu) %*% (Lambda.inv.fitc * train$f),
         f.var=pmax(0, diag(Kstarstar - Qstarstar + Sigma.starstar)),
         method='FITC')
#
# Combine all inferences
posterior <- bind_rows(exact, sor, dtc, fitc)
gp <- ggplot(posterior,
                   aes(x=input, y=f, ymin=f-sqrt(f.var), ymax=f+sqrt(f.var), fill=method, colour=method)) +
  geom_point(data=train, size=5) +
  geom_vline(xintercept=inputs.inducing, alpha=.4, linetype='dashed') +
  geom_line(size=3) +
  geom_ribbon(aes(ymin=f-f.var, ymax=f+f.var), alpha=.3)
print(gp)
