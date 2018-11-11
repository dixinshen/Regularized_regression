library(survival)
library(glmnet)
library(CVXR)

data(CoxExample)
x <- x[1:50, 1:3]
y <- y[1:50,]
time <- y[,1]
event <- y[,2]

# estimate using survival
dat <- data.frame(y=y,x=x)
fit <- coxph(Surv(y.time, y.status) ~ ., data = dat)
coef_surv <- fit$coefficients; coef_surv

# estimate with CVXR for unpernalized Cox
t_event <- time[event==1]
d <- outer(t_event, time, function(a,b) {a==b})
M <- outer(t_event, time, function(a,b) {a<=b})

logsum <- function(beta, x, M) {
    logsum <- 0
    for (i in 1:nrow(M)) {
        logsum <- logsum + log_sum_exp((x%*%beta)[M[i,]])
    }
    return (logsum)
}

p <- ncol(x)
beta <- Variable(p)
obj <- sum(d %*% x %*% beta) - logsum(beta, x, M)
prob <- Problem(Maximize(obj))
res <- solve(prob)

print(coef_surv)
print(res$getValue(beta))


# extimate using glmnet
data(CoxExample)
x <- x[1:50, 1:3]
y <- y[1:50,]
time <- y[,1]
event <- y[,2]
t_event <- time[event==1]
d <- outer(t_event, time, function(a,b) {a==b})
M <- outer(t_event, time, function(a,b) {a<=b})

fit_glmnet <- glmnet(x, y, family="cox" ,alpha=1)
lambda <- fit_glmnet$lambda[3]
coef_glmnet <- fit_glmnet$beta[,3]
# test with CVXR for LASSO
elastic_reg <- function(theta, lambda = 0, alpha = 0) {
    ridge <- (1 - alpha) * sum(theta^2) / 2
    lasso <- alpha * p_norm(theta, 1)
    lambda * (lasso + ridge)
}
p <- ncol(x)
n <- nrow(x)
v <- rep(1/n, n)
xm <- colMeans(x)
x_norm <- sweep(x, 2, xm, "-")
xs <- drop(sqrt(crossprod(v, x_norm^2)))
x_norm <- sweep(x_norm, 2, xs, "/")
beta <- Variable(p)
loss <- (sum(d %*% x_norm %*% beta) - logsum(beta, x_norm, M))/n
obj <- loss - elastic_reg(beta, lambda, 1)
prob <- Problem(Maximize(obj))
res <- solve(prob)

print(coef_glmnet)
print(res$getValue(beta))


# test with regularized cox with external information
source("/Users/dixinshen/Dropbox/Convex_Regularization/Cox Lasso/Regularized_Cox.R")
load("/Users/dixinshen/Dropbox/Convex_Regularization/Cox Lasso/simCox.RData")

# corDesc_Cox function
ext1_corCox <- corDesc_Cox(y1, x1, z=z1, prior = TRUE)
coef_ext1Cox <- ext1_corCox$coef
lambda_grid <- ext1_corCox$lambda
lambda <- lambda_grid[,3]

# CVXR 
x <- x1
y <- y1
z <- z1
n <- nrow(y)
p <- ncol(x)
q <- ncol(z)
v <- rep(1/n, n)
lambda <- lambda_grid[,3]
xm <- colMeans(x)
x_norm <- sweep(x, 2, xm, "-")
xs <- drop(sqrt(crossprod(v, x_norm^2)))
x_norm <- sweep(x_norm, 2, xs, "/")
q <- ncol(z)
xz <- x %*% z
xzm <- colMeans(xz)
xz <- sweep(xz, 2, xzm, "-")
xzs <- drop(sqrt(crossprod(v, xz^2)))
xz <- sweep(xz, 2, xzs, "/")

x_prime <- cbind(x_norm, xz)
theta <- Variable(nrow(x_prime))
time <- y[,1]
event <- y[,2]
t_event <- time[event==1]
d <- outer(t_event, time, function(a,b) {a==b})
M <- outer(t_event, time, function(a,b) {a<=b})

theta <- Variable(p+q)
loss <- (1/n)*sum(d %*% x_prime %*% theta) - (1/n)*logsum(theta, x_prime, M)

obj <- loss - elastic_reg(theta[1:p], lambda[1], 0) - elastic_reg(theta[(p+1):(p+q)], lambda[2], 1) 
prob <- Problem(Maximize(obj))
res <- solve(prob)
theta <- as.vector(res$getValue(theta))
theta[1:p] <- theta[1:p] + drop(z %*% theta[(p+1):(p+q)])

print(coef_ext1Cox[43,])
print(theta)
