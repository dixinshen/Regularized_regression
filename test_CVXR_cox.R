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
fit <- coxph(Surv(y[,1],y[,2]) ~ x)
coef_surv <- fit$coefficients; coef_surv

# estimate with CVXR for unpenalized Cox
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

cat("Compare survival package and CVXR, Cox model without penalization", "\n")
res_1 <- rbind(coef_surv, as.vector(res$getValue(beta)))
row.names(res_1) <- c("Survival_coxph", "CVXR")
print(res_1)
cat("\n")


# extimate using glmnet
data(CoxExample)
x <- x[1:50, 1:3]
y <- y[1:50,]
time <- y[,1]
event <- y[,2]
t_event <- time[event==1]
d <- outer(t_event, time, function(a,b) {a==b})
M <- outer(t_event, time, function(a,b) {a<=b})

fit_glmnet <- glmnet(x, y, family="cox", alpha=1, thresh = 1e-15)
lambda <- fit_glmnet$lambda[3]
coef_glmnet <- coef(fit_glmnet)[,3]
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

cat("Compare glmnet Cox model with LASSO penalty and CVXR", "\n")
res_2 <- rbind(coef_glmnet, as.vector(res$getValue(beta) / xs))
row.names(res_2) <- c("glmnet_LASSO", "CVXR")
print(res_2)
cat("\n")


# test with regularized cox with external information
source("fast_regularized_cox.R")
Rcpp::sourceCpp('cdl_cox_rcpp.cpp')
load("simCox.RData")

# cdl_Cox, linearized function
start <- Sys.time()
ext1_cdl <- cdl_Cox(y1, x1, z=z1, external = TRUE, thresh = 1e-10)
coef_ext1cdl <- ext1_cdl$coef
end <- Sys.time()
runtime_cdl <- end - start
lambda_grid <- ext1_cdl$lambda
lambda <- lambda_grid[,c(1,3,10,20)]

# cdlcox Rcpp function
start <- Sys.time()
ext1_rcpp <- cdlcoxRcpp(y1, x1, z=z1, thresh = 1e-10)
coef_rcpp <- ext1_rcpp$beta
end <- Sys.time()
runtime_rcpp <- end - start


# CVXR 
x <- x1
y <- y1
z <- z1
n <- nrow(y)
p <- ncol(x)
q <- ncol(z)
w <- rep(1/n, n)
xm <- colMeans(x)
x_norm <- sweep(x, 2, xm, "-")
xs <- drop(sqrt(crossprod(w, x_norm^2)))
x_norm <- sweep(x_norm, 2, xs, "/")
q <- ncol(z)
zm <- colMeans(z)
z_norm <- sweep(z, 2, zm, "-")
zs <- drop(sqrt(crossprod(rep(1/p, p), z_norm^2)))
z_norm <- sweep(z, 2, zs, "/")
xz <- x_norm %*% z_norm
xzs <- drop(sqrt(crossprod(w, sweep(x%*%z, 2, colMeans(x%*%z), "-")^2)))
xzs_norm <- drop(sqrt(crossprod(w, sweep(xz, 2, colMeans(xz), "-")^2)))
# xz <- x %*% z
# xzm <- colMeans(xz)
# xz <- sweep(xz, 2, xzm, "-")
# xzs <- drop(sqrt(crossprod(w, xz^2)))
# xz <- sweep(xz, 2, xzs, "/")

x_prime <- cbind(x_norm, xz)
time <- y[,1]
event <- y[,2]
t_event <- time[event==1]
d <- outer(t_event, time, function(a,b) {a==b})
M <- outer(t_event, time, function(a,b) {a<=b})

thetas <- mat.or.vec(ncol(lambda), p+q)
start <- Sys.time()
for (i in 1:ncol(lambda)) {
    theta <- Variable(p+q)
    loss <- (1/n)*sum(d %*% x_prime %*% theta) - (1/n)*logsum(theta, x_prime, M)
    obj <- loss - elastic_reg(theta[1:p], lambda[1,i], 0) - elastic_reg(theta[(p+1):(p+q)], lambda[2,i], 1) 
    prob <- Problem(Maximize(obj))
    res <- solve(prob)
    theta <- as.vector(res$getValue(theta))
    theta[1:p] <- theta[1:p] + drop(z_norm %*% theta[(p+1):(p+q)])
    theta <- theta / c(xs, xzs/xzs_norm)
    thetas[i,] <- theta
}
end <- Sys.time()
runtime_CVXR <- end - start

cat("Compare hierarchical regularized Cox model with CVXR", "\n")
res_3 <- rbind(coef_ext1cdl[,1,1], coef_rcpp[,1],thetas[1,], coef_ext1cdl[,3,3], coef_rcpp[,43],thetas[2,],
               coef_ext1cdl[,10,10], coef_rcpp[,190],thetas[3,], coef_ext1cdl[,20,20],coef_rcpp[,400],thetas[4,])
row.names(res_3) <- c("cdl_lam_1_1", "cdl_rcpp_lam_1_1","CVXR_lam_1_1", "cdl_lam_3_3", "cdl_rcpp_lam_3_3","CVXR_lam_3_3",
                      "cdl_lam_10_10", "cdl_rcpp_lam_10_10","CVXR_lam_10_10", 
                      "cdl_lam_20_20", "cdl_rcpp_20_20", "CVXR_lam_20_20")
print(res_3)
print(cbind(runtime_cdl, runtime_rcpp, runtime_CVXR))
cat("\n")