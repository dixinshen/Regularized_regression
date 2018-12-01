library(survival)
library(glmnet)
library(CVXR)

data(BinomialExample)
x <- x[1:50, 1:4]
y <- y[1:50]

# estimate using glm_binomial
fit <- glm(y ~ x, family = "binomial")
coef_logit <- fit$coefficients; coef_logit

# estimate with CVXR
p <- ncol(x)
x_1 <- cbind(1, x)
beta <- Variable(p+1)
obj <- -sum(y*(x_1 %*% beta)) + sum(logistic(x_1 %*% beta))
prob <- Problem(Minimize(obj))
res <- solve(prob)

cat("Compare glm logistic and CVXR, logistic regression without penalization", "\n")
res_1 <- rbind(coef_logit, as.vector(res$getValue(beta)))
row.names(res_1) <- c("glm_binomial", "CVXR")
print(res_1)
cat("\n")



# estimate using glmnet LASSO
fit_glmnet <- glmnet(x, y, family="binomial", alpha=1, thresh = 1e-15)
lambda <- fit_glmnet$lambda[12]
coef_glmnet <- coef(fit_glmnet)[,12]; coef_glmnet

# estimate using CVXR for LASSO
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
x_1 <- cbind(1, x_norm)

beta <- Variable(p+1)
loss <- ( -sum(y*(x_1 %*% beta)) + sum(logistic(x_1 %*% beta)) )/n
obj <- loss + elastic_reg(beta[-1], lambda, 1)
prob <- Problem(Minimize(obj))
res <- solve(prob)

cat("Compare glmnet logistic with LASSO penalty and CVXR", "\n")
res <- as.vector(res$getValue(beta))
res[1] <- res[1] - drop(crossprod(res[-1]/xs, xm))
res_2 <- rbind(coef_glmnet, res / c(1,xs))
row.names(res_2) <- c("glmnet_LASSO", "CVXR")
print(res_2)
cat("\n")



# test with hierarchical regularized logistic with external information
source("Regularized_Logistic.R")
load("simLogit.RData")

# corDesc_logit function
ext_corLogit <- corDesc_logit(y, x, z=z, external = TRUE, thresh = 1e-15)
coef_extLogit <- ext_corLogit$coef
lambda_grid <- ext_corLogit$lambda
lambda <- lambda_grid[,c(1,3,10,20)]

# CVXR
n <- nrow(x)
p <- ncol(x)
q <- ncol(z)
v <- rep(1/n, n)
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
x_prime <- cbind(1, x_norm, xz)

thetas <- mat.or.vec(ncol(lambda), 1+p+q)
for (i in 1:ncol(lambda)) {
    theta <- Variable(1+p+q)
    loss <- ( -sum(y*(x_prime %*% theta)) + sum(logistic(x_prime %*% theta)) )/n
    obj <- loss + elastic_reg(theta[2:(p+1)], lambda[1,i], 0) + elastic_reg(theta[(p+2):(p+q+1)], lambda[2,i], 1) 
    prob <- Problem(Minimize(obj))
    res <- solve(prob)
    theta <- as.vector(res$getValue(theta))
    theta[1] <- theta[1] - drop(crossprod(theta[2:(p+1)]/xs, xm)) - drop(crossprod(theta[(p+2):(p+q+1)]/xzs, xzm))
    theta[2:(p+1)] <- theta[2:(p+1)] + drop(z %*% theta[(p+2):(p+q+1)])
    theta <- theta / c(1,xs,xzs)
    thetas[i,] <- theta
}

cat("Compare hierarchical regularized logistic regression with CVXR", "\n")
res_3 <- rbind(coef_extLogit[1,],thetas[1,],coef_extLogit[43,],thetas[2,],coef_extLogit[190,],thetas[3,],coef_extLogit[400,],thetas[4,])
row.names(res_3) <- c("hierr_logit_lam_1_1","CVXR_lam_1_1","hierr_logit_lam_3_3","CVXR_lam_3_3",
                      "hierr_logit_lam_10_10","CVXR_lam_10_10", "hierr_logit_lam_20_20", "CVXR_lam_20_20")
print(res_3)
cat("\n")



