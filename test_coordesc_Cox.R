library(survival)
library(glmnet)
library(MASS)   # gehan dataset is from MASS package

source("Cox_Est.R")
source("Regularized_Cox.R")
source("fast_regularized_cox.R")
Rcpp::sourceCpp('coxphRcpp.cpp')
Rcpp::sourceCpp('cdl_cox_rcpp.cpp')


cat("Gehan dataset in survival package, with ties", "\n")
gehan$treat <- relevel(gehan$treat, ref="control")
cat("Result from coxph: ", coef(coxph(Surv(time, cens) ~ treat, data=gehan, ties="breslow")), "\n")
y <- gehan[, 2:3]
x <- as.matrix(ifelse(gehan$treat=="control", 0, 1))
cat("Result from Breslow: ", irls_cox(y, x, method = "Breslow"), "\n")
cat("Result from after linearization: ", irls_cox(y, x, method = "Linearized"), "\n")
y <- data.matrix(y, rownames.force = F)
cat("Result from coxphRcpp: ", coxphRcpp(y, x), "\n")


data(CoxExample) # no ties
# x <- x[1:50, 1:10]
# y <- y[1:50,]

# cox model without penalty, survival package
dat <- data.frame(y=y,x=x)
start <- Sys.time()
fit <- coxph(Surv(y.time, y.status) ~ ., data = dat)
end <- Sys.time()
runtime_surv <- end - start
coef_surv <- fit$coefficients

start <- Sys.time()
coef_bres <- irls_cox(y, x, method = "Breslow")
end <- Sys.time()
runtime_bres <- end - start

start <- Sys.time()
coef_linear <- irls_cox(y, x, method = "Linearized")
end <- Sys.time()
runtime_linear <- end - start

start <- Sys.time()
coef_rcpp <- coxphRcpp(y, x)
end <- Sys.time()
runtime_rcpp <- end - start

print(rbind(coef_surv, coef_bres, coef_linear, coef_rcpp))
print(cbind(runtime_surv, runtime_bres, runtime_linear, runtime_rcpp))


# compare with glmnet
start <- Sys.time()
fit_cdlCox <- cdl_Cox(y, x, external = FALSE, thresh = 1e-7)
coef_cdlCox <- fit_cdlCox$coef
end <- Sys.time()
runtime_cdl <- end - start

start <- Sys.time()
fit_corCox <- corDesc_Cox(y, x, external = FALSE, thresh = 1e-7)
coef_corCox <- fit_corCox$coef
end <- Sys.time()
runtime_cd <- end - start

start <- Sys.time()
fit_glmnet <- glmnet(x, y, family = "cox", alpha = 1)
coef_glmnet <- coef(fit_glmnet)
end <- Sys.time()
runtime_glmnet <- end - start

start <- Sys.time()
fit_rcpp <- cdlcoxRcpp(y, x, z=matrix(nrow=0, ncol=0))
coef_rcpp <- fit_rcpp$beta
end <- Sys.time()
runtime_rcpp <- end - start


print(rbind(coef_cdlCox[10,], coef_corCox[10,], coef_glmnet[,10], coef_rcpp[,10]))
print(rbind(coef_cdlCox[49,], coef_corCox[49,], coef_glmnet[,49], coef_rcpp[,49]))
print(cbind(runtime_cdl, runtime_cd, runtime_glmnet, runtime_rcpp))

# print(rbind(fit_glmnet$lambda, fit_corCox$lambda[1:49]))


