library(survival)
library(glmnet)
library(MASS)   # gehan dataset is from MASS package

source("Cox_Est.R")
source("Regularized_Cox.R")
source("fast_regularized_cox.R")


cat("Gehan dataset in survival package, with ties", "\n")
gehan$treat <- relevel(gehan$treat, ref="control")
cat("Result from coxph: ", coef(coxph(Surv(time, cens) ~ treat, data=gehan, ties="breslow")), "\n")
y <- gehan[, 2:3]
x <- as.matrix(ifelse(gehan$treat=="control", 0, 1))
cat("Result from Breslow: ", irls_cox(y, x, method = "Breslow"), "\n")
cat("Result from after linearization: ", irls_cox(y, x, method = "Linearized"), "\n")


data(CoxExample) # no ties
# x <- x[1:50, 1:10]
# y <- y[1:50,]

# cox model without penalty, survival package
start <- Sys.time()
coef_bres <- irls_cox(y, x, method = "Breslow")
end <- Sys.time()
runtime_bres <- end - start

dat <- data.frame(y=y,x=x)
start <- Sys.time()
fit <- coxph(Surv(y.time, y.status) ~ ., data = dat)
end <- Sys.time()
runtime_surv <- end - start
coef_surv <- fit$coefficients

start <- Sys.time()
coef_linear <- irls_cox(y, x, method = "Linearized")
end <- Sys.time()
runtime_linear <- end - start

print(rbind(coef_bres, coef_surv, coef_linear))
print(cbind(runtime_bres, runtime_surv, runtime_linear))


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


print(rbind(coef_cdlCox[2,], coef_corCox[2,], coef_glmnet[,2]))
print(rbind(coef_cdlCox[30,], coef_corCox[30,], coef_glmnet[,30]))
print(cbind(runtime_cdl, runtime_cd, runtime_glmnet))

# print(rbind(fit_glmnet$lambda, fit_corCox$lambda[1:49]))


