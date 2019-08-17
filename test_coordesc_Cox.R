library(survival)
library(glmnet)

source("Regularized_Cox.R")
source("fast_regularized_cox.R")
load("simCox.RData")

data(CoxExample) # no ties
# x <- x[1:50, 1:10]
# y <- y[1:50,]

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

#test external information model
ext1_corCox <- corDesc_Cox(y1, x1, z=z1, external = TRUE)
coef_ext1Cox <- ext1_corCox$coef
lambda_grid <- ext1_corCox$lambda

