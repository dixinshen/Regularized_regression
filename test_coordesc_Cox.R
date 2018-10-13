library(survival)
library(glmnet)
library(CVXR)
source("/Users/dixinshen/Dropbox/Convex_Regularization/Cox Lasso/Regularized_Cox.R")
load("/Users/dixinshen/Dropbox/Convex_Regularization/Cox Lasso/simCox.RData")

data(CoxExample) # no ties

fit_corCox <- corDesc_Cox(y, x, prior = FALSE)
coef_corCox <- fit_corCox$coef

fit_glmnet <- glmnet(x, y, family = "cox", alpha = 1)
coef_glmnet <- coef(fit_glmnet)

print(rbind(coef_corCox[2,], coef_glmnet[,2]))
print(rbind(coef_corCox[30,], coef_glmnet[,30]))

# print(rbind(fit_glmnet$lambda, fit_corCox$lambda[1:49]))

#test external information model
ext1_corCox <- corDesc_Cox(y1, x1, z=z1, prior = TRUE)
coef_ext1Cox <- ext1_corCox$coef
lambda_grid <- ext1_corCox$lambda
