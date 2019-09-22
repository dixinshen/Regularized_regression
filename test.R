library(glmnet)
Rcpp::sourceCpp('cdl_cox_rcpp.cpp')
source("fast_regularized_cox.R")

data(CoxExample)
x <- x[1:50, 1:10]
y <- y[1:50,]
ordered <- order(y[,1])
y <- y[ordered, ]
x <- as.matrix(x[ordered, ])
z <- c(1,1,0,0,1,1,1,0,1,1)
z <- cbind(z, ifelse(z==1, 0, 1))

w <- rep(1/50, 50)

xm <- colMeans(x)
x_norm <- sweep(x, 2, xm, "-")
xs <- drop(sqrt(crossprod(w, x_norm^2)))
x_norm <- sweep(x_norm, 2, xs, "/")

zm <- colMeans(z)
z_norm <- sweep(z, 2, zm, "-")
zs <- drop(sqrt(crossprod(rep(1/10, 10), z_norm^2)))
z_norm <- sweep(z, 2, zs, "/")

xz <- x_norm %*% z_norm

xzs <- drop(sqrt(crossprod(w, sweep(x%*%z, 2, colMeans(x%*%z), "-")^2)))
xzs_norm <- drop(sqrt(crossprod(w, sweep(xz, 2, colMeans(xz), "-")^2)))

# xz_norm <- x %*% z
# xzm <- colMeans(xz_norm)
# xz_norm <- sweep(xz_norm, 2, xzm, "-")
# xzs <- drop(sqrt(crossprod(w, xz_norm^2)))
# xz_norm <- sweep(xz_norm, 2, xzs, "/")


fit_glmnet <- glmnet(x, y, family = "cox")
fit_cdl <- cdl_Cox(y, x, external = F)
fit_rcpp <- cdlcoxRcpp(y, x, z = matrix(nrow=0, ncol=0))
rbind(fit_rcpp$beta[,10], fit_cdl$coef[10,], fit_glmnet$beta[,10])
rbind(fit_rcpp$beta[,30], fit_cdl$coef[30,], fit_glmnet$beta[,30])


ext_rcpp <- cdlcoxRcpp(y, x, z)
coef_extrcpp <- ext_rcpp$beta

ext1_cdl <- cdl_Cox(y, x, z, external = TRUE, thresh = 1e-7)
coef_ext1cdl <- ext1_cdl$coef

rbind(coef_extrcpp[,1], coef_ext1cdl[,1,1])
rbind(coef_extrcpp[,43], coef_ext1cdl[,3,3])
rbind(coef_extrcpp[,190], coef_ext1cdl[,10,10])

# dat <- cbind(y, x)
# write.table(dat, file = "/Users/dixinshen/Dropbox/hierr_cox_logit/PlayEigen/CoxExample.csv", 
#             row.names = F, col.names = F, sep = ",")

