# Regularized logistic function

corDesc_logit <- function(y, x, z = NULL, standardize = TRUE, external = TRUE, 
                          alpha = 1, alpha1 = 0, alpha2 = 1, thresh = 1e-5) {
    n <- NROW(y)
    p <- NCOL(x)
    v <- rep(1/n, n)
    nlam <- 100
    
    # standardize x
    if (standardize == TRUE) {
        xm <- colMeans(x)
        x_norm <- sweep(x, 2, xm, "-")
        xs <- drop(sqrt(crossprod(v, x_norm^2)))
        x_norm <- sweep(x_norm, 2, xs, "/")
    } else {
        x_norm <- x
    }
    
    if (external == FALSE) {
    # initialize beta, weights, working response
        b_prime <- rep(0, p + 1)
        betaHats <- mat.or.vec(nr=nlam, nc=p)
        intercept <- numeric()
        
        p1 <- mean(y)
        w <- rep(p1*(1 - p1), n)
        r <- (y-p1)/n
        intercept_prime <- log(p1/(1 - p1))
        b_prime[1] <- intercept_prime
        
        # compute penalty path
        if (alpha > 0) {
            lambdaMax <- max( abs(crossprod(x_norm, r)) ) / alpha
        } else if (alpha == 0) {
            lambdaMax <- 1000 * max( abs(crossprod(x_norm, r)) )
        }
        if (n >= p) {
            lambdaMin <- 0.0001*lambdaMax
            logLambda <- seq(log(lambdaMax), log(lambdaMin), length.out = 100)
            lambdas <- exp(logLambda)
        } else {
            lambdaMin <- 0.01*lambdaMax
            logLambda <- seq(log(lambdaMax), log(lambdaMin), length.out = 100)
            lambdas <- exp(logLambda)
        }
        
        # calculate null deviance
        logLikelihood <- function(y, x, intercept, beta) {
            p0 <- 1/(1 + exp(drop(x %*% beta[-1]) + beta[1]))
            p <- ifelse(y==1, (1-p0), p0)
            ll <- sum(log(p))
            return(ll)
        }
        deviance_null <- -2*logLikelihood(y,x,intercept,b_prime)
        
        b_current <- b_prime
        
        # penalty path loop
        for(i in 1:nlam) {
            lambda_current <- lambdas[i] 
            iter <- 0
            
            # outer reweighted least square loop
            converge_outer <- FALSE
            while(!converge_outer) {
                b_old <- b_current
                x2w <- drop(crossprod(x_norm^2, w/n))
                
                converge_inner <- FALSE
                while (!converge_inner) {
                    iter <- iter + 1
                    dev_in <- 0.0
                    # inner coordinate descent loop
                    for(j in 1:p) {
                        bj <- b_current[j+1]
                        rj <- r + w*(x_norm[,j]*b_current[j+1])/n
                        wls <- sum(x_norm[,j]*rj)
                        arg <- abs(wls) - alpha*lambda_current
                        if (arg > 0.0) {
                            b_current[j+1] <- sign(wls)*arg / (x2w[j] + (1-alpha)*lambda_current)
                        } else {
                            b_current[j+1] <- 0
                        }
                        del <- b_current[j+1] - bj
                        if (abs(del) > 0.0) {
                            dev_in <- max(dev_in, abs(del))
                            r <- r - del*(x_norm[,j])*w/n
                        }
                    }
                    # update intercept
                    del <- sum(r) / sum(w/n)
                    if (abs(del) > 0.0) {
                        b_current[1] <- b_current[1] + del
                        dev_in <- max(dev_in, abs(del))
                        r <- r - del * (w/n)
                    }
                    if (dev_in < thresh*deviance_null) {converge_inner <- TRUE}
                } # inner while loop
                
                # update weights and working responses
                xbeta <- b_current[1] + drop(x_norm %*% b_current[-1])
                m <- ifelse(xbeta>0, xbeta, 0)
                p1 <- exp(xbeta-m) / (exp(0-m) + exp(xbeta-m))
                w <- p1*(1-p1)
                
                if (sum(w/n) < 1e-9) {
                    converge_outer <- TRUE
                } else {
                    r <- (y-p1)/n
                    dev_out <- max( abs(b_current-b_old) )
                    if (dev_out < thresh*deviance_null) {converge_outer <- TRUE}
                }
            } # outer while loop
            
            # unstanrdardize
            if (standardize == TRUE) {
                betaHats[i, ] <- b_current[-1]/xs
                intercept[i] <- b_current[1] - drop(crossprod(b_current[-1]/xs, xm))
            } else {
                betaHats[i, ] <- b_current[-1]
                intercept[i] <- b_current[1]
            }
        }
        
        return(list(coef=cbind(intercept, betaHats), lambda=lambdas))
    }
    
    # Incorporating external information
    if (external == TRUE) {
        q <- ncol(z)
        # standardize xz
        if (standardize == TRUE) {
            xz <- x %*% z
            xzm <- colMeans(xz)
            xz <- sweep(xz, 2, xzm, "-")
            xzs <- drop(sqrt(crossprod(v, xz^2)))
            xz <- sweep(xz, 2, xzs, "/")
        } else {
            xz <- x %*% z
        }
        
        g_prime <- rep(0, p+1)
        a_prime <- rep(0, q)
        nlam <- 20
        betaHats <- array(0, dim = c(1+p+q, nlam, nlam))
        
        p1 <- mean(y)
        w <- rep(p1*(1 - p1), n)
        r_l11 <- (y-p1)/n
        intercept_prime <- log(p1/(1 - p1))
        g_prime[1] <- intercept_prime
        
        # compute penalty path
        if (alpha2 > 0) {
            lambda2_max <- max( abs(crossprod(xz, r_l11)) ) / alpha2
        } else if (alpha2 == 0) {
            lambda2_max <- 1000 * max( abs(crossprod(xz, r_l11)) ) 
        }
        if (alpha1 > 0 ) {
            lambda1_max <- max( abs(crossprod(x_norm, r_l11)) ) / alpha1
        } else if (alpha1 == 0) {
            lambda1_max <- 1000 * max( abs(crossprod(x_norm, r_l11)) ) 
        }
        if (n >= (p+q)) {
            lambda2_min <- 0.0001*lambda2_max
            lambda1_min <- 0.0001*lambda1_max
            lambda2 <- exp( seq(log(lambda2_max), log(lambda2_min), length.out = nlam) )
            lambda1 <- exp( seq(log(lambda1_max), log(lambda1_min), length.out = nlam) )
        } else {
            lambda2_min <- 0.01*lambda2_max
            lambda1_min <- 0.01*lambda1_max
            lambda2 <- exp( seq(log(lambda2_max), log(lambda2_min), length.out = nlam) )
            lambda1 <- exp( seq(log(lambda1_max), log(lambda1_min), length.out = nlam) )
        }
        
        # calculate null deviance
        logLikelihood <- function(y, x, intercept, beta) {
            p0 <- 1/(1 + exp(drop(x %*% beta[-1]) + beta[1]))
            p <- ifelse(y==1, (1-p0), p0)
            ll <- sum(log(p))
            return(ll)
        }
        deviance_null <- -2*logLikelihood(y,x,intercept_prime,g_prime)
        
        g_current <- g_prime
        a_current <- a_prime
        # penalty path loop
        for (l2 in 1:nlam) {
            lambda2_current <- lambda2[l2]
            r <- r_l11
            
            for (l1 in 1:nlam) {
                lambda1_current <- lambda1[l1]
                dev_out <- 1e10
                iter <- 0
                
                # keep the residual of lambda1 at l1=1 for warm start
                if (l1==2) {r_l11 <- r}
                
                # outer reweighted least square loop
                converge_outer <- FALSE
                while (!converge_outer) {
                    g_old <- g_current
                    a_old <- a_current
                    x2w <- drop(crossprod(x_norm^2, w/n))
                    xz2w <- drop(crossprod(xz^2, w/n))
                    
                    # inner coordinate descent loop
                    converge_inner <- FALSE
                    while (!converge_inner) {
                        iter <- iter + 1
                        dev_in <- 0.0
                        for (j in 1:p) {
                            gj <- g_current[j+1]
                            rj <- r + w*(x_norm[,j]*g_current[j+1])/n
                            wls <- sum(x_norm[,j]*rj)
                            arg <- abs(wls) - alpha1*lambda1_current
                            if (arg > 0.0) {
                                g_current[j+1] <- sign(wls)*arg / (x2w[j] + (1-alpha1)*lambda1_current)
                            } else {
                                g_current[j+1] <- 0
                            }
                            del <- g_current[j+1] - gj
                            if (abs(del) > 0.0) {
                                dev_in <- max(dev_in, abs(del))
                                r <- r - del*(x_norm[,j])*w/n
                            }
                        }
                        for (k in 1:q) {
                            ak <- a_current[k]
                            rk <- r + w*(a_current[k]*xz[,k])/n
                            wls <- sum(xz[,k]*rk)
                            arg <- abs(wls) - alpha2*lambda2_current
                            if(arg > 0.0) {
                                a_current[k] <- sign(wls)*arg / (xz2w[k] + (1-alpha2)*lambda2_current)
                            } else {
                                a_current[k] <- 0
                            }
                            del <- a_current[k] - ak
                            if(abs(del) > 0.0) {
                                dev_in <- max(dev_in, abs(del))
                                r <- r - del*(xz[,k])*w/n
                            }
                        }
                        # update intercept
                        del <- sum(r) / sum(w/n)
                        if (abs(del) > 0.0) {
                            g_current[1] <- g_current[1] + del
                            dev_in <- max(dev_in, abs(del))
                            r <- r - del * (w/n)
                        }
                        if (dev_in < thresh*deviance_null) {converge_inner <- TRUE}
                    } # inner while loop
                    
                    # update weights and working responses
                    xbeta <- g_current[1] + drop(x_norm %*% g_current[-1]) + drop(xz %*% a_current)
                    m <- ifelse(xbeta>0, xbeta, 0)
                    p1 <- exp(xbeta-m) / (exp(0-m) + exp(xbeta-m))
                    w <- p1*(1-p1)
                    
                    if (sum(w/n) < 1e-9) {
                        converge_outer <- TRUE
                    } else {
                        r <- (y-p1)/n
                        dev_out <- max( abs(c(g_current,a_current) - c(g_old,a_old)) )
                        if (dev_out < thresh*deviance_null) {converge_outer <- TRUE}
                    }
                } # outer while loop
                
                if (standardize == TRUE) {
                    betaHats[1,l1,l2] <- g_current[1] - drop(crossprod(g_current[-1]/xs, xm)) - drop(crossprod(a_current/xzs, xzm))
                    betaHats[2:(1+p+q),l1,l2] <- c((g_current[-1] + drop(z %*% a_current)), a_current) / c(xs, xzs) 
                } else {
                    betaHats[1,l1,l2] <- g_current[1]
                    betaHats[2:(1+p+q),l1,l2] <- c(g_current[-1] + drop(z %*% a_current), a_current) 
                }
            } # inner for lambda1 loop
            
            # create estimate names, column names
            est_names <- c("intercept", paste("beta_", 1:p, sep=""), paste("alpha_", 1:q, sep=""))
            dimnames(betaHats)[[1]] <- est_names
            
        } # outer for lambda2 loop
        
        return(list(coef=betaHats, lambda=rbind(lambda1, lambda2)))
    }
}



# Cross Validation
CV_cor_logit <- function(y, x, z, fold = 5, standardize = TRUE, alpha1 = 0, alpha2 = 1) {
    library(hierr)
    library(bigmemory)
    library(pROC)
    
    n <- NROW(y)
    p <- NCOL(x)
    v <- rep(1/n, n)
    nlam <- 20
    p1 <- mean(y)
    r_l11 <- (y-p1)/n
    if (standardize == TRUE) {
        xm <- colMeans(x)
        x_norm <- sweep(x, 2, xm, "-")
        xs <- drop(sqrt(crossprod(v, x_norm^2)))
        x_norm <- sweep(x_norm, 2, xs, "/")
        
        zm <- colMeans(z)
        z_norm <- sweep(z, 2, zm, "-")
        zs <- drop(sqrt(crossprod(rep(1/p,p), z_norm^2)))
        z <- sweep(z_norm, 2, zs, "/")
        
        xz <- x_norm %*% z_norm
        
    } else {
        x_norm <- x
        xz <- x %*% z
    }
    
    # compute penalty path
    if (alpha2 > 0) {
        lambda2_max <- max( abs(crossprod(xz, r_l11)) ) / alpha2
    } else if (alpha2 == 0) {
        lambda2_max <- 1000 * max( abs(crossprod(xz, r_l11)) ) 
    }
    if (alpha1 > 0 ) {
        lambda1_max <- max( abs(crossprod(x_norm, r_l11)) ) / alpha1
    } else if (alpha1 == 0) {
        lambda1_max <- 1000 * max( abs(crossprod(x_norm, r_l11)) ) 
    }
    if (n >= (p+q)) {
        lambda2_min <- 0.0001*lambda2_max
        lambda1_min <- 0.0001*lambda1_max
        lambda2 <- exp( seq(log(lambda2_max), log(lambda2_min), length.out = nlam) )
        lambda1 <- exp( seq(log(lambda1_max), log(lambda1_min), length.out = nlam) )
    } else {
        lambda2_min <- 0.01*lambda2_max
        lambda1_min <- 0.01*lambda1_max
        lambda2 <- exp( seq(log(lambda2_max), log(lambda2_min), length.out = nlam) )
        lambda1 <- exp( seq(log(lambda1_max), log(lambda1_min), length.out = nlam) )
    }
    
    # randomly split data to k folds
    folds <- mat.or.vec(fold, n/fold)
    idx <- 1:n
    for (i in 1:fold) {
        folds[i, ] <- sample(idx, n/fold)
        idx <- idx[-folds[i, ]]
    }
    
    myPenalty <- definePenalty(penalty_type = 0,
                               penalty_type_ext = 1,
                               user_penalty = lambda1,
                               user_penalty_ext = lambda2)
    
    auc_val <- mat.or.vec(20, 20)
    for (i in 1:fold) {
        xtrain <- as.big.matrix(x[-folds[i, ], ])
        ytrain <- y[-folds[i, ]]
        yval <- y[folds[i, ]]
        fit <- hierr(x = xtrain, 
                     y = ytrain,
                     external = z,
                     intercept = c(T, F),
                     standardize = c(T, T),
                     family = "binomial", 
                     penalty = myPenalty, 
                     control = list(tolerance = 1e-7))
        intercept <- fit$beta0
        beta <- fit$betas
        xbeta <- apply(beta, 2:3, function(a) drop(x[folds[i, ], ] %*% a))
        xbeta <- sweep(xbeta, 2:3, intercept, "+")
        pred <- ifelse(1/(1+exp(-xbeta)) > 0.5, 1, 0)
        auc_val <- auc_val + apply(pred, 2:3, function(x, y) auc(y, x), y=yval)
    }
    auc_mean <- auc_val / fold
    lam_inx <- which(auc_mean==max(auc_mean), arr.ind = T)
    
}

# simulation
simLogitExt <- function(N, p, q, true_auc, snr_z, nSim, cov_x = "auto", alpha1 = 0, alpha2 = 1, 
                        nonzero_a = 0.2, external=NULL) {
    # fix alphas
    a <- rep(0, q)
    a[c(1:ceiling(nonzero_a*q))] <- 0.2
    
    if ( is.null(external) ) {
        z <- rbinom(p*q, size = 1, prob = 0.05)
        z <- matrix(z, nrow = p, ncol = q)
    } else {
        z <- external
    }
    
    # simulate betas
    var_za <- var(drop(z %*% a))  # empirical variance
    var_b <- var_za / snr_z
    b <- drop(z %*% a) + sqrt(var_b)*rnorm(p)
    
    # set parameters
    if (cov_x == "auto") {
        corr_x <- 0.5
        sigma <- matrix(NA, nrow = p, ncol = p)
        for (j in 1:p) {
            for (i in 1:p) {
                sigma[i, j] <- corr_x^abs(i - j)
            }
        }
    } else if (cov_x == "cs") {
        block <- matrix( rep(.5, (p/q)^2), nrow=p/q )
        diag(block) <- 1
        sigma <- as.matrix(bdiag(replicate(q, block, simplify=F)))
    }
    
    # set true AUC, generate means by LDA linear regression connection
    f <- function(constant) {
        pnorm(0, mean = drop(crossprod(drop(sigma %*% b)/constant, b)), 
              sd = sqrt(2*drop(crossprod(b, sigma) %*% b)), lower.tail = F) - true_auc
    }
    constant <- uniroot(f, c(.001,1000), tol = 0.0001)$root
    u0 <- rep(0, p)
    u1 <- u0 + drop(sigma %*% b)/constant
    
    N1 <- N/2
    N0 <- N - N1
    b0 <- -1/2 * drop(crossprod((u1 + u0), b)) + log(N1/N0)   # from ESL page 109
    y <- c(rep(1, N1), rep(0, N0))
    
    # create validation set
    y_val <- c(rep(1, N1), rep(0, N0))
    x1_val <- mvrnorm(N1, u1, sigma)
    x0_val <- mvrnorm(N0, u0, sigma)
    x_val <- rbind(x1_val, x0_val)
    
    # create test set
    y_test <- c(rep(1, 5000), rep(0, 5000))
    x1_test <- mvrnorm(5000, u1, sigma)
    x0_test <- mvrnorm(5000, u0, sigma)
    x_test <- rbind(x1_test, x0_test)
    
 
    myPenalty <- define_penalty(penalty_type = alpha1,
                                penalty_type_ext = alpha2)
    # myPenalty_noext <- define_penalty(penalty_type = 0,
    #                                   num_penalty = 100) 
    auc_test <- numeric()
    coefs <- mat.or.vec(nSim, 1+p+q)
    auc_test_noext <- numeric()
    coefs_noext <- mat.or.vec(nSim, 1+p)
    for (j in 1:nSim) {
        x1 <- mvrnorm(N1, u1, sigma)
        x0 <- mvrnorm(N0, u0, sigma)
        x <- rbind(x1, x0)

        fit <- xrnet(x = x, 
                     y = y,
                     external = z,
                     intercept = c(T, F),
                     standardize = c(T, T),
                     family = "binomial", 
                     penalty = myPenalty, 
                     control = list(tolerance = 1e-7, max_iterations = 1e7))

        intercept <- fit$beta0
        beta <- fit$betas
        xbeta <- apply(beta, 2:3, function(a) drop(x_val %*% a))
        xbeta <- sweep(xbeta, 2:3, intercept, "+")
        auc_val <- apply(xbeta, 2:3, function(x, y) auc(y, x), y=y_val)
        idx <- which(auc_val == max(auc_val), arr.ind = TRUE)

        coefs[j, ] <- c(intercept[idx[1,1], idx[1,2]], beta[, idx[1,1], idx[1,2]], fit$alphas[, idx[1,1], idx[1,2]])
        pred <- coefs[j, 1] + drop(x_test %*% coefs[j, 2:(p+1)])
        auc_test[j] <- auc(y_test, pred)
        
        
        fit_glmnet <- glmnet(x = x, y = y, alpha = 0, family = "binomial")
        intercept_noext <- fit_glmnet$a0
        beta_noext <- as.matrix(fit_glmnet$beta)
        xbeta_noext <- apply(beta_noext, 2, function(a) drop(x_val %*% a)) 
        xbeta_noext <- sweep(xbeta_noext, 2, intercept_noext, "+")
        auc_val_noext <- apply(xbeta_noext, 2, function(x, y) auc(y, x), y=y_val)
        idx_noext <- which.max(auc_val_noext)
        
        coefs_noext[j, ] <- c(intercept_noext[idx_noext], beta_noext[, idx_noext])
        pred_noext <- coefs_noext[j, 1] + drop(x_test %*% coefs_noext[j, -1]) 
        auc_test_noext[j] <- auc(y_test, pred_noext)

    }
    return(list(auc_test = auc_test, auc_test_noext = auc_test_noext, 
                coefs = coefs, coefs_noext = coefs_noext))
}



#Regular logistic regression irls method
# irls <- function(y, x, tolerance=1e-5, iterMax=50) {
#     x <- model.matrix(y~x)
#     iter <- 0
#     deviation <- 100
#     betaHat <- rep(0, ncol(x))
#     
#     while(deviation >= tolerance && iter <= iterMax) {
#         betaOld <- betaHat
#         p0 <- 1/(1 + exp(x%*%betaOld))
#         p1 <- 1 - p0
#         W <- diag(as.vector(p1*p0))
#         # z <- x%*%betaOld + solve(W)%*%(y - p1) 
#         gradient <- t(x) %*% (y - p1)
#         hessian <- -t(x) %*% W %*% x
#         betaHat <- betaOld - solve(hessian, gradient)
#         iter <- iter + 1
#         deviation <- sum(abs(betaHat-betaOld))
#     }
# if (iter > iterMax) warning("Maximum number of iterations reached.")
# return(as.vector(betaHat))    
# }
# 
# betas <- irls(y,x[,1:15])
# 
# dat <- data.frame(y=y, x=x)
# x_used <- xnam <- paste("x.", 1:15, sep="")
# fmla <- as.formula(paste("y ~ ", paste(x_used, collapse= "+")))
# fit <- glm(fmla, family="binomial", data=dat)
# coefs <- coef(fit)
# 
# rbind(betas, coefs)






