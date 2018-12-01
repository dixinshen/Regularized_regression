# Regularized logistic function

corDesc_logit <- function(y, x, z = NULL, external = TRUE, alpha = 1, alpha1 = 0, alpha2 = 1, thresh = 1e-5) {
    n <- NROW(y)
    p <- NCOL(x)
    v <- rep(1/n, n)
    nlam <- 100
    
    # standardize x
    xm <- colMeans(x)
    x_norm <- sweep(x, 2, xm, "-")
    xs <- drop(sqrt(crossprod(v, x_norm^2)))
    x_norm <- sweep(x_norm, 2, xs, "/")
    
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
            
            betaHats[i, ] <- b_current[-1]/xs
            intercept[i] <- b_current[1] - drop(crossprod(b_current[-1]/xs, xm))
        }
        
        return(list(coef=cbind(intercept, betaHats), lambda=lambdas))
    }
    
    # Incorporating external information
    if (external == TRUE) {
        q <- ncol(z)
        # standardize xz
        xz <- x %*% z
        xzm <- colMeans(xz)
        xz <- sweep(xz, 2, xzm, "-")
        xzs <- drop(sqrt(crossprod(v, xz^2)))
        xz <- sweep(xz, 2, xzs, "/")
        
        g_prime <- rep(0, p+1)
        a_prime <- rep(0, q)
        nlam <- 20
        betaHats <- mat.or.vec(nlam^2, p+q)
        intercept <- numeric()
        
        p1 <- mean(y)
        w <- rep(p1*(1 - p1), n)
        r_l11 <- (y-p1)/n
        intercept_prime <- log(p1/(1 - p1))
        g_prime[1] <- intercept_prime
        
        # compute penalty pathm 
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
            betaHats_inner <- mat.or.vec(nlam, p+q)
            intercept_inner <- numeric()
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
                
                intercept_inner[l1] <- g_current[1] - drop(crossprod(g_current[-1]/xs, xm)) - drop(crossprod(a_current/xzs, xzm))
                betaHats_inner[l1, ] <- c((g_current[-1] + drop(z %*% a_current)), a_current) / c(xs, xzs) 
            } # inner for lambda1 loop
            
            betaHats[((l2-1)*20+1):(l2*20), ] <- betaHats_inner 
            intercept <- c(intercept, intercept_inner)
            # create lambda names, row names
            lambda_names <- {}
            for(i in 1:nlam) {
                lambda_names <- c(lambda_names, paste("lambda_", 1:20, "_", i, sep=""))
            }
            row.names(betaHats) <- lambda_names
            # create estimate names, column names
            est_names <- c(paste("beta_", 1:p, sep=""), paste("alpha_", 1:q, sep=""))
            colnames(betaHats) <- est_names
            
        } # outer for lambda2 loop
        
        return(list(coef=cbind(intercept, betaHats), lambda=rbind(lambda1, lambda2)))
    }
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






