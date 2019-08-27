# Regular Cox proportional hazard model

irls_cox <- function(y, x, tolerance = 1e-7, iterMax = 1000, method="Linearized") {
    
    if(method == "Linearized") {
        ordered <- order(y[,1])
        y <- y[ordered, ]
        x <- as.matrix(x[ordered, ])

        n <- NROW(y)
        p <- NCOL(x)
        iter <- 0
        deviation <- 1e10
        betaHat <- rep(0, p)
        
        D <- -1
        d <- numeric()   # number of events at time ti
        idx <- 2
        for (i in 1:n) {
            if (y[i, 2]==1 && y[i, 1]!=D[idx-1]) {
                D[idx] <- y[i, 1]
                d[idx-1] <- 1
                idx <- idx + 1
            } else if (y[i, 2]==1 && y[i, 1]==D[idx-1]) {
                d[idx-2] <- d[idx-2] + 1
            }
        }
        D <- D[-1]   # unique event time
        m <- length(D) # number of unique times
        
        Ck_prime <- 0
        Ck <- numeric()
        Ri <- numeric()
        for (k in 1:n) {
            Ck[k] <- Ck_prime
            for (j in (Ck_prime+1):m) {
                if (D[j] <= y[k,1]) {
                    Ck[k] <- Ck[k] + 1
                    Ri[Ck[k]] <- n - k + 1
                } else {
                    break
                }
                Ck_prime <- Ck[k]
            }
            if (Ck_prime==m) {break}
        }
        Ck[(k+1):n] <- Ck_prime
        Ri <- c(n, Ri)
        Ck <- c(0, Ck)
        
        while(deviation >= tolerance && iter <= iterMax) {
            betaPrime <- betaHat
            eta <- as.vector(x %*% betaPrime)
            exp_eta <- exp(eta)
            sum_exp_eta_prime <- sum(exp_eta)
            sum_exp_eta <- numeric()
            for (i in 1:m) {
                if (Ri[i]==Ri[i+1]) {
                    sum_exp_eta[i] <- sum_exp_eta_prime
                } else {
                    sum_exp_eta[i] <- sum_exp_eta_prime - sum( exp(eta[(n-Ri[i]+1):(n-Ri[i+1])]) )
                    sum_exp_eta_prime <- sum_exp_eta[i]
                }
            }
            u_prime <- 0
            u2_prime <- 0
            W <- numeric()   # weights
            wr <- numeric()  # weights * working response
            for (k in 1:n) {
                if (Ck[k+1] - Ck[k] == 0) {
                    W[k] <- exp_eta[k] * u_prime - exp_eta[k] * exp_eta[k] * u2_prime 
                    wr[k] <- W[k] * eta[k] + y[k, 2] - exp_eta[k] * u_prime
                } else {
                    u_prime <- u_prime + d[Ck[k+1]] / sum_exp_eta[Ck[k+1]]
                    u2_prime <- u2_prime + d[Ck[k+1]] / (sum_exp_eta[Ck[k+1]])^2
                    W[k] <- exp_eta[k] * u_prime - exp_eta[k] * exp_eta[k] * u2_prime 
                    wr[k] <- W[k] * eta[k] + y[k, 2] - exp_eta[k] * u_prime
                }
            }
            
            xW <- sweep(t(x), 2, W, "*")
            betaHat <- solve( xW %*% x, t(x) %*% wr )
            iter <- iter + 1
            deviation <- max(abs(betaHat-betaPrime))
        }
        if (iter > iterMax) warning("Maximum number of iterations reached.")
        return(as.vector(betaHat))
    }
    
    # use Breslow's likelihood, handling ties, when no ties, reduce to Cox's partial likelihood
    if(method == "Breslow") {
        n <- NROW(y)
        p <- NCOL(x)
        iter <- 0
        deviation <- 1e10
        betaHat <- rep(0, p)
        
        t_event <- y[,1][y[,2]==1]
        D <- unique(t_event)
        indicator <- function(a,b) {as.numeric(a>=b)}
        M <- outer(y[,1], D, indicator)
        d <- numeric()
        for(i in 1:length(D)) {
            d[i] <- sum(t_event==D[i])
        }
    
        while(deviation >= tolerance && iter <= iterMax) {
            betaPrime <- betaHat
            eta <- x %*% betaPrime
            haz <- numeric()
            cumHaz <- numeric()
            for(i in 1:length(d)) {
                haz[i] <- d[i] / sum( (y[,1] >= D[i])*exp(eta) )
            }
            for(i in 1:n) {
                cumHaz[i] <- sum( (D <= y[i,1])*haz )
            }
            u <- cumHaz*exp(eta)
            #w3m <- sweep(M, 1, exp(eta), "*")
            W2 <- numeric()
            for(i in 1:n) {
                W2[i] <- crossprod( ((exp(eta))^2)[i] * (haz^2/d), M[i,] )
            }
            W <- as.vector(u) - W2
            # W <- diag(as.vector(u)) - w3m %*% diag(haz^2/d) %*% t(w3m)
            # W <- diag(W)
            #w <- eta + solve(W, (y[,2] - u)) #ginv(W) %*% (y[,2] - u)
            wr <- W*eta + y[,2] - u
            xW <- sweep(t(x), 2, W, "*")
            #wStar <- W %*% eta + (y[,2] - u)
            #betaHat <- solve( t(x) %*% W %*% x, t(x) %*% W %*% w )
            betaHat <- solve( xW %*% x, t(x) %*% wr )
            #betaHat <- solve( xW %*% x, xW %*% w)
            #betaHat <- lm.wfit(x, w, W)$coefficients
            iter <- iter + 1
            deviation <- max(abs(betaHat-betaPrime))
        }
        
        if (iter > iterMax) warning("Maximum number of iterations reached.")
        return(as.vector(betaHat))
    } 
    
    # Cox's partial likelihood with no ties
    if(method == "Cox") {
        orders <- order(y[,1])
        y <- y[orders, ]
        x <- x[orders, ]
        n <- NROW(y)
        p <- NCOL(x)
        iter <- 0
        deviation <- 1e10
        betaHat <- rep(0, p)
        
        while(deviation >= tolerance && iter <= iterMax) {
            betaPrime <- betaHat
            eta <- as.vector(x %*% betaPrime)
            haz <- exp(eta)
            rsk <- rev(cumsum(rev(haz)))
            p <- outer(haz, rsk, "/")
            p[upper.tri(p)] <- 0 # this is because at times larger than event time, the event probability of this subject is 0
            gradient <- t(x) %*% (y[,2] - p %*% y[,2])
            
            W <- -p %*% diag(y[,2]) %*% t(p) 
            diag(W) <- diag(p %*% diag(y[,2]) %*% t(1-p))
            hessian <- -t(x) %*% W %*% x
            betaHat <- betaPrime - solve(hessian, gradient)
            iter <- iter + 1
            deviation <- max(abs(betaHat-betaPrime))
        }
        
        if (iter > iterMax) warning("Maximum number of iterations reached.")
        return(as.vector(betaHat))
    }
}



