library(Matrix)
library(MASS)
library(smoothmest)

set.seed(2018)

# Simulate regularized Cox model with prior information
N <- 50
p1 <- 4
p2 <- 100
q1 <- 2
q2 <- 20


# Simulate features with MVN distribution, correlation matrix is block diagnal, with q blocks each has compound symmetry structure
block <- matrix( rep(.5, 25), nrow=5 )
diag(block) <- 1
sigma1 <- bdiag(replicate(q1, rbind(c(1,0.5),c(0.5,1)), simplify=F))
sigma2 <- bdiag(replicate(q2, block, simplify=F))
x1 <- mvrnorm(N, mu=rep(0,p1), Sigma=sigma1)
x2 <- mvrnorm(N, mu=rep(0,p2), Sigma=sigma2)

# Simulate alphas with iid double exponential distribution
a1 <- rdoublex(q1, lambda = 1)
a2 <- rdoublex(q2, lambda = 1)
a1[2] <- 0
a2[seq(1,20,2)] <- 0

# Simulate betas with MVN(zalpha, var*I)
z1 <- cbind(c(1,1,0,0),c(0,0,1,1))
z2 <- mat.or.vec(p2,q2)
for(i in 1:ncol(z2)) {
    z2[(5*i-4):(5*i), i] <- 1
}

b1 <- mvrnorm(mu = z1%*%a1, Sigma = diag(p1))
b2 <- mvrnorm(mu = z2%*%a2, Sigma = diag(p2))

# Simulate time and censoring status
xPrime1 <- cbind(x1, x1%*%z1)
xPrime2 <- cbind(x2, x2%*%z2)
theta1 <- c((b1-z1%*%a1), a1)
theta2 <- c((b2-z2%*%a2), a2)
lambdaT1 <- 0.6 #baseline hazard T1
lambdaT2 <- 0.6 # baseline hazard T2
lambdaC <- 0.2 # hazard of censoring
T1 <- rweibull(N, shape=1, scale=(1/lambdaT1)*exp(-xPrime1%*%theta1))
C <- rweibull(N, shape=1, scale=1/lambdaC)
T2 <- rweibull(N, shape=1, scale=(1/lambdaT2)*exp(-xPrime2%*%theta2))

time1 <- pmin(T1, C)
event1 <- as.numeric(time1==T1)
time2 <- pmin(T2, C)
event2 <- as.numeric(time2==T2)

y1 <- cbind(time1, event1)
y2 <- cbind(time2, event2)

save(y1, y2, x1, x2, z1, z2, b1, b2, a1, a2, file = "/Users/dixinshen/Dropbox/Convex_Regularization/Cox Lasso/simCox.RData")

