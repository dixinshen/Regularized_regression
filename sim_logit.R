library(Matrix)
library(MASS)
library(smoothmest)

set.seed(201811)

N <- 50
p <- 4
q <- 2

sigma <- bdiag(replicate(q, rbind(c(1,0.5),c(0.5,1)), simplify=F))
x <- mvrnorm(N, mu=rep(0,p), Sigma=sigma)
a <- c(1,0)
z <- cbind(c(1,1,0,0),c(0,0,1,1))
b <- mvrnorm(mu = z%*%a, Sigma = diag(p))

sigmoid <- 1 / (1 + exp(-drop(x%*%b)))
y <- ifelse(sigmoid>0.5, 1, 0)



save(x, y, z, a, b, file = "simLogit.RData")
