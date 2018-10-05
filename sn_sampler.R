# Sample from SN using stochastic representation
library(truncnorm)

# Parameters
n <- 1000
xi <- 0 # location
alpha <- -2 # skew
omega <- 2 # scale
delta <- alpha/sqrt(1+alpha^2)
psi <- omega * delta
sig2 <- (omega^2)/(1+alpha^2)

# Samples 
z <- rtruncnorm(n,0,Inf,0,1)
eps <- rnorm(n)
y <- xi + psi*z + sqrt(sig2)*eps # y ~ SN(xi,omega^2,alpha)