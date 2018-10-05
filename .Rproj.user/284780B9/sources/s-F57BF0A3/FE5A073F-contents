# Skew normal bayesian sampler
#
# Carter Allen
# Fall 2018

library(sn) # for rsn
library(truncnorm) # for rtruncnorm
library(mvtnorm) # for rmvtnorm

# Generate Data
# set.seed(1801)
# Generate Data
set.seed(1801)
n <- 100
x <- rnorm(n)
eps <- rnorm(n)
z <- rtruncnorm(n,0,Inf,0,1)
X <- cbind(1,x,z)
alpha <- -1 # skew
omega <- 2 # scale
delta <- alpha/sqrt(1+alpha^2)
psi <- omega * delta
beta <- c(0,1,psi) # adding psi to beta vector
sig2 <- (omega^2)/(1+alpha^2)
y <- X%*%beta + sqrt(sig2)*eps # skew normal data
X <- X[,-3]

# Priors
p <- ncol(X) + 1 # dimensions of X plus and column for Z
betastar0 <- rep(0,p) # prior mean for beta with psi
T0 <- diag(0.01,p) # prior precision for beta with psi
a <- b <- 0.001 # gamma hyper parameters for tau

# Initialize
tau <- 1 # initial precision
psi <- 0 # initial psi
betastar <- c(0,0,psi) # initial beta0, beta1, psi
z <- rep(0,n) # empty z storage

# Sample storage
nsim<-1000 # Number of MCMC Iterations
thin<-1	# Thinning interval
burn<-0	# Burnin
lastit<-(nsim-burn)/thin # Last stored value
Betastar<-matrix(0,lastit,p) # storage for each betastar
Sigma2<-rep(0,lastit) # storage for each 1/tau
Alpha <- rep(0,lastit) # storage for each alpha
Omega <- rep(0,lastit) # storage for each omega
Resid<-matrix(0,nsim,n)  # Store resids
Dy<-matrix(0,lastit,512) # Store density values for residual density plot
Qy<-matrix(0,lastit,100) # Store quantiles for QQ plot

# Gibbs
tmp<-proc.time()
for(i in 1:nsim)
{
  # update z
  for(j in 1:n)
  {
    v_j <- 1/(tau * psi^2 + 1) # full conditional variance of z
    m_j <- tau * v_j * betastar[3] * (y[j]-t(X[j,]) %*% betastar[-3]) # full conditional mean of z
    z[j] <- rtruncnorm(1,0,Inf,m_j,sqrt(v_j)) # update jth component of z
  }
  
  # Form Xstar, augmented design matrix of X and z
  Xstar <- cbind(X,z)
  
  # Update betastar
  v<-solve(T0+tau*crossprod(Xstar)) # full conditional var of betastar
  m<-v%*%(T0%*%betastar0+tau*crossprod(Xstar,y)) # full conditional mean of betastar
  betastar<-c(rmvnorm(1,m,v)) 
  
  # extract psi
  psi <- betastar[3]
  
  # Update tau
  tau <- rgamma(1,a+n/2,b+crossprod(y-Xstar%*%betastar)/2)
  
  # Store Results
  if (i> burn & i%%thin==0) {
    j<-(i-burn)/thin
    Betastar[j,]<-betastar 
    Sigma2[j]<-1/tau
    Alpha[j] <- betastar[3]*sqrt(tau)
    Omega[j] <- sqrt(1/tau + betastar[3]^2)
    Resid[j,]<-resid<-y-Xstar%*%betastar                    # Raw Resid
    Dy[j,]<-density(resid/sd(resid),from=-5,to=5)$y  # Density of Standardized Resids
    Qy[j,]<-quantile(resid/sd(resid),probs=seq(.001,.999,length=100)) # Quantiles for QQ Plot
  }
  if (i%%100==0) print(i)
}

run.time<-proc.time()-tmp
