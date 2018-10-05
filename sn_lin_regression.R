# Skew normal bayesian sampler
#
# Carter Allen
# Fall 2018

# For rsn()
# library(sn)
# library(ggplot2)
library(truncnorm)
library(mvtnorm)

# Generate Data
set.seed(1801)
n <- 1000
x <- rnorm(n)
eps <- rnorm(n)
z <- rtruncnorm(n,0,Inf,0,1)
X <- cbind(1,x,z)
p <- ncol(X)
alpha <- -1 # skew
omega <- 2 # scale
delta <- alpha/sqrt(1+alpha^2)
psi <- omega * delta
beta <- c(0,1,psi) # adding psi to beta vector
sig2 <- (omega^2)/(1+alpha^2)
y <- X%*%beta + sqrt(sig2)*eps # skew normal data

# Priors
beta0<-rep(0,p) # Prior mean of beta
T0<-diag(0.01,p) # Prior Precision of beta
a<-b<-0.001 # Gamma hyperparms for tau
taue<-1 # Prior error precision

# Sample storage
nsim<-1000 # Number of MCMC Iterations
thin<-1	# Thinning interval
burn<-nsim/2	# Burnin
lastit<-(nsim-burn)/thin # Last stored value
Beta<-matrix(0,lastit,p)
Sigma2<-rep(0,lastit)
Alpha <- rep(0,lastit)
Omega <- rep(0,lastit)
Resid<-matrix(0,nsim,n)  # Store resids
Dy<-matrix(0,lastit,512) # Store density values for residual density plot
Qy<-matrix(0,lastit,100) # Store quantiles for QQ plot

# Gibbs sampler
tmp<-proc.time()

for(i in 1:nsim)
{
  # Update beta
  v<-solve(T0+taue*crossprod(X))
  m<-v%*%(T0%*%beta0+taue*crossprod(X,y))
  beta<-c(rmvnorm(1,m,v))
  
  # Update tau
  taue<-rgamma(1,a+n/2,b+crossprod(y-X%*%beta)/2)
  
  #################
  # Store Results #
  #################
  if (i> burn & i%%thin==0) {
    j<-(i-burn)/thin
    Beta[j,]<-beta 
    Sigma2[j]<-1/taue
    Alpha[j] <- beta[3]*sqrt(taue)
    Omega[j] <- sqrt(1/taue + beta[3]^2)
    Resid[j,]<-resid<-y-X%*%beta                     # Raw Resid
    Dy[j,]<-density(resid/sd(resid),from=-5,to=5)$y  # Density of Standardized Resids
    Qy[j,]<-quantile(resid/sd(resid),probs=seq(.001,.999,length=100)) # Quantiles for QQ Plot
    
  }
  if (i%%100==0) print(i)
}
run.time<-proc.time()-tmp

# Results
mbeta<-colMeans(Beta)
sbeta<-apply(Beta,2,sd)
qbeta<-apply(Beta,2,quantile,prob=c(0.025,0.975))
msigma2<-mean(Sigma2)
ssigma2<-sd(Sigma2)
qsigma2<-quantile(Sigma2,prob=c(0.025,0.975))