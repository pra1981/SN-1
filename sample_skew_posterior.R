sample_skew_posterior <- function(y,nsim = 1000,burn = 100)
{
  library(tidyverse)
  library(sn)
  
  x = rep(0,length(y))
  X = cbind(0,x)
  n = length(y)
  
  # Priors
  p <- ncol(X) + 1 # dimensions of X plus and column for Z
  betastar0 <- rep(0,p) # prior mean for beta with psi
  T0 <- diag(0.01,p) # prior precision for beta with psi
  a <- b <- 0.001 # gamma hyper parameters for tau
  
  # Initialize
  tau <- 1 # initial precision
  psi <- 0 # initial psi
  betastar <- rep(0,p) # initial beta0, beta1, psi
  z <- rep(0,n) # empty z storage
  
  # Sample storage
  thin<-1	# Thinning interval
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
      z[j] <- rtruncnorm(n = 1,a = 0,b = Inf,mean = m_j,sd = sqrt(v_j)) # update jth component of z
    }
    
    # Form Xstar, augmented design matrix of X and z
    Xstar <- cbind(X,z)
    
    # Update betastar
    v<-solve(T0+tau*crossprod(Xstar)) # full conditional var of betastar
    m<-v%*%(T0%*%betastar0+tau*crossprod(Xstar,y)) # full conditional mean of betastar
    betastar<-c(rmvnorm(1,m,v)) # draw new betastar
    
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
  
  return(Alpha)
}