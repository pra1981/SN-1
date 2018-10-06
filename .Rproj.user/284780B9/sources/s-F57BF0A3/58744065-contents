library(tidyverse)
library(sn)

# data
nurt <- read.csv("for_sas_no_miss.csv") %>%
  filter(timepoint == 3) %>%
  mutate(fdsec_status_0_low = if_else(fdsec_status_0 == 2,1,0),
         fdsec_status_0_vlow = if_else(fdsec_status_0 == 3,1,0))
x = cbind(nurt$fdsec_status_0_low,nurt$fdsec_status_0_vlow)
y = nurt$zwfl
n = length(y)

# null model
fit.null <- lm(zwfl ~ fdsec_status_0_low + fdsec_status_0_vlow,data = nurt)

# mle model
fit.mle <- selm(zwfl ~ fdsec_status_0_low + fdsec_status_0_vlow,data = nurt)

# Sampler
X <- cbind(1,x) # design matrix w/ intercept

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
nsim<-11000 # Number of MCMC Iterations
thin<-1	# Thinning interval
burn<-1000	# Burnin
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

theme_set(theme_minimal(base_family = "serif"))
p0 <- ggplot() + 
  geom_histogram(aes(x = rsn(n= 1000,alpha = 2)),fill = "#617AA3",bins = 35) + 
  xlab("Residuals of classic regression model")
p1 <- ggplot() + 
  geom_line(aes(x = 1:lastit, y = Alpha),
            alpha = 0.40) + 
  geom_hline(yintercept = mean(Alpha))
p2 <- ggplot() + 
  geom_line(aes(x = 1:lastit, y = Betastar[,1])) + 
  geom_hline(yintercept = mean(Betastar[,1]))
p3 <- ggplot() + 
  geom_line(aes(x = 1:lastit, y = Betastar[,2])) + 
  geom_hline(yintercept = mean(Betastar[,2]))
p4 <- ggplot() + 
  geom_histogram(aes(x = Alpha)) + 
  xlab("Posterior of Alpha")
p5 <- ggplot() + 
  geom_histogram(aes(x = Betastar[,2])) + 
  xlab("Posterior of low security coefficient")
p6 <- ggplot() + 
  geom_histogram(aes(x = Betastar[,3])) + 
  xlab("Posterior of very low security coefficient")
p7 <- ggplot() + 
  geom_histogram(aes(x = Sigma2)) + 
  xlab("Posterior of Sigma^2")