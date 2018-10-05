###################################
# Linear Regression.r             #
# Linear Regression Gibbs Sampler #
# Includes Residual Diagnostics   #
# Outfiles: 1) I:\\Brian\\Bayesian Course\\Lectures\\Lecture Slides\\fig\\resid_density.pdf
#           2) I:\\Brian\\Bayesian Course\\Lectures\\Lecture Slides\\fig\\qq.pdf
#           3) I:\\Brian\\Bayesian Course\\SAS\\Data\\linreg.txt (for SAS)
# July 15, 2018                   #
###################################
library(mvtnorm)
#library(coda)            # For MCMC diagnostics

#################
# Generate data #
#################
set.seed(071518)
n<-1000                  # Sample size
x<-rnorm(n)              # Covariate
X<-cbind(1,x)            # Design matrix
p<-ncol(X)               # Number of regression coefs including intercept
beta<-c(-1,1)            # Reg Coefs
mu<-X%*%beta             # Mean
sigma2<-4                # Error variance
y<-rnorm(n,mu,sqrt(sigma2))
fit<-summary(lm(y~x))    # MLE fit for comparison

##########
# Priors #
##########
beta0<-rep(0,p)          # Prior mean of beta
T0<-diag(0.01,p)         # Prior Precision of beta
a<-b<-0.001              # Gamma hyperparms for tau

##########
# Inits  #
##########
taue<-1                  # Prior error precision

#################
# Store Samples #
#################
nsim<-1000               # Number of MCMC Iterations
thin<-1				          # Thinning interval
burn<-nsim/2	            # Burnin
lastit<-(nsim-burn)/thin	# Last stored value

Beta<-matrix(0,lastit,p)
Sigma2<-rep(0,lastit)
Resid<-matrix(0,nsim,n)  # Store resids
Dy<-matrix(0,lastit,512) # Store density values for residual density plot
Qy<-matrix(0,lastit,100) # Store quantiles for QQ plot

#########
# Gibbs #
#########
tmp<-proc.time()

for (i in 1:nsim){
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
    Resid[j,]<-resid<-y-X%*%beta                     # Raw Resid
    Dy[j,]<-density(resid/sd(resid),from=-5,to=5)$y  # Density of Standardized Resids
    Qy[j,]<-quantile(resid/sd(resid),probs=seq(.001,.999,length=100)) # Quantiles for QQ Plot
    
  }
  if (i%%100==0) print(i)
}
run.time<-proc.time()-tmp # MCMC run time 

###########
# Results #
###########
mbeta<-colMeans(Beta)
sbeta<-apply(Beta,2,sd)
qbeta<-apply(Beta,2,quantile,prob=c(0.025,0.975))
msigma2<-mean(Sigma2)
ssigma2<-sd(Sigma2)
qsigma2<-quantile(Sigma2,prob=c(0.025,0.975))

# Compare MLEs and 
fit
cat("mbeta","\n",mbeta,"\n","sbeta","\n",sbeta,"\n","msigma2","\n",msigma2,"\n")

##########################################  
# Density plot of standardized residuals #
##########################################
# pdf(file="I:\\Brian\\Bayesian Course\\Lectures\\Lecture Slides\\fig\\resids_normal.pdf")
par(mfrow=c(2,1))
resids<-rstandard(lm(y~x))        # Standardized resids
plot(density(resids),col="blue4",xlim=c(-5,5),lwd=2,main="Density Plot of Standardized Residuals", xlab="Quantile")

dy<-colMeans(Dy)
lines(seq(-5,5,length=512),dy,col="darkgreen",lty=2,lwd=2)
legend("topleft",col=c("blue4","darkgreen"),lty=c(1,2),legend=c("MLE","Bayes"))

#####################################  
# QQ Plot of standardized residuals #
#####################################

qx<-qnorm(seq(.001,.999,length=100))
qy<-colMeans(Qy)

qqplot(qx,rstandard(lm(y~x)),col="darkred",lwd=2,ylim=c(-3,3),main="QQ Plot of Standardized Residuals",xlab="Normal Quantile",
       ylab="Sample Quantiles") 
points(qx,qy,col="darkgreen",lwd=2)
abline(0,1,col="blue4",lwd=2)
legend("topleft",col=c("darkred","darkgreen"),pch=1,lwd=2,lty=c(NA,NA),legend=c("MLE","Bayes"))

# dev.off()

##############
# Traceplots #
##############
par(mfrow=c(2,1))

plot(1:lastit,Beta[,2],type="l",col="lightgreen")
abline(h=mbeta[2],col="blue4")

plot(1:lastit,Sigma2,type="l",col="lightgreen")
abline(h=msigma2,col="blue4")

################# 
# Output to SAS #
#################
data<-cbind(y,x)
#write.table(data,"I:\\Brian\\Bayesian Course\\SAS\\Data\\linreg.txt",col.names=F,quote=F,row.names=F)