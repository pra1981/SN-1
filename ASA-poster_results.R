# Read complete case data
nurt <- read.csv("for_sas_no_miss.csv")

# Fit MLE model 
fit.mle <- selm(zwfl ~ as.factor(fdsec_status_0),data = nurt)

# Extract MLE alpha
alpha.mle <- extractSECdistr(selm(zwfl ~ as.factor(fdsec_status_0),data = nurt))@dp[3]

# Extract MLE betas
betas.mle <- summary(selm(zwfl ~ as.factor(fdsec_status_0),data = nurt))@param.table[1:3,]

