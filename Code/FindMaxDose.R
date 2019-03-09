# Script for find the posterior of the time until maximum in quadratic regression
# Author: Mattias Villani, Stockholm and Linkoping University. http://mattiasvillani.com

# Loading a package with multivariate normal distribution
#install.packages("mvtnorm")
library(mvtnorm)

# Just some figure settings
#install.packages("RColorBrewer")
library("RColorBrewer")
plotColors = brewer.pal(12, "Paired")
pointColor = plotColors[5] # Color for single dots
lwdDef = 8                 # Default line thickness
lwdThin = 6
lwdThinner = 3
pointSizeDef = 4
cexLabDef = 1.5            # Default scaling of font size labels
cexAxisDef = 1.5           # Default scaling of tick labels


# Defining a function that simulates from the scaled inverse Chi-square distribution
rScaledInvChi2 <- function(n, df, scale){
  return((df*scale)/rchisq(n,df=df))
}

BayesLinReg <- function(y, X, mu_0, Omega_0, v_0, sigma2_0, nIter){
  # Direct sampling from a Gaussian linear regression with conjugate prior:
  #
  # beta | sigma2 ~ N(mu_0, sigma2*inv(Omega_0))
  # sigma2 ~ Inv-Chi2(v_0,sigma2_0)
  # 
  # Author: Mattias Villani, IDA, Linkoping University. http://mattiasvillani.com
  #
  # INPUTS:
  #   y - n-by-1 vector with response data observations
  #   X - n-by-nCovs matrix with covariates, first column should be ones if you want an intercept.
  #   mu_0 - prior mean for beta
  #   Omega_0  - prior precision matrix for beta
  #   v_0      - degrees of freedom in the prior for sigma2
  #   sigma2_0 - location ("best guess") in the prior for sigma2
  #   nIter - Number of samples from the posterior (iterations)
  #
  # OUTPUTS:
  #   results$betaSample     - Posterior sample of beta.     nIter-by-nCovs matrix
  #   results$sigma2Sample   - Posterior sample of sigma2.   nIter-by-1 vector
  
  # Compute posterior hyperparameters
  n = length(y) # Number of observations
  nCovs = dim(X)[2] # Number of covariates
  XX = t(X)%*%X
  betaHat <- solve(XX,t(X)%*%y)
  Omega_n = XX + Omega_0
  mu_n = solve(Omega_n,XX%*%betaHat+Omega_0%*%mu_0)
  v_n = v_0 + n
  sigma2_n = as.numeric((v_0*sigma2_0 + ( t(y)%*%y + t(mu_0)%*%Omega_0%*%mu_0 - t(mu_n)%*%Omega_n%*%mu_n))/v_n)
  invOmega_n = solve(Omega_n)
  
  # The actual sampling
  sigma2Sample = rep(NA, nIter)
  betaSample = matrix(NA, nIter, nCovs)
  for (i in 1:nIter){
    
    # Simulate from p(sigma2 | y, X)
    sigma2 = rScaledInvChi2(n=1, df = v_n, scale = sigma2_n)
    sigma2Sample[i] = sigma2
    
    # Simulate from p(beta | sigma2, y, X)
    beta_ = rmvnorm(n=1, mean = mu_n, sigma = sigma2*invOmega_n)
    betaSample[i,] = beta_
    
  }
  return(results = list(sigma2Sample = sigma2Sample, betaSample=betaSample))
}

# Simuting some data from a quadratic model
set.seed(221)
n = 20 # Number of observations
x = seq(0,60,length = n)
beta0 = 0
beta1 = 10
beta2 = -0.2
sigmaEps = 50
y = beta0 + beta1*x + beta2*x^2 + sigmaEps*rnorm(n)
X = cbind(1,x,x^2)

# Plotting the data
png('FindingMax.png')
par(cex.lab=cexLabDef, cex.axis = cexAxisDef)
par(mfrow = c(2,3))
plot(x,y, xlab = 'time', ylab = 'pain relief', main ='data + fit', 
     axes=FALSE, lwd = 2)
axis(side = 1, at = seq(0, max(x), by = 10))
axis(side = 2, at = seq(-200, 200, by = 100))

# Simuting from the posterior of beta and sigma2
mu_0 = c(0,0,0)
Omega_0 = diag(c(0.01,2,4))
v_0 = 4
sigma2_0 = 0.2
nIter = 10000 # number of simulations from posterior 
simRes = BayesLinReg(y, X, mu_0, Omega_0, v_0, sigma2_0, nIter)

meanSigma = mean(sqrt(simRes$sigma2Sample))
meanBeta0 = mean(simRes$betaSample[,1])
meanBeta1 = mean(simRes$betaSample[,2])
meanBeta2 = mean(simRes$betaSample[,3])

xGrid = seq(min(x),max(x), length = 1000)
fit = meanBeta0 + meanBeta1*xGrid + meanBeta2*xGrid^2
lines(xGrid, fit, col ="red", lwd = lwdThinner)

hist(sqrt(simRes$sigma2Sample), 30, freq = F, xlab = '', ylab = "", main = expression(sigma))
hist(simRes$betaSample[,1], 30, freq = F, xlab = '', ylab = "", main = expression(beta[0]))
hist(simRes$betaSample[,2], 30, freq = F, xlab = '', ylab = "", main =expression(beta[1]))
hist(simRes$betaSample[,3], 30, freq = F, xlab = '', ylab = "", main =expression(beta[2]))

maxPoint = -simRes$betaSample[,2]/(2*simRes$betaSample[,3])
hist(maxPoint, 30, freq = F, xlab = 'time', ylab = "", main ='time until max')

dev.off()

