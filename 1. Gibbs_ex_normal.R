
## settings
M = 3000 # total number of MCMC samples
m = 500 # burn-in

n = 50 # number of obs.
mu0 = 10 # true mean
sig20 = 25 # true \sigma^2

delta = 0 # hyperparameters for \mu
tau2 = 10
a=0.1 # hyperparameters for \sigma^2
b=0.1


## data generation
y = rnorm(n, mean=mu0, sd=sqrt(sig20)) 


## Gibbs sampler 
theta.mat = matrix(nrow=M ,ncol=2) # MCMC samples for (\mu, \sigma^2)
sigsq = var(y) 	#initial value of \sigma^2
for(nsim in 1:M){
  # 1. \mu ~ full conditional posterior
  mu = rnorm(1, mean=(delta/tau2 + sum(y)/sig2)/(1/tau2 +n/sig2), sd=sqrt(1/(1/tau2 +n/sig2)))
  
  # 2. \sigma^2 ~ full conditional posterior
  sig2 =1/rgamma(1, n/2 + a, 1/2*sum((y-mu)^2) + b)
  
  theta.mat[nsim, ] = c(mu, sig2)
}


## results
par(mfrow=c(1,2))
plot(theta.mat[,1], type="l", cex.lab=1.5, xlab=expression(mu), ylab="")
plot(theta.mat[,2], type="l", cex.lab=1.5, xlab=expression(sigma^2), ylab="")
par(mfrow=c(1,1))

par(mfrow=c(1,2))
plot(density(theta.mat[m:M,1]), cex.lab=1.5, xlab=expression(mu),ylab="marginal posterior",main="");abline(v=mu0, col=2)
plot(density(theta.mat[m:M,2]), cex.lab=1.5, xlab=expression(sigma^2),ylab="marginal posterior", main="");abline(v=sig20, col=2)
par(mfrow=c(1,1))

