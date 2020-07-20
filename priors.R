# VISUALIZE PRIORS
 
n=10000

##  BINOMIAL

mu0 <- rnorm(n,0,1)
sigma0 <- extraDistr::rhnorm(n,1)
#sigma0 <- runif(n,0,2)
#sigma <- runif(n,0,2)

mu <- rnorm(n,mu0,sigma0)

alpha <- rnorm(n, mean = mu, sd = sigma0)

#r <- rgamma(n ,shape = 0.001, rate = 0.001)
#r <- runif(n, 0, 100)

#lambda <- exp(mu)

par(mfrow=c(2,3))
hist(mu0,breaks=50)
hist(sigma0,breaks=50)

hist(mu,breaks=50)
#hist(r,breaks=50)

#par(mfrow=c(1,1))
plot(density(mu0))
lines(density(mu),lty=2)
 lines(density(alpha),lty=3)
 
theta <- boot::inv.logit(alpha)
hist(theta,breaks=50)

# hist(lambda,breaks=50)
 

## NEGATIVE BINOMIAL

n=10000

mu0 <- rnorm(n,0,1)
sigma0 <- extraDistr::rhnorm(n,.5)

mu <- rnorm(n,mu0,sigma0)

lambda <- exp(mu)

# https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#story-when-the-generic-prior-fails-the-case-of-the-negative-binomial
#r <- rgamma(n ,shape = 0.001, rate = 0.001)
#r <- runif(n, 0, 100)

#lambda <- exp(mu)

par(mfrow=c(2,2))
hist(mu0,breaks=50)
hist(sigma0,breaks=50)

hist(mu,breaks=50)
hist(lambda,breaks=50)


par(mfrow=c(2,2))
k.inv <- extraDistr::rhnorm(n,1)
hist(k.inv,breaks=50)
phi = k.inv^(-1)
hist(phi,breaks=50)

var.nb = lambda + (lambda^2)*k.inv

hist(var.nb,breaks=250)

var.nb = lambda + (lambda^2)/phi
hist(var.nb,breaks=250)

#hist(lambda^2/(var.nb-lambda),breaks=50)


## HIERARCHICAL POISSON

n=10000

mu0 <- rnorm(n,0,1)
sigma0 <- extraDistr::rhnorm(n,1)

mu <- rnorm(n,mu0,sigma0)

#inv.gamma<-truncdist::rtrunc(n,"gamma",a=0.1,shape=0.001,scale=0.001)
inv.gamma<-runif(n,0,1)

# https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#story-when-the-generic-prior-fails-the-case-of-the-negative-binomial

gamma<-exp(mu)

lambda <- rlnorm(n,log(gamma),inv.gamma)

par(mfrow=c(2,2))
hist(mu0,breaks=50)
hist(sigma0,breaks=50)

hist(mu,breaks=50)
hist(lambda,breaks=50)


par(mfrow=c(2,2))
k.inv <- extraDistr::rhnorm(n,1)
hist(k.inv,breaks=50)
phi = k.inv^(-1)
hist(phi,breaks=50)

var.nb = lambda + (lambda^2)*k.inv

hist(var.nb,breaks=250)

var.nb = lambda + (lambda^2)/phi
hist(var.nb,breaks=250)

#hist(lambda^2/(var.nb-lambda),breaks=50)
