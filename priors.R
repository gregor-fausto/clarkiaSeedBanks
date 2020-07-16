# VISUALIZE PRIORS
 

## NEGATIVE BINOMIAL

mu0 <- rnorm(n,0,1)
#sigma0 <- extraDistr::rhnorm(n,1)
sigma0 <- runif(n,0,2)
sigma <- runif(n,0,2)

mu <- rnorm(n,mu0,sigma)

#r <- rgamma(n ,shape = 0.001, rate = 0.001)
r <- runif(n, 0, 100)

lambda <- exp(mu)

par(mfrow=c(3,2))
hist(mu0,breaks=50)
hist(sigma0,breaks=50)

hist(mu,breaks=50)
hist(r,breaks=50)

#par(mfrow=c(1,1))
plot(density(mu0))
lines(density(mu),lty=2)
# lines(density(alpha),lty=3)

hist(lambda,breaks=50)