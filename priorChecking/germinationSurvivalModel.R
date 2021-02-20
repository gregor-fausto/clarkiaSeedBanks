n=10000
# generate Uniform(0,1) prior on germination
g=runif(n,0,1)

# plot Uniform(0,1) histogram
par(mfrow=c(2,1))
hist(g,breaks=100)

# the product of uniform distributions can be calculated as
# https://mathworld.wolfram.com/UniformProductDistribution.html
par(mfrow=c(3,1))
hist(1-g,breaks=100)
hist((1-g)*(1-g),breaks=100)
hist((1-g)*(1-g)*(1-g),breaks=100)

# uniform priors on the probability of germination at age 1, 2, and 3
# generates a not-uniform prior on the probability that a seed has not-germinated at age 2, 3
f= function(x,n){
  (((-1)^(n-1))/factorial(n-1))*(log(x)^(n-1))
}

par(mfrow=c(3,1))
hist(1-g,breaks=100)
hist((1-g)*(1-g),breaks=100,freq=FALSE)
lines(seq(0,1,.001),f(seq(0,1,.001),2),lwd=2,col='red')
hist((1-g)*(1-g)*(1-g),breaks=100,freq=FALSE)
lines(seq(0,1,.001),f(seq(0,1,.001),3),lwd=2,col='red')

#
s = runif(n,0,1)

# if the probability of survival at each time point is uniform,
# the resulting distribution is a uniform product with n+1
par(mfrow=c(3,1))
hist((1-g)*s,breaks=100)
hist((1-g)*(1-g)*s,breaks=100)
hist((1-g)*(1-g)*(1-g)*s,breaks=100)

## the other consideration is whether the risk of mortality varies with age
## the exponential model assumes it does not.