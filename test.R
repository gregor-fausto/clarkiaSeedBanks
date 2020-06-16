model{
  # priors
# mu ~ dnorm(0, 0.001)
# sigma ~ dunif(0,1.5)
# tau <- 1/(sigma*sigma)
theta ~ dbeta(1,1)
# p ~ dbeta(1, 1)
lambda ~ dgamma(0.01, 0.01)
lambda2 ~ dgamma(0.01, 0.01)

# likelihood
for(i in 1:n){

 # undercount[i] ~ dbern(p)
  seedlingNumber[i] ~ dpois(lambda)
  fruitplNumber[i] ~ dpois(lambda2)

  # u[i] =  seedlingNumber[i]# +  n_unobs[i]

 # alpha ~ dnorm(mu, tau)

  # logit(theta) <- alpha
  # fruitplNumber[i] ~ dbinom(theta, seedlingNumber[i])
  
}

  p = lambda2/lambda

}
