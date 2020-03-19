model { 

  # priors
  
   mu ~ dnorm(0,.001)
   sigma ~ dunif(0,20)
  
  # Likelihood
  for(i in 1:N){
    alpha.std[i] ~ dnorm(0,1)
    alpha[i]<- mu + sigma*alpha.std[i]
    logit(theta[i]) <- alpha[i]
    y[i] ~ dbinom(theta[i], n[i])
  }

    tau = 1/(sigma*sigma)
    p.i ~ dnorm(mu,tau)
    logit(p) <- p.i
  
}