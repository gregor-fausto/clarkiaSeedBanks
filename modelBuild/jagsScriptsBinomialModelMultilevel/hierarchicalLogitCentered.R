model { 

  # priors
  
   mu ~ dnorm(0,.001)
   sigma ~ dunif(0,20)
   tau <- 1/(sigma*sigma)
   
  
  # Likelihood
  for(i in 1:N){
    alpha[i] ~ dnorm(mu,tau)
    logit(theta[i]) <- alpha[i]
    y[i] ~ dbinom(theta[i], n[i])
    }
  
    p.i ~ dnorm(mu,tau)
    logit(p) <- p.i

}

