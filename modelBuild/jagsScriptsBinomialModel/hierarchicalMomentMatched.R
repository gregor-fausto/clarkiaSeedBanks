model { 
  
   mu ~ dnorm(0,.001)
   sigma ~ dt(0, .04, 1)I(0,)
   alpha <- (mu^2-mu^3-mu*sigma^2)/(sigma^2)
   beta <- (mu-2*mu^2+mu^3-sigma^2+mu*sigma^2)/(sigma^2)
  
  # Likelihood
  for(i in 1:N){
    theta[i] ~ dbeta(alpha+0.01,beta+0.01)
    y[i] ~ dbinom(theta[i], n[i])
  }
    
    p ~ dbeta(alpha,beta)

  
}

