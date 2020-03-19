model { 

  # priors

   vfreq ~ dbeta(1,1)
   vcount ~ dgamma(1,1/20)
   alpha <- vfreq*vcount
   beta <- (1-vfreq)*vcount
   
  # Likelihood
  for(i in 1:N){
    theta[i] ~ dbeta(alpha,beta)
    y[i] ~ dbinom(theta[i], n[i])
  }
    
    p ~ dbeta(alpha,beta)
  
}
