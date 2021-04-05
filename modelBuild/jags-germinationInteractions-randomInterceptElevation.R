model {
  
  # Likelihood --------------------------------------------------------------
  
  # seed germination --------------------------------------------------------------
  for(i in 1:n){
    
    # LIKELIHOOD
    seedlingJan[i] ~ dbinom(g[i], totalJan[i])
    
    logit(g[i]) <- alpha[siteGermination[i]] + beta[1]*clim.std[i] #+ beta[2]*p.std[i] + beta[3]*t.std[i]*p.std[i]
    
  }
  
  for(i in 1:n_siteGermination){
    alpha[i] ~ dnorm(mu.alpha,tau.alpha)
  }
  
  tau.alpha <- 1/sigma.alpha*sigma.alpha
  sigma.alpha ~ dnorm(0,.001) T(0,)
  mu.alpha ~ dnorm(0,1)
  
  # Priors --------------------------------------------------------------
  for (i in 1){
    beta[i] ~ dnorm(0, 1)  
  }
  
}



