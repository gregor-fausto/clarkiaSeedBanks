model {
  
  # Likelihood --------------------------------------------------------------
  
  # seed germination --------------------------------------------------------------
  for(i in 1:n){
    
    # LIKELIHOOD
    seedlingJan[i] ~ dbinom(g[i], totalJan[i])
    
    logit(g[i]) <- alpha[siteGermination[i]] + beta[1]*clim.std[i]# + beta[2]*p.std[i] + beta[3]*t.std[i]*p.std[i]
    
  }
  
  # Weakly informative prior (Rosenbaum et al. 2019)
  #sigma_g[i] ~ dnorm(0,1) T(0,)
  #tau_g[i] <- 1/(sigma_g[i]*sigma_g[i])
  
  for(i in 1:n_siteGermination){
    alpha[i] ~ dnorm(0,1)
  }
  
  # Priors --------------------------------------------------------------
  for (i in 1){
    beta[i] ~ dnorm(0, 1)  
  }
  
}



