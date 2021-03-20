model { 
  
  # viability ---------------------------------------------------------------
  
  for(k in 1:n_siteViab){
    
    # hyperpriors -------------------------------------------------------------
    
    # theta 1
    mu0_g[k] ~  dnorm(0, 0.001)
    sigma0_g[k] ~ dnorm(0, 0.3) T(0,)
    tau0_g[k] <- 1/(sigma0_g[k]*sigma0_g[k])
    
    # theta 2
    mu0_v[k] ~  dnorm(0, 0.001)
    sigma0_v[k] ~ dnorm(0, 0.3) T(0,)
    tau0_v[k] <- 1/(sigma0_v[k]*sigma0_v[k])
    
    for(i in 1:n_yearViab){
      
      
      # priors ------------------------------------------------------------------
      
      # theta 1
      mu_g[k,i] ~ dnorm(mu0_g[k], tau0_g[k])
      sigma_g[k,i] ~ dnorm(0, 0.3) T(0,)
      tau_g[k,i] <- 1/(sigma_g[k,i]*sigma_g[k,i])
      
      # theta 2
      mu_v[k,i] ~ dnorm(mu0_v[k], tau0_v[k])
      sigma_v[k,i] ~ dnorm(0, 0.3) T(0,)
      tau_v[k,i] <- 1/(sigma_v[k,i]*sigma_v[k,i])
    }
    
  }
  
  # likelihood (viability) --------------------------------------------------------------
  
  for(i in 1:n_bag){
    
    # alpha
    alpha_g[i] ~ dnorm(mu_g[siteViab[i],yearViab[i]],tau_g[siteViab[i],yearViab[i]])
    alpha_v[i] ~ dnorm(mu_v[siteViab[i],yearViab[i]],tau_v[siteViab[i],yearViab[i]])
    
    # logit 
    logit(theta_g[i]) <- alpha_g[i]
    logit(theta_v[i]) <- alpha_v[i]
    
    germCount[i] ~ dbinom( theta_g[i] , germStart[i] )
    viabStain[i] ~ dbinom( theta_v[i] , viabStart[i] )
    
  }
  
  # derived quantities --------------------------------------------------------------
  
  for(i in 1:n_bag){

    nu_bag[i] = theta_g[i] + theta_v[i]*(1-theta_g[i])
    
  }
  
  for(i in 1:n_siteBags){
    
    p.i_g[i] ~ dnorm(mu0_g[i],tau0_g[i])
    logit(p_g[i]) <- p.i_g[i]
    
    p.i_v[i] ~ dnorm(mu0_v[i],tau0_v[i])
    logit(p_v[i]) <- p.i_v[i]
    
    nu_1[i] = p_g[i] + p_v[i]*(1-p_g[i])
    
    for(k in 1:n_yearBags){
      
      p.i0_g[i,k] ~ dnorm(mu_g[i,k],tau_g[i,k])
      logit(p0_g[i,k]) <- p.i0_g[i,k]
      
      p.i0_v[i,k] ~ dnorm(mu_v[i,k],tau_v[i,k])
      logit(p0_v[i,k]) <- p.i0_v[i,k]
      
      nu0_1[i,k] = p0_g[i,k] + p0_v[i,k]*(1-p0_g[i,k])
      
    }
  }
  
} 

