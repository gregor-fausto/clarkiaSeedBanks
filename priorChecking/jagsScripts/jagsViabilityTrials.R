model { 
  
  # viability ---------------------------------------------------------------
  
  for(i in 1:n_siteViab){
    
    # hyperpriors -------------------------------------------------------------
    for(k in 1:2){
      
      # theta 1
      mu0_g[i,k] ~  dnorm(0, 1)
      sigma0_g[i,k] ~ dnorm(0, 1) T(0,)
      tau0_g[i,k] <- 1/(sigma0_g[i,k]*sigma0_g[i,k])
      
      # theta 2
      mu0_v[i,k] ~  dnorm(0, 1)
      sigma0_v[i,k] ~ dnorm(0, 1) T(0,)
      tau0_v[i,k] <- 1/(sigma0_v[i,k]*sigma0_v[i,k])
      
    } 
      
      for(j in 1:3){
        
        # priors ------------------------------------------------------------------
        
        # theta 1
        mu_g[i,j,1] ~ dnorm(mu0_g[i,1], tau0_g[i,1])
        sigma_g[i,j,1] ~ dnorm(0, 1) T(0,)
        tau_g[i,j,1] <- 1/(sigma_g[i,j,1]*sigma_g[i,j,1])
        
        # theta 2
        mu_v[i,j,1] ~ dnorm(mu0_v[i,1], tau0_v[i,1])
        sigma_v[i,j,1] ~ dnorm(0, 1) T(0,)
        tau_v[i,j,1] <- 1/(sigma_v[i,j,1]*sigma_v[i,j,1])
      }
    
    for(j in 1:2){
      
      # priors ------------------------------------------------------------------
      
      # theta 1
      mu_g[i,j,2] ~ dnorm(mu0_g[i,2], tau0_g[i,2])
      sigma_g[i,j,2] ~ dnorm(0, 1) T(0,)
      tau_g[i,j,2] <- 1/(sigma_g[i,j,2]*sigma_g[i,j,2])
      
      # theta 2
      mu_v[i,j,2] ~ dnorm(mu0_v[i,2], tau0_v[i,2])
      sigma_v[i,j,2] ~ dnorm(0, 1) T(0,)
      tau_v[i,j,2] <- 1/(sigma_v[i,j,2]*sigma_v[i,j,2])
    }
    
    for(j in 1){
      
      # priors ------------------------------------------------------------------
      
      # theta 1
    #  mu_g[i,j,3] ~ dnorm(mu0_g[i,3], tau0_g[i,3])
      mu_g[i,j,3] ~ dnorm(0, 1)
      sigma_g[i,j,3] ~ dnorm(0, 1) T(0,)
      tau_g[i,j,3] <- 1/(sigma_g[i,j,3]*sigma_g[i,j,3])
      
      # theta 2
     # mu_v[i,j,3] ~ dnorm(mu0_v[i,3], tau0_v[i,3])
       mu_v[i,j,3] ~ dnorm(0, 1)
      sigma_v[i,j,3] ~ dnorm(0, 1) T(0,)
      tau_v[i,j,3] <- 1/(sigma_v[i,j,3]*sigma_v[i,j,3])
    }
  }
  
  
  # likelihood (viability) --------------------------------------------------------------
  
  for(i in 1:n){
    
    # alpha
    alpha_g[i] ~ dnorm(mu_g[siteViab[i],yearViab[i],ageViab[i]],tau_g[siteViab[i],yearViab[i],ageViab[i]])
    alpha_v[i] ~ dnorm(mu_v[siteViab[i],yearViab[i],ageViab[i]],tau_v[siteViab[i],yearViab[i],ageViab[i]])
    
    # logit 
    logit(theta_g[i]) <- alpha_g[i]
    logit(theta_v[i]) <- alpha_v[i]
    
    germCount[i] ~ dbinom( theta_g[i] , germStart[i] )
    viabStain[i] ~ dbinom( theta_v[i] , viabStart[i] )
    
  }
  
  # derived quantities --------------------------------------------------------------
  
  # for(i in 1:n_siteBags){
  # 
  #   p.i_g[i] ~ dnorm(mu0_g[i],tau0_g[i])
  #   logit(p_g[i]) <- p.i_g[i]
  # 
  #   p.i_v[i] ~ dnorm(mu0_v[i],tau0_v[i])
  #   logit(p_v[i]) <- p.i_v[i]
  #   
  #   nu_1[i] = p_g[i] + p_v[i]*(1-p_g[i])
  #   
  #   
  #   for(k in 1:n_yearBags){
  #     
  #     p.i0_g[i,k] ~ dnorm(mu_g[i,k],tau_g[i,k])
  #     logit(p0_g[i,k]) <- p.i0_g[i,k]
  #     
  #     p.i0_v[i,k] ~ dnorm(mu_v[i,k],tau_v[i,k])
  #     logit(p0_v[i,k]) <- p.i0_v[i,k]
  #     
  #     nu0_1[i,k] = p0_g[i,k] + p0_v[i,k]*(1-p0_g[i,k])
  #     
  #   }
  # }
  # 
  
} 

