model { 
  
  # viability ---------------------------------------------------------------
  
  for(j in 1:n_siteViab){
    
    for(k in 1:n_ageViab){
      
      # hyperpriors -------------------------------------------------------------
      
      # theta 1
      mu0_g[j,k] ~  dnorm(0, 0.001)
      sigma0_g[j,k] ~ dnorm(0, 0.3) T(0,)
      tau0_g[j,k] <- 1/(sigma0_g[j,k]*sigma0_g[j,k])
      
      # theta 2
      mu0_v[j,k] ~  dnorm(0, 0.001)
      sigma0_v[j,k] ~ dnorm(0, 0.3) T(0,)
      tau0_v[j,k] <- 1/(sigma0_v[j,k]*sigma0_v[j,k])
      
    }
    
    for(l in 1:3){
      
      # priors ------------------------------------------------------------------
      
      # theta 1
      mu_g[j,l] ~ dnorm(mu0_g[j,1], tau0_g[j,1])
      sigma_g[j,l] ~ dnorm(0, 0.3) T(0,)
      tau_g[j,l] <- 1/(sigma_g[j,1]*sigma_g[j,1])
      
      # theta 2
      mu_v[j,l] ~ dnorm(mu0_v[j,1], tau0_v[j,1])
      sigma_v[j,l] ~ dnorm(0, 0.3) T(0,)
      tau_v[j,l] <- 1/(sigma_v[j,1]*sigma_v[j,1])
    }
    
    for(l in 4:5){
      
      # priors ------------------------------------------------------------------
      
      # theta 1
      mu_g[j,l] ~ dnorm(mu0_g[j,2], tau0_g[j,2])
      sigma_g[j,l] ~ dnorm(0, 0.3) T(0,)
      tau_g[j,l] <- 1/(sigma_g[j,2]*sigma_g[j,2])
      
      # theta 2
      mu_v[j,l] ~ dnorm(mu0_v[j,2], tau0_v[j,2])
      sigma_v[j,l] ~ dnorm(0, 0.3) T(0,)
      tau_v[j,l] <- 1/(sigma_v[j,2]*sigma_v[j,2])
    }
    
    for(l in 6){
      
      # priors ------------------------------------------------------------------
      
      # theta 1
      mu_g[j,l] ~ dnorm(mu0_g[j,3], tau0_g[j,3])
      sigma_g[j,l] ~ dnorm(0, 0.3) T(0,)
      tau_g[j,l] <- 1/(sigma_g[j,3]*sigma_g[j,3])
      
      # theta 2
      mu_v[j,l] ~ dnorm(mu0_v[j,3], tau0_v[j,3])
      sigma_v[j,l] ~ dnorm(0, 0.3) T(0,)
      tau_v[j,l] <- 1/(sigma_v[j,3]*sigma_v[j,3])
    }
    
  }
  
  # likelihood (viability) --------------------------------------------------------------
  
  for(i in 1:n_bag){
    
    # alpha
    alpha_g[i] ~ dnorm(mu_g[siteViab[i],round[i]],tau_g[siteViab[i],round[i]])
    alpha_v[i] ~ dnorm(mu_v[siteViab[i],round[i]],tau_v[siteViab[i],round[i]])
    
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

  for(j in 1:n_siteViab){
    
    for(k in 1:n_ageViab){

    p.i_g[j,k] ~ dnorm(mu0_g[j,k],tau0_g[j,k])
    logit(p_g[j,k]) <- p.i_g[j,k]

    p.i_v[j,k] ~ dnorm(mu0_v[j,k],tau0_v[j,k])
    logit(p_v[j,k]) <- p.i_v[j,k]

    nu_1[j,k] = p_g[j,k] + p_v[j,k]*(1-p_g[j,k])

    }
    
    for(l in 1:6){
      
      p.i0_g[j,l] ~ dnorm(mu_g[j,l],tau_g[j,l])
      logit(p0_g[j,l]) <- p.i0_g[j,l]

      p.i0_v[j,l] ~ dnorm(mu_v[j,l],tau_v[j,l])
      logit(p0_v[j,l]) <- p.i0_v[j,l]

      nu0_1[j,l] = p0_g[j,l] + p0_v[j,l]*(1-p0_g[j,l])

    }
    
    # for(l in 4:5){
    #   
    #   p.i0_g[j,l] ~ dnorm(mu_g[j,l],tau_g[j,l])
    #   logit(p0_g[j,l]) <- p.i0_g[j,l]
    #   
    #   p.i0_v[j,l] ~ dnorm(mu_v[j,l],tau_v[j,l])
    #   logit(p0_v[j,l]) <- p.i0_v[j,l]
    #   
    #   nu0_1[j,l] = p0_g[j,l] + p0_v[j,l]*(1-p0_g[j,l])
    #   
    # }
    # 
    # for(l in 6){
    #   
    #   p.i0_g[j,l] ~ dnorm(mu_g[j,l],tau_g[j,l])
    #   logit(p0_g[j,l]) <- p.i0_g[j,l]
    #   
    #   p.i0_v[j,l] ~ dnorm(mu_v[j,l],tau_v[j,l])
    #   logit(p0_v[j,l]) <- p.i0_v[j,l]
    #   
    #   nu0_1[j,l] = p0_g[j,l] + p0_v[j,l]*(1-p0_g[j,l])
    #   
    # }
    
  }

} 

