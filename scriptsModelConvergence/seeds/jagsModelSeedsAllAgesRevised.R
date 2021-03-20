model {
  
  # seed bags ---------------------------------------------------------------
  
  for(j in 1:n_siteBags){
    
    for(k in 1:n_ageBags){
      
      # hyperpriors -------------------------------------------------------------
      
      # theta 1
      mu0_1[j,k] ~  dnorm(0, .5)
      sigma0_1[j,k] ~ dnorm(0, 1) T(0,)
      tau0_1[j,k] <- 1/(sigma0_1[j,k]*sigma0_1[j,k])
      
      # theta 2
      mu0_2[j,k] ~  dnorm(0, .5)
      sigma0_2[j,k] ~ dnorm(0, 1) T(0,)
      tau0_2[j,k] <- 1/(sigma0_2[j,k]*sigma0_2[j,k])
      
      # theta 3
      mu0_3[j,k] ~  dnorm(0, .5)
      sigma0_3[j,k] ~ dnorm(0, 1) T(0,)
      tau0_3[j,k] <- 1/(sigma0_3[j,k]*sigma0_3[j,k])
      
      # theta 4
      mu0_4[j,k] ~  dnorm(0, .5)
      sigma0_4[j,k] ~ dnorm(0, 1) T(0,)
      tau0_4[j,k] <- 1/(sigma0_4[j,k]*sigma0_4[j,k])
    }
    
    for(l in 1:3){
      
      # priors ------------------------------------------------------------------
      
      # theta 1
      mu_1[j,l] ~ dnorm(mu0_1[j,1], tau0_1[j,1])
      sigma_1[j,l] ~ dnorm(0, 1) T(0,)
      tau_1[j,l] <- 1/(sigma_1[j,l]*sigma_1[j,l])
      
      # theta 2
      mu_2[j,l] ~ dnorm(mu0_2[j,1], tau0_2[j,1])
      sigma_2[j,l] ~ dnorm(0, 1) T(0,)
      tau_2[j,l] <- 1/(sigma_2[j,l]*sigma_2[j,l])
      
      # theta 3
      mu_3[j,l] ~ dnorm(mu0_3[j,1], tau0_3[j,1])
      sigma_3[j,l] ~ dnorm(0, 1) T(0,)
      tau_3[j,l] <- 1/(sigma_3[j,l]*sigma_3[j,l])
      
      # theta 4
      mu_4[j,l] ~ dnorm(mu0_4[j,1], tau0_4[j,1])
      sigma_4[j,l] ~ dnorm(0, 1) T(0,)
      tau_4[j,l] <- 1/(sigma_4[j,l]*sigma_4[j,l])
      
    }
    
    for(l in 4:5){
      
      # priors ------------------------------------------------------------------
      
      # theta 1
      mu_1[j,l] ~ dnorm(mu0_1[j,2], tau0_1[j,2])
      sigma_1[j,l] ~ dnorm(0, 1) T(0,)
      tau_1[j,l] <- 1/(sigma_1[j,l]*sigma_1[j,l])
      
      # theta 2
      mu_2[j,l] ~ dnorm(mu0_2[j,2], tau0_2[j,2])
      sigma_2[j,l] ~ dnorm(0, 1) T(0,)
      tau_2[j,l] <- 1/(sigma_2[j,l]*sigma_2[j,l])
      
      # theta 3
      mu_3[j,l] ~ dnorm(mu0_3[j,2], tau0_3[j,2])
      sigma_3[j,l] ~ dnorm(0, 1) T(0,)
      tau_3[j,l] <- 1/(sigma_3[j,l]*sigma_3[j,l])
      
      # theta 4
      mu_4[j,l] ~ dnorm(mu0_4[j,2], tau0_4[j,2])
      sigma_4[j,l] ~ dnorm(0, 1) T(0,)
      tau_4[j,l] <- 1/(sigma_4[j,l]*sigma_4[j,l])
      
    }
    
    for(l in 6){
      
      # priors ------------------------------------------------------------------
      
      # theta 1
      mu_1[j,l] ~ dnorm(mu0_1[j,3], tau0_1[j,3])
      sigma_1[j,l] ~ dnorm(0, 1) T(0,)
      tau_1[j,l] <- 1/(sigma_1[j,l]*sigma_1[j,l])
      
      # theta 2
      mu_2[j,l] ~ dnorm(mu0_2[j,3], tau0_2[j,3])
      sigma_2[j,l] ~ dnorm(0, 1) T(0,)
      tau_2[j,l] <- 1/(sigma_2[j,l]*sigma_2[j,l])
      
      # theta 3
      mu_3[j,l] ~ dnorm(mu0_3[j,3], tau0_3[j,3])
      sigma_3[j,l] ~ dnorm(0, 1) T(0,)
      tau_3[j,l] <- 1/(sigma_3[j,l]*sigma_3[j,l])
      
      # theta 4
      mu_4[j,l] ~ dnorm(mu0_4[j,3], tau0_4[j,3])
      sigma_4[j,l] ~ dnorm(0, 1) T(0,)
      tau_4[j,l] <- 1/(sigma_4[j,l]*sigma_4[j,l])
      
    }
    
  }
  
  # likelihood (seed bags) --------------------------------------------------------------
  
  for(i in 1:n){
    
    # alpha
    alpha_1[i] ~ dnorm(mu_1[siteBags[i],roundBags[i]],tau_1[siteBags[i],roundBags[i]])
    alpha_2[i] ~ dnorm(mu_2[siteBags[i],roundBags[i]],tau_2[siteBags[i],roundBags[i]])
    alpha_3[i] ~ dnorm(mu_3[siteBags[i],roundBags[i]],tau_3[siteBags[i],roundBags[i]])
    alpha_4[i] ~ dnorm(mu_4[siteBags[i],roundBags[i]],tau_4[siteBags[i],roundBags[i]])
    
    # logit
    logit(theta_1[i]) <- alpha_1[i]
    logit(theta_2[i]) <- alpha_2[i]
    logit(theta_3[i]) <- alpha_3[i]
    logit(theta_4[i]) <- alpha_4[i]
    
    # likelihood
    totalJan[i] ~ dbinom(theta_1[i], seedStart[i])
    seedlingJan[i] ~ dbinom(theta_2[i], totalJan[i])
    intactJan[i] = totalJan[i]-seedlingJan[i]
    intactOct[i] ~ dbinom(theta_3[i], intactJan[i])
    totalJan[i] ~ dbinom(p0_1[siteBags[i],roundBags[i]]*p0_3[siteBags[i],roundBags[i]]*theta_4[i], seedStart[i])
    
  }

  
  # derived quantities --------------------------------------------------------------
  
  for(j in 1:n_siteBags){
    
    for(k in 1:n_ageBags){
      
      p.i_1[j,k] ~ dnorm(mu0_1[j,k],tau0_1[j,k])
      logit(p_1[j,k]) <- p.i_1[j,k]
      
      p.i_2[j,k] ~ dnorm(mu0_2[j,k],tau0_2[j,k])
      logit(p_2[j,k]) <- p.i_2[j,k]
      
      p.i_3[j,k] ~ dnorm(mu0_3[j,k],tau0_3[j,k])
      logit(p_3[j,k]) <- p.i_3[j,k]
      
      p.i_4[j,k] ~ dnorm(mu0_4[j,k],tau0_4[j,k])
      logit(p_4[j,k]) <- p.i_4[j,k]
      
    }
    
    for(l in 1:6){
      
      p.i0_1[j,l] ~ dnorm(mu_1[j,l],tau_1[j,l])
      logit(p0_1[j,l]) <- p.i0_1[j,l]
      
      p.i0_2[j,l] ~ dnorm(mu_2[j,l],tau_2[j,l])
      logit(p0_2[j,l]) <- p.i0_2[j,l]
      
      p.i0_3[j,l] ~ dnorm(mu_3[j,l],tau_3[j,l])
      logit(p0_3[j,l]) <- p.i0_3[j,l]
      
      p.i0_4[j,l] ~ dnorm(mu_4[j,l],tau_4[j,l])
      logit(p0_4[j,l]) <- p.i0_4[j,l]
      
    }
    
    
  }
  
}
