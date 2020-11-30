model { 
  
  # seed bags ---------------------------------------------------------------
  
  for(k in 1:n_siteBags){
    
    # hyperpriors -------------------------------------------------------------
    
    # theta 1
    mu0_1[k] ~  dnorm(0, 0.001)
    sigma0_1[k] ~ dnorm(0, 0.3) T(0,)
    tau0_1[k] <- 1/(sigma0_1[k]*sigma0_1[k])
    
    # theta 2
    mu0_2[k] ~  dnorm(0, 0.001)
    sigma0_2[k] ~ dnorm(0, 0.3) T(0,)
    tau0_2[k] <- 1/(sigma0_2[k]*sigma0_2[k])
    
    # theta 3
    mu0_3[k] ~  dnorm(0, 0.001)
    sigma0_3[k] ~ dnorm(0, 0.3) T(0,)
    tau0_3[k] <- 1/(sigma0_3[k]*sigma0_3[k])
    
    # theta 4
    mu0_4[k] ~  dnorm(0, 0.001)
    sigma0_4[k] ~ dnorm(0, 0.3) T(0,)
    tau0_4[k] <- 1/(sigma0_4[k]*sigma0_4[k])
    
    # theta 5
    mu0_5[k] ~  dnorm(0, 0.001)
    sigma0_5[k] ~ dnorm(0, 0.3) T(0,)
    tau0_5[k] <- 1/(sigma0_5[k]*sigma0_5[k])

    for(i in 1:n_yearBags){
      
      # priors ------------------------------------------------------------------
      
      # theta 1
      mu_1[k,i] ~ dnorm(mu0_1[k], tau0_1[k])
      sigma_1[k,i] ~ dnorm(0, 0.3) T(0,)
      tau_1[k,i] <- 1/(sigma_1[k,i]*sigma_1[k,i])
      
      # theta 2
      mu_2[k,i] ~ dnorm(mu0_2[k], tau0_2[k])
      sigma_2[k,i] ~ dnorm(0, 0.3) T(0,)
      tau_2[k,i] <- 1/(sigma_2[k,i]*sigma_2[k,i])
      
      # theta 3
      mu_3[k,i] ~ dnorm(mu0_3[k], tau0_3[k])
      sigma_3[k,i] ~ dnorm(0, 0.3) T(0,)
      tau_3[k,i] <- 1/(sigma_3[k,i]*sigma_3[k,i])
      
      # theta 4
      mu_4[k,i] ~ dnorm(mu0_4[k], tau0_4[k])
      sigma_4[k,i] ~ dnorm(0, 0.3) T(0,)
      tau_4[k,i] <- 1/(sigma_4[k,i]*sigma_4[k,i])
      
      # theta 5
      mu_5[k,i] ~ dnorm(mu0_5[k], tau0_5[k])
      sigma_5[k,i] ~ dnorm(0, 0.3) T(0,)
      tau_5[k,i] <- 1/(sigma_5[k,i]*sigma_5[k,i])
      
    }
    
  }
 
  # likelihood (seed bags) --------------------------------------------------------------
  
  for(i in 1:n){
    
    # alpha
    alpha_1[i] ~ dnorm(mu_1[siteBags[i],yearBags[i]],tau_1[siteBags[i],yearBags[i]])
    alpha_2[i] ~ dnorm(mu_2[siteBags[i],yearBags[i]],tau_2[siteBags[i],yearBags[i]])
    alpha_3[i] ~ dnorm(mu_3[siteBags[i],yearBags[i]],tau_3[siteBags[i],yearBags[i]])
    alpha_4[i] ~ dnorm(mu_4[siteBags[i],yearBags[i]],tau_4[siteBags[i],yearBags[i]])
    alpha_5[i] ~ dnorm(mu_5[siteBags[i],yearBags[i]],tau_5[siteBags[i],yearBags[i]])
    
    # logit 
    logit(theta_1[i]) <- alpha_1[i]
    logit(theta_2[i]) <- alpha_2[i]
    logit(theta_3[i]) <- alpha_3[i]
    logit(theta_4[i]) <- alpha_4[i]
    logit(theta_5[i]) <- alpha_5[i]

    # likelihood
    totalJan[i] ~ dbinom(theta_1[i], seedStart[i])
    seedlingJan[i] ~ dbinom(theta_2[i], totalJan[i])
    #seedlingJan[i] ~ dbinom(theta_1[1]*theta_2[i], seedStart[i])
    seedlingJan2[i] ~ dbinom(theta_4[i], seedStart[i])
  
    intactJan[i] = totalJan[i]-seedlingJan[i]
    intactOct[i] ~ dbinom(theta_3[i], intactJan[i])
    #intactOct[i] ~ dbinom(theta_1[1]*(1-theta_2[i])*theta_3[i], seedStart[i])
    intactOct2[i] ~ dbinom(theta_5[i], seedStart[i])
    
    # y_total[i] ~ dbinom(theta_1[i], seedStart[i])
    # y_seedling[i] ~ dbinom(theta_2[i], totalJan[i])
    # y_october[i] ~ dbinom(theta_3[i], intactJan[i])
    # y_seedling2[i] ~ dbinom(theta_4[i], seedStart[i])
    # y_october2[i] ~ dbinom(theta_5[i], seedStart[i])

  }
  
  # derived quantities --------------------------------------------------------------
  
  for(i in 1:n_siteBags){
    
    p.i_1[i] ~ dnorm(mu0_1[i],tau0_1[i])
    logit(p_1[i]) <- p.i_1[i]
    
    p.i_2[i] ~ dnorm(mu0_2[i],tau0_2[i])
    logit(p_2[i]) <- p.i_2[i]
    
    p.i_3[i] ~ dnorm(mu0_3[i],tau0_3[i])
    logit(p_3[i]) <- p.i_3[i]
    
    p.i_4[i] ~ dnorm(mu0_4[i],tau0_4[i])
    logit(p_4[i]) <- p.i_4[i]
    
    p.i_5[i] ~ dnorm(mu0_5[i],tau0_5[i])
    logit(p_5[i]) <- p.i_5[i]
    
    for(k in 1:n_yearBags){
      
      p.i0_1[i,k] ~ dnorm(mu_1[i,k],tau_1[i,k])
      logit(p0_1[i,k]) <- p.i0_1[i,k]
      
      p.i0_2[i,k] ~ dnorm(mu_2[i,k],tau_2[i,k])
      logit(p0_2[i,k]) <- p.i0_2[i,k]
      
      p.i0_3[i,k] ~ dnorm(mu_3[i,k],tau_3[i,k])
      logit(p0_3[i,k]) <- p.i0_3[i,k]
      
      p.i0_4[i,k] ~ dnorm(mu_4[i,k],tau_4[i,k])
      logit(p0_4[i,k]) <- p.i0_4[i,k]
      
      p.i0_5[i,k] ~ dnorm(mu_5[i,k],tau_5[i,k])
      logit(p0_5[i,k]) <- p.i0_5[i,k]
      
    }
  }
 
} 

