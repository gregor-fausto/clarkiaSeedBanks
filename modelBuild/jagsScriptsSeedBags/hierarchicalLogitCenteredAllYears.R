model { 

  # viability ---------------------------------------------------------------
  for(i in 1:n_ageBags){
    for(j in 1:n_siteViab){
    
    # hyperpriors -------------------------------------------------------------
    
    # theta 1
    mu0_g[i,j] ~  dnorm(0, 0.001)
    sigma0_g[i,j] ~ dnorm(0, 0.001) T(0,)
    tau0_g[i,j] <- 1/(sigma0_g[i,j]*sigma0_g[i,j])
    
    # theta 2
    mu0_v[i,j] ~  dnorm(0, 0.001)
    sigma0_v[i,j] ~ dnorm(0, 0.001) T(0,)
    tau0_v[i,j] <- 1/(sigma0_v[i,j]*sigma0_v[i,j])
    
    for(k in 1:n_yearViab){
     

    # priors ------------------------------------------------------------------

      # theta 1
      mu_g[i,j,k] ~ dnorm(mu0_g[i,j], tau0_g[i,j])
      sigma_g[i,j,k] ~ dnorm(0, 0.001) T(0,)
      tau_g[i,j,k] <- 1/(sigma_g[i,j,k]*sigma_g[i,j,k])
      
      # theta 2
      mu_v[i,j,k] ~ dnorm(mu0_v[i,j], tau0_v[i,j])
      sigma_v[i,j,k] ~ dnorm(0, 0.001) T(0,)
      tau_v[i,j,k] <- 1/(sigma_v[i,j,k]*sigma_v[i,j,k])
    }
  }
  }

# seed bags ---------------------------------------------------------------
  for(i in 1:n_ageBags){
    
  for(j in 1:n_siteBags){

  # hyperpriors -------------------------------------------------------------

    # theta 1
    mu0_1[i,j] ~  dnorm(0, 0.001)
    sigma0_1[i,j] ~ dnorm(0, 0.001) T(0,)
    tau0_1[i,j] <- 1/(sigma0_1[i,j]*sigma0_1[i,j])
    
    # theta 2
    mu0_2[i,j] ~  dnorm(0, 0.001)
    sigma0_2[i,j] ~ dnorm(0, 0.001) T(0,)
    tau0_2[i,j] <- 1/(sigma0_2[i,j]*sigma0_2[i,j])
    
    # theta 3
    mu0_3[i,j] ~  dnorm(0, 0.001)
    sigma0_3[i,j] ~ dnorm(0, 0.001) T(0,)
    tau0_3[i,j] <- 1/(sigma0_3[i,j]*sigma0_3[i,j])
    
    for(k in 1:n_yearBags){
   
    # priors ------------------------------------------------------------------
         
      # theta 1
      mu_1[i,j,k] ~ dnorm(mu0_1[i,j], tau0_1[i,j])
      sigma_1[i,j,k] ~ dnorm(0, 0.001) T(0,)
      tau_1[i,j,k] <- 1/(sigma_1[i,j,k]*sigma_1[i,j,k])
      
      # theta 2
      mu_2[i,j,k] ~ dnorm(mu0_2[i,j], tau0_2[i,j])
      sigma_2[i,j,k] ~ dnorm(0, 0.001) T(0,)
      tau_2[i,j,k] <- 1/(sigma_2[i,j,k]*sigma_2[i,j,k])
      
      # theta 3
      mu_3[i,j,k] ~ dnorm(mu0_3[i,j], tau0_3[i,j])
      sigma_3[i,j,k] ~ dnorm(0, 0.001) T(0,)
      tau_3[i,j,k] <- 1/(sigma_3[i,j,k]*sigma_3[i,j,k])
    }
  }
  }
  
  
  # likelihood (viability) --------------------------------------------------------------
  
  for(i in 1:n_bag){

    # alpha
    alpha_g[i] ~ dnorm(mu_g[ageViab[i],siteViab[i],yearViab[i]], tau_g[ageViab[i],siteViab[i],yearViab[i]])
    alpha_v[i] ~ dnorm(mu_v[ageViab[i],siteViab[i],yearViab[i]], tau_v[ageViab[i],siteViab[i],yearViab[i]])
    
    # logit 
    logit(theta_g[i]) <- alpha_g[i]
    logit(theta_v[i]) <- alpha_v[i]
    
    germCount[i] ~ dbinom( theta_g[i] , germStart[i] )
    viabStain[i] ~ dbinom( theta_v[i] , viabStart[i] )
    
  }
  
  # likelihood (seed bags) --------------------------------------------------------------
  
  for(i in 1:n){
    
    # alpha
    alpha_1[i] ~ dnorm(mu_1[ageBags[i],siteBags[i],yearBags[i]],tau_1[ageBags[i],siteBags[i],yearBags[i]])
    alpha_2[i] ~ dnorm(mu_2[ageBags[i],siteBags[i],yearBags[i]],tau_2[ageBags[i],siteBags[i],yearBags[i]])
    alpha_3[i] ~ dnorm(mu_3[ageBags[i],siteBags[i],yearBags[i]],tau_3[ageBags[i],siteBags[i],yearBags[i]])
    
    # logit 
    logit(theta_1[i]) <- alpha_1[i]
    logit(theta_2[i]) <- alpha_2[i]
    logit(theta_3[i]) <- alpha_3[i]
    
    # likelihood
    totalJan[i] ~ dbinom(theta_1[i], seedStart[i])
    seedlingJan[i] ~ dbinom(theta_2[i], totalJan[i])
    intactJan[i] = totalJan[i]-seedlingJan[i]
    intactOct[i] ~ dbinom(theta_3[i], intactJan[i])
    
  }
  
  # derived quantities --------------------------------------------------------------
  
  
  
} 

