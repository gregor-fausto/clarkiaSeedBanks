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
    
    # theta 1 - age 2
    mu0_g2[k] ~  dnorm(0, 0.001)
    sigma0_g2[k] ~ dnorm(0, 0.3) T(0,)
    tau0_g2[k] <- 1/(sigma0_g2[k]*sigma0_g2[k])
    
    # theta 2 - age 2
    mu0_v2[k] ~  dnorm(0, 0.001)
    sigma0_v2[k] ~ dnorm(0, 0.3) T(0,)
    tau0_v2[k] <- 1/(sigma0_v2[k]*sigma0_v2[k])
    
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
    
    for(i in 1:n_yearViab2){
      
      # priors ------------------------------------------------------------------
      
      # theta 1
      mu_g2[k,i] ~ dnorm(mu0_g2[k], tau0_g2[k])
      sigma_g2[k,i] ~ dnorm(0, 0.3) T(0,)
      tau_g2[k,i] <- 1/(sigma_g2[k,i]*sigma_g2[k,i])
      
      # theta 2
      mu_v2[k,i] ~ dnorm(mu0_v2[k], tau0_v2[k])
      sigma_v2[k,i] ~ dnorm(0, 0.3) T(0,)
      tau_v2[k,i] <- 1/(sigma_v2[k,i]*sigma_v2[k,i])
    }
  }

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
    
    # theta 6
    mu0_6[k] ~  dnorm(0, 0.001)
    sigma0_6[k] ~ dnorm(0, 0.3) T(0,)
    tau0_6[k] <- 1/(sigma0_6[k]*sigma0_6[k])
    
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
    }
    
    for(i in 1:(n_yearBags2)){
      
      # theta 4
      mu_4[k,i] ~ dnorm(mu0_4[k], tau0_4[k])
      sigma_4[k,i] ~ dnorm(0, 0.3) T(0,)
      tau_4[k,i] <- 1/(sigma_4[k,i]*sigma_4[k,i])
      
      # theta 5
      mu_5[k,i] ~ dnorm(mu0_5[k], tau0_5[k])
      sigma_5[k,i] ~ dnorm(0, 0.3) T(0,)
      tau_5[k,i] <- 1/(sigma_5[k,i]*sigma_5[k,i])
      
      # theta 4
      mu_6[k,i] ~ dnorm(mu0_6[k], tau0_6[k])
      sigma_6[k,i] ~ dnorm(0, 0.3) T(0,)
      tau_6[k,i] <- 1/(sigma_6[k,i]*sigma_6[k,i])
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
  
  for(i in 1:n_bag2){
    
    # alpha
    alpha_g2[i] ~ dnorm(mu_g2[siteViab2[i],yearViab2[i]],tau_g2[siteViab2[i],yearViab2[i]])
    alpha_v2[i] ~ dnorm(mu_v2[siteViab2[i],yearViab2[i]],tau_v2[siteViab2[i],yearViab2[i]])
    
    # logit 
    logit(theta_g2[i]) <- alpha_g2[i]
    logit(theta_v2[i]) <- alpha_v2[i]
    
    germCount2[i] ~ dbinom( theta_g2[i] , germStart2[i] )
    viabStain2[i] ~ dbinom( theta_v2[i] , viabStart2[i] )
    
  }
  
  # likelihood (seed bags) --------------------------------------------------------------
  
  for(i in 1:n){
    
    # alpha
    alpha_1[i] ~ dnorm(mu_1[siteBags[i],yearBags[i]],tau_1[siteBags[i],yearBags[i]])
    alpha_2[i] ~ dnorm(mu_2[siteBags[i],yearBags[i]],tau_2[siteBags[i],yearBags[i]])
    alpha_3[i] ~ dnorm(mu_3[siteBags[i],yearBags[i]],tau_3[siteBags[i],yearBags[i]])

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
  
  for(i in 1:n2){
    
    # alpha
    alpha_4[i] ~ dnorm(mu_4[siteBags2[i],yearBags2[i]],tau_4[siteBags2[i],yearBags2[i]])
    alpha_5[i] ~ dnorm(mu_5[siteBags2[i],yearBags2[i]],tau_5[siteBags2[i],yearBags2[i]])
    alpha_6[i] ~ dnorm(mu_6[siteBags2[i],yearBags2[i]],tau_6[siteBags2[i],yearBags2[i]])
    
    # logit 
    logit(theta_4[i]) <- alpha_4[i]
    logit(theta_5[i]) <- alpha_5[i]
    logit(theta_6[i]) <- alpha_6[i]
    
    # likelihood
    totalJan2[i] ~ dbinom(theta_4[i]*(s1.0[siteBags[i],yearBags[i]]*(1-g1.0[siteBags[i],yearBags[i]])*s2.0[siteBags[i],yearBags[i]]), seedStart2[i])
    seedlingJan2[i] ~ dbinom(theta_5[i], totalJan2[i])
    intactJan2[i] = totalJan2[i]-seedlingJan2[i]
    intactOct2[i] ~ dbinom(theta_6[i], intactJan2[i])
    
  }
  
  # derived quantities --------------------------------------------------------------
  
  for(i in 1:n_siteBags){

    p.i_g[i] ~ dnorm(mu0_g[i],tau0_g[i])
    logit(p_g[i]) <- p.i_g[i]

    p.i_v[i] ~ dnorm(mu0_v[i],tau0_v[i])
    logit(p_v[i]) <- p.i_v[i]

    nu_1[i] = p_g[i] + p_v[i]*(1-p_g[i])

    p.i_1[i] ~ dnorm(mu0_1[i],tau0_1[i])
    logit(p_1[i]) <- p.i_1[i]

    p.i_2[i] ~ dnorm(mu0_2[i],tau0_2[i])
    logit(p_2[i]) <- p.i_2[i]

    p.i_3[i] ~ dnorm(mu0_3[i],tau0_3[i])
    logit(p_3[i]) <- p.i_3[i]

    s1[i] = p_1[i]*(p_2[i] + (1-p_2[i])*(nu_1[i])^(1/3))
    g1[i] = p_2[i]/(1-(1-(nu_1[i]^(1/3)))*(1-p_2[i]))
    s2[i] = p_3[i]*(nu_1[i]^(2/3))


    for(k in 1:n_yearBags){

      p.i0_g[i,k] ~ dnorm(mu_g[i,k],tau_g[i,k])
      logit(p0_g[i,k]) <- p.i0_g[i,k]

      p.i0_v[i,k] ~ dnorm(mu_v[i,k],tau_v[i,k])
      logit(p0_v[i,k]) <- p.i0_v[i,k]

      nu0_1[i,k] = p0_g[i,k] + p0_v[i,k]*(1-p0_g[i,k])

      p.i0_1[i,k] ~ dnorm(mu_1[i,k],tau_1[i,k])
      logit(p0_1[i,k]) <- p.i0_1[i,k]

      p.i0_2[i,k] ~ dnorm(mu_2[i,k],tau_2[i,k])
      logit(p0_2[i,k]) <- p.i0_2[i,k]

      p.i0_3[i,k] ~ dnorm(mu_3[i,k],tau_3[i,k])
      logit(p0_3[i,k]) <- p.i0_3[i,k]

      s1.0[i,k] = p0_1[i,k]*(p0_2[i,k] + (1-p0_2[i,k])*(nu0_1[i,k])^(1/3))
      g1.0[i,k] = p0_2[i,k]/(1-(1-(nu0_1[i,k]^(1/3)))*(1-p0_2[i,k]))
      s2.0[i,k] = p0_3[i,k]*(nu0_1[i,k]^(2/3))


    }
  }
  
} 

