model {

  # seed bags ---------------------------------------------------------------
  
  for(i in 1:n_siteSurvival){
    
    # Smits 2015 PNAS: assume the shape is the same but that
    # the scale varies by year
    # saying that a population has a characteristic pattern through time
    # but that the magnitude at which that pattern proceeds varies
    
    # shape parameter: population parameter
    # controls rate of change of age-specific mortality rate
    # change prior to be half-cauchy
    
    a[i] ~ dscaled.gamma(2.5,1)
    #dgamma(0.001,0.001)
    
    # scale parameter: population*year parameter
    #b[i] ~ dgamma(0.001,0.001)
    mu0_s[i] ~ dnorm(0, 0.001)
    sigma0_s[i] ~ dnorm(0, 0.3) T(0,)
    tau0_s[i] <- 1/(sigma0_s[i]*sigma0_s[i])
    
    for(j in 1:n_yearSurvival){
      mu_s[i,j] ~ dnorm(mu0_s[i], tau0_s[i])
      sigma_s[i,j] ~ dnorm(0,0.3) T(0,)
      tau_s[i,j] <- 1/(sigma_s[i,j]*sigma_s[i,j])
    }
    
    # 2nd term of the beta distribution
    for(j in 1:6){
      beta[i,j] ~ dgamma(1,0.001)
    }
    
    # # population*year parameter, hierarchical logit
    # for(k in 1:3){
    #   mu0_g[i,k] ~ dnorm(0, 0.001)
    #   sigma0_g[i,k] ~ dnorm(0, 0.3) T(0,)
    #   tau0_g[i,k] <- 1/(sigma0_g[i,k]*sigma0_g[i,k])
    # }
    # 
    # for(j in refGerm){
    #   # i indexes site; j indexes year; k indexes germination variable
    #   mu_g[i,yearRefGerm[j],indexRefGerm[j]] ~ dnorm(mu0_g[i,indexRefGerm[j]], tau0_g[i,indexRefGerm[j]])
    #   sigma_g[i,yearRefGerm[j],indexRefGerm[j]] ~ dnorm(0,0.3) T(0,)
    #   tau_g[i,yearRefGerm[j],indexRefGerm[j]] <- 1/(sigma_g[i,yearRefGerm[j],indexRefGerm[j]]*sigma_g[i,yearRefGerm[j],indexRefGerm[j]])
    # }
    
  }
  # likelihood (seed bags) --------------------------------------------------------------

  for(i in 1:n1){
    
    # # get probability of germination for each value of index
    # alpha_g1[i] ~ dnorm(mu_g[siteSurvival[i],yearSurvival[i],gIndexSurvival[i]],tau_g[siteSurvival[i],yearSurvival[i],gIndexSurvival[i]])
    # logit(theta_g1[i]) <- alpha_g1[i]
    # 
    # alpha_g2[i] ~ dnorm(mu_g[siteSurvival[i],yearSurvival[i],gIndexSurvival[i]],tau_g[siteSurvival[i],yearSurvival[i],gIndexSurvival[i]])
    # logit(theta_g2[i]) <- alpha_g2[i]
    # 
    # alpha_g3[i] ~ dnorm(mu_g[siteSurvival[i],yearSurvival[i],gIndexSurvival[i]],tau_g[siteSurvival[i],yearSurvival[i],gIndexSurvival[i]])
    # logit(theta_g3[i]) <- alpha_g3[i]

    eta_surv[i] ~ dnorm(mu_s[siteSurvival[i],yearSurvival[i]],tau_s[siteSurvival[i],yearSurvival[i]])
    b[i] <- exp(-(eta_surv[i])/(a[siteSurvival[i]]))

    # mu[i] <- (((1-theta_g1[i])^(g_1[i]))*((1-theta_g2[i])^(g_2[i]))*((1-theta_g3[i])^(g_3[i])))*exp(-(months[i]/b[i])^a[siteSurvival[i]])
    # 
    # # beta is the 2nd term of beta dist.
    # beta_mod[i] <- beta[siteSurvival[i],betaIndex[i]]+.01
    # 
    # # use the beta term to write 1st term of beta dist.
    # alpha_s[i] <- (mu[i]*beta_mod[i])/(1-mu[i])
    # 
    # theta_s[i] ~ dbeta(alpha_s[i]+.01,beta_mod[i]) T(0.001,0.999)
    # 
    # # likelihood
    # y[i] ~ dbinom(theta_s[i], seedStart[i])
    #
    # # log-likelihood
    # logLik_y[i] <- logdensity.bin(y[i],theta_s[i],seedStart[i])
  }

  # for(i in 1:n2){
  #   
  #   alpha_g[i] ~ dnorm(mu_g[siteGermination[i],yearGermination[i],gIndex[i]],tau_g[siteGermination[i],yearGermination[i],gIndex[i]])
  #   
  #   # logit
  #   logit(g[i]) <- alpha_g[i]
  #   
  #   # likelihood
  #   seedlingJan[i] ~ dbinom(g[i], totalJan[i])
  # 
  #   # log-likelihood
  #   logLik_g[i] <- logdensity.bin(seedlingJan[i],g[i] ,totalJan[i])
  # }


  # derived quantities --------------------------------------------------------------



}
