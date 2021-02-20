model {

  # https://github.com/kionaogle/EcoApps_Random_Effects/blob/v1.0/Code%20A3.R
  #https://github.com/CCheCastaldo/CheCastaldo_etal_2019_EcolMongr/blob/master/ModelBuild/GlobalModel/i-009f9c75a9206a0ac/Implementation.R
  # https://github.com/colemonnahan/gradmcmc/tree/v1.0/models/wildflower
  # seed bags ---------------------------------------------------------------
  
  for(i in 1:n_siteSurvival){
    
    # # Smits 2015 PNAS: assume the shape is the same but that
    # # the scale varies by year
    # # saying that a population has a characteristic pattern through time
    # # but that the magnitude at which that pattern proceeds varies
    # 
    # # shape parameter: population parameter
    # # controls rate of change of age-specific mortality rate
    # # change prior to be half-cauchy
    # 
    # a[i] ~ dscaled.gamma(2.5,1)
    # #dgamma(0.001,0.001)
    # 
    # # scale parameter: population*year parameter
    # #b[i] ~ dgamma(0.001,0.001)
    # mu0_s[i] ~ dnorm(0, 0.001)
    # sigma0_s[i] ~ dnorm(0, 0.3) T(0,)
    # tau0_s[i] <- 1/(sigma0_s[i]*sigma0_s[i])
    # 
    # for(j in 1:n_yearSurvival){
    #   mu_s[i,j] ~ dnorm(mu0_s[i], tau0_s[i])
    #   sigma_s[i,j] ~ dnorm(0,0.3) T(0,)
    #   tau_s[i,j] <- 1/(sigma_s[i,j]*sigma_s[i,j])
    # }
    # 
    # # 2nd term of the beta distribution
    # for(j in 1:6){
    #   beta[i,j] ~ dgamma(1,0.001)
    # }
    
    # population*year parameter, hierarchical logit
    for(k in 1:3){
      mu0_g[i,k] ~ dnorm(0, 1)
      
     # sigma0_g[i,k] ~ dt(0,1/.5^2,4) I(0,)
     # sigma0_g[i,k] ~ dunif(0,5)
      
    # sigma0_g[i,k] ~ dt(0, 10^-2, 1)T(0,)
      
      # Weakly informative prior (Rosenbaum et al. 2019)
      sigma0_g[i,k] ~ dnorm(0,1) T(0,)
      tau0_g[i,k] <- 1/(sigma0_g[i,k]*sigma0_g[i,k])
 
      mu0_pred[i,k] ~ dnorm(0,1)
      sigma0_pred[i,k] ~ dnorm(0,1) T(0,)
      tau0_pred[i,k] <- 1/(sigma0_pred[i,k]*sigma0_pred[i,k])     
      # sigma_g[i,yearRefGerm[j],indexRefGerm[j]] ~ dnorm(0,2) T(0,)
      # tau_g[i,yearRefGerm[j],indexRefGerm[j]] <- 1/(sigma_g[i,yearRefGerm[j],indexRefGerm[j]]*sigma_g[i,yearRefGerm[j],indexRefGerm[j]])
      # 
      
    }

    for(j in refGerm){
      # i indexes site; j indexes year; k indexes germination variable
      mu_g[i,yearRefGerm[j],indexRefGerm[j]] ~ dnorm(mu0_g[i,indexRefGerm[j]], tau0_g[i,indexRefGerm[j]])

      mu_pred[i,yearRefGerm[j],indexRefGerm[j]] ~ dnorm(mu0_pred[i,indexRefGerm[j]], tau0_pred[i,indexRefGerm[j]])
      
           # sigma_g[i,yearRefGerm[j],indexRefGerm[j]] ~ dnorm(0,1) T(0,)
     # sigma_g[i,yearRefGerm[j],indexRefGerm[j]] ~ dt(0,1/(5*5),1)I(0,)

      ## OPTION 1
      # tau_g[i,yearRefGerm[j],indexRefGerm[j]] ~ dgamma(1.7, 1.7)
      # sigma_g[i,yearRefGerm[j],indexRefGerm[j]] <-  1/sqrt(tau_g[i,yearRefGerm[j],indexRefGerm[j]])
      
      ## OPTION 2
      #  sigma_g[i,yearRefGerm[j],indexRefGerm[j]] ~ dunif(0,5)
       # tau_g[i,yearRefGerm[j],indexRefGerm[j]] <- 1/(sigma_g[i,yearRefGerm[j],indexRefGerm[j]]*sigma_g[i,yearRefGerm[j],indexRefGerm[j]])
       
      ## OPTION 3
      # Weakly informative prior (Rosenbaum et al. 2019)
      sigma_g[i,yearRefGerm[j],indexRefGerm[j]] ~ dnorm(0,1) T(0,)
      tau_g[i,yearRefGerm[j],indexRefGerm[j]] <- 1/(sigma_g[i,yearRefGerm[j],indexRefGerm[j]]*sigma_g[i,yearRefGerm[j],indexRefGerm[j]])
      
      sigma_pred[i,yearRefGerm[j],indexRefGerm[j]] ~ dnorm(0,1) T(0,)
      tau_pred[i,yearRefGerm[j],indexRefGerm[j]] <- 1/(sigma_pred[i,yearRefGerm[j],indexRefGerm[j]]*sigma_pred[i,yearRefGerm[j],indexRefGerm[j]])
      
      ## OPTION 4
     # sigma_g[i,yearRefGerm[j],indexRefGerm[j]] ~ dt(0, 10^-2, 1)T(0,)
    # tau_g[i,yearRefGerm[j],indexRefGerm[j]] <- 1/(sigma_g[i,yearRefGerm[j],indexRefGerm[j]]*sigma_g[i,yearRefGerm[j],indexRefGerm[j]])
      
      
      }
    
  }
  # likelihood (seed bags) --------------------------------------------------------------

  # for(i in 1:n1){
  #   
  #   # # get probability of germination for each value of index
  #   # alpha_g1[i] ~ dnorm(mu_g[siteSurvival[i],yearSurvival[i],gIndexSurvival[i]],tau_g[siteSurvival[i],yearSurvival[i],gIndexSurvival[i]])
  #   # logit(theta_g1[i]) <- alpha_g1[i]
  #   # 
  #   # alpha_g2[i] ~ dnorm(mu_g[siteSurvival[i],yearSurvival[i],gIndexSurvival[i]],tau_g[siteSurvival[i],yearSurvival[i],gIndexSurvival[i]])
  #   # logit(theta_g2[i]) <- alpha_g2[i]
  #   # 
  #   # alpha_g3[i] ~ dnorm(mu_g[siteSurvival[i],yearSurvival[i],gIndexSurvival[i]],tau_g[siteSurvival[i],yearSurvival[i],gIndexSurvival[i]])
  #   # logit(theta_g3[i]) <- alpha_g3[i]
  # 
  #   eta_surv[i] ~ dnorm(mu_s[siteSurvival[i],yearSurvival[i]],tau_s[siteSurvival[i],yearSurvival[i]])
  #   b[i] <- exp(-(eta_surv[i])/(a[siteSurvival[i]]))
  # 
  #   # mu[i] <- (((1-theta_g1[i])^(g_1[i]))*((1-theta_g2[i])^(g_2[i]))*((1-theta_g3[i])^(g_3[i])))*exp(-(months[i]/b[i])^a[siteSurvival[i]])
  #   # 
  #   # # beta is the 2nd term of beta dist.
  #   # beta_mod[i] <- beta[siteSurvival[i],betaIndex[i]]+.01
  #   # 
  #   # # use the beta term to write 1st term of beta dist.
  #   # alpha_s[i] <- (mu[i]*beta_mod[i])/(1-mu[i])
  #   # 
  #   # theta_s[i] ~ dbeta(alpha_s[i]+.01,beta_mod[i]) T(0.001,0.999)
  #   # 
  #   # # likelihood
  #   # y[i] ~ dbinom(theta_s[i], seedStart[i])
  #   #
  #   # # log-likelihood
  #   # logLik_y[i] <- logdensity.bin(y[i],theta_s[i],seedStart[i])
  # }

  for(i in 1:n2){

    alpha_g[i] ~ dnorm(mu_g[siteGermination[i],yearGermination[i],gIndex[i]],tau_g[siteGermination[i],yearGermination[i],gIndex[i]])
  #  alpha_g[i] ~ dnorm(mu0_g[siteGermination[i],gIndex[i]],tau_g[siteGermination[i],yearGermination[i],gIndex[i]])
    #alpha_g[i] <- mu_g[siteGermination[i],yearGermination[i],gIndex[i]]
    # logit
    logit(g[i]) <- alpha_g[i]

    # likelihood
    seedlingJan[i] ~ dbinom(g[i], totalJan[i])

    y_sim[i] ~ dbinom(g[i], totalJan[i])
    
    # prior predictive
    alpha_pred[i] ~ dnorm(mu_pred[siteGermination[i],yearGermination[i],gIndex[i]],tau_pred[siteGermination[i],yearGermination[i],gIndex[i]])
    logit(g_pred[i]) <- alpha_pred[i]
    
    y_pred[i] ~ dbinom(g_pred[i], totalJan[i])
    
    # log-likelihood
   # logLik_g[i] <- logdensity.bin(seedlingJan[i],g[i] ,totalJan[i])
  }


  # derived quantities --------------------------------------------------------------



}
