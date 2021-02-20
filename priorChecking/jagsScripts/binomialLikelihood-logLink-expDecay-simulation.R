model {

  # https://github.com/kionaogle/EcoApps_Random_Effects/blob/v1.0/Code%20A3.R
  #https://github.com/CCheCastaldo/CheCastaldo_etal_2019_EcolMongr/blob/master/ModelBuild/GlobalModel/i-009f9c75a9206a0ac/Implementation.R
  # https://github.com/colemonnahan/gradmcmc/tree/v1.0/models/wildflower
  # seed bags ---------------------------------------------------------------

  for(i in 1:n_siteSurvival){

    # Smits 2015 PNAS: assume the shape is the same but that
    # the scale varies by year
    # saying that a population has a characteristic pattern through time
    # but that the magnitude at which that pattern proceeds varies

    # shape parameter: population parameter
    # controls rate of change of age-specific mortality rate
    # change prior to be half-cauchy

    # a[i] ~ dscaled.gamma(2.5,1)
    # a[i] ~ dlnorm(0,10) T(0,)
    # a[i]~ dt(0,1/2.5^2,1) I(0,)

    # scale parameter: population*year parameter
    # b[i] ~ dnorm(0,.001) T(0,)
    # b[i] ~ dscaled.gamma(2.5,1)
    # b[i] ~ dgamma(0.001,0.001)

    ## half-cauchy seems to have a good property
    #  b[i] ~ dt(0,1/2.5^2,1) I(0,)
    
    #mu~n(0,1);sigma(0,10)t(0,) covers [0,150] for lambda
    #but heavy tailed on lower end
    
    
    ## Weakly informative prior
  #  mu0_s[i] ~ dt(0,10,4) I(0,)
    mu0_s[i] ~ dnorm(0,1)
    
    # WEAKLY INFORMATIVE
   sigma0_s[i] ~ dnorm(0,1) T(0,)
   tau0_s[i] <- 1/(sigma0_s[i]*sigma0_s[i])

    # HALF CAUCHY
     # sigma0_s[i] ~ dt(0,1/5,7) I(0,)
     # tau0_s[i] <- 1/(sigma0_s[i]*sigma0_s[i])

    for(j in 1:n_yearSurvival){
      mu_s[i,j] ~ dnorm(mu0_s[i], tau0_s[i])

      # WEAKLY INFORMATIVE
     sigma_s[i,j] ~ dnorm(0,1) T(0,)
     tau_s[i,j] <- 1/(sigma_s[i,j]*sigma_s[i,j])

      # HALF CAUCHY
       # sigma_s[i,j] ~ dt(0,1/5,7) I(0,)
       # tau_s[i,j] <- 1/(sigma_s[i,j]*sigma_s[i,j])
    }

   # 2nd term of the beta distribution
   #  for(j in 1:6){
   # #beta[i,j] ~ dt(0,1/2.5^2,1) I(0,)
   # beta[i,j] ~ dnorm(0,.001) T(0,)
   # # beta[i,j] ~ dgamma(.001,0.001)
   #   }

    # population*year parameter, hierarchical logit
    for(k in 1:3){
      mu0_g[i,k] ~ dnorm(0, 1)

      # sigma0_g[i,k] ~ dt(0,1/.5^2,4) I(0,)
      # sigma0_g[i,k] ~ dunif(0,5)

      # Weakly informative prior (Rosenbaum et al. 2019)
      sigma0_g[i,k] ~ dnorm(0,1) T(0,)
      tau0_g[i,k] <- 1/(sigma0_g[i,k]*sigma0_g[i,k])

      mu0_pred[i,k] ~ dnorm(0,1)
      sigma0_pred[i,k] ~ dnorm(0,1) T(0,)
      tau0_pred[i,k] <- 1/(sigma0_pred[i,k]*sigma0_pred[i,k])

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
      
      # # POPULATION*YEAR MARGINAL GERMINATION RATES FOR SURVIVAL
      # alpha_g1[i,yearRefGerm[j],indexRefGerm[j]] ~ dnorm(mu_g[i,yearRefGerm[j],indexRefGerm[j]],tau_g[i,yearRefGerm[j],indexRefGerm[j]])
      # logit(theta_g1[i,yearRefGerm[j],indexRefGerm[j]]) <- alpha_g1[i,yearRefGerm[j],indexRefGerm[j]]

    }

    # theta_c[i,1] = 1                   # jan - year/round 1 - age 1
    # theta_c[i,2] = (1-theta_g1[i,1,1]) # oct - year/round 1 - age 1
    # theta_c[i,3] = (1-theta_g1[i,1,1]) # jan - year/round 1 - age 2
    # theta_c[i,4] = (1-theta_g1[i,1,1])*(1-theta_g1[i,1,2]) # oct - year/round 1 - age 2
    # theta_c[i,5] = (1-theta_g1[i,1,1])*(1-theta_g1[i,1,2]) # jan - year/round 1 - age 3
    # theta_c[i,6] = (1-theta_g1[i,1,1])*(1-theta_g1[i,1,2])*(1-theta_g1[i,1,3]) # oct - year/round 1 - age 3
    # theta_c[i,7] =  1                   # jan - year/round 2 - age 1
    # theta_c[i,8] = (1-theta_g1[i,2,1]) # oct - year/round 2 - age 1
    # theta_c[i,9] = (1-theta_g1[i,2,1]) # jan - year/round 2 - age 2
    # theta_c[i,10] = (1-theta_g1[i,2,1])*(1-theta_g1[i,2,2]) # oct - year/round 2 - age 2
    # theta_c[i,11] =  1                   # jan - year/round 3 - age 1
    # theta_c[i,12] = (1-theta_g1[i,3,1]) # oct - year/round 3 - age 1
  }

  # likelihood (seed bags) --------------------------------------------------------------

  for(i in 1:n1){
    
    # POPULATION*YEAR MARGINAL GERMINATION RATES FOR SURVIVAL
    alpha_g1[i,1,1] ~ dnorm(mu_g[siteSurvival[i],1,1],tau_g[siteSurvival[i],1,1])
    logit(theta_g1[i,1,1]) <- alpha_g1[i,1,1]
    
    alpha_g1[i,1,2] ~ dnorm(mu_g[siteSurvival[i],1,2],tau_g[siteSurvival[i],1,2])
    logit(theta_g1[i,1,2]) <- alpha_g1[i,1,2]
    
    alpha_g1[i,1,3] ~ dnorm(mu_g[siteSurvival[i],1,3],tau_g[siteSurvival[i],1,3])
    logit(theta_g1[i,1,3]) <- alpha_g1[i,1,3]
    
    
    # alpha_g1[i,2,1] ~ dnorm(mu_g[siteSurvival[i],2,1],tau_g[siteSurvival[i],2,1])
    # logit(theta_g1[i,2,1]) <- alpha_g1[i,2,1]
    # 
    # alpha_g1[i,2,2] ~ dnorm(mu_g[siteSurvival[i],2,2],tau_g[siteSurvival[i],2,2])
    # logit(theta_g1[i,2,2]) <- alpha_g1[i,2,2]
    # 
    # 
    # alpha_g1[i,3,1] ~ dnorm(mu_g[siteSurvival[i],3,1],tau_g[siteSurvival[i],3,1])
    # logit(theta_g1[i,3,1]) <- alpha_g1[i,3,1]
    
    theta_c[i,1] = 1                   # jan - year/round 1 - age 1
    theta_c[i,2] = (1-theta_g1[i,1,1]) # oct - year/round 1 - age 1
    theta_c[i,3] = (1-theta_g1[i,1,1]) # jan - year/round 1 - age 2
    theta_c[i,4] = (1-theta_g1[i,1,1])*(1-theta_g1[i,1,2]) # oct - year/round 1 - age 2
    theta_c[i,5] = (1-theta_g1[i,1,1])*(1-theta_g1[i,1,2]) # jan - year/round 1 - age 3
    theta_c[i,6] = (1-theta_g1[i,1,1])*(1-theta_g1[i,1,2])*(1-theta_g1[i,1,3]) # oct - year/round 1 - age 3
    # theta_c[i,7] =  1                   # jan - year/round 2 - age 1
    # theta_c[i,8] = (1-theta_g1[i,2,1]) # oct - year/round 2 - age 1
    # theta_c[i,9] = (1-theta_g1[i,2,1]) # jan - year/round 2 - age 2
    # theta_c[i,10] = (1-theta_g1[i,2,1])*(1-theta_g1[i,2,2]) # oct - year/round 2 - age 2
    # theta_c[i,11] =  1                   # jan - year/round 3 - age 1
    # theta_c[i,12] = (1-theta_g1[i,3,1]) # oct - year/round 3 - age 1
    

   # # EXPONENTIAL MODEL
    eta_surv[i] ~ dnorm(mu_s[siteSurvival[i],yearSurvival[i]],tau_s[siteSurvival[i],yearSurvival[i]])
   ## lambda parameter
    b[i] <- exp(eta_surv[i])
  # # model without germination pieces
    mu_survival[i] <- exp(-b[i]*months[i])
    mu[i] <- theta_c[i,compIndex[i]]*exp(-b[i]*months[i])

    # WEIBULL MODEL
    # eta_surv[i] ~ dnorm(mu_s[siteSurvival[i],yearSurvival[i]],tau_s[siteSurvival[i],yearSurvival[i]])
    # inv.b[i] <- exp(-(eta_surv[i])/(a[siteSurvival[i]]))
    # # model without germination pieces
    # mu_survival[i] <- exp(-(months[i]/inv.b[i])^a[siteSurvival[i]])
    # mu[i] <- theta_c[siteSurvival[i],compIndex[i]]*exp(-(months[i]/inv.b[i])^a[siteSurvival[i]])

    # EXPONENTIAL AND WEIBULL CODE
    # mu[i] <- (((1-theta_g1[i])^(g_1[i]))*((1-theta_g2[i])^(g_2[i]))*((1-theta_g3[i])^(g_3[i])))*exp(-b[i]*months[i])
    # mu[i] <- (((1-theta_g1[i])^(g_1[i]))*((1-theta_g2[i])^(g_2[i]))*((1-theta_g3[i])^(g_3[i])))*exp(-(months[i]/b[i])^a[siteSurvival[i]])

    # beta is the 2nd term of beta dist.
    # beta_mod[i] <- beta[siteSurvival[i],betaIndex[i]] +.01
    # # # use the beta term to write 1st term of beta dist.
    # alpha_s[i] <- (mu[i]*beta_mod[i])/(1-mu[i])
    # theta_s[i] ~ dbeta(alpha_s[i]+.01,beta_mod[i]) #T(0.001,0.999)

  #   ## LIKELIHOOD
  #   y[i] ~ dbinom(theta_s[i], seedStart[i])
     y[i] ~ dbinom(mu[i], seedStart[i])
    
  #   ## POSTERIOR PREDICTIVE
  #   # y_surv[i] ~ dbinom(theta_s[i], seedStart[i])
  # 
  #   ## LOG-LIKELIHOOD
      logLik_y[i] <- logdensity.bin(y[i],mu[i],seedStart[i])
  #   
  #   # PRIOR PREDICTIVE
  #   y_prior.pred[i] ~ dbinom(mu[i], seedStart[i])
  }
  

  # for(i in 1:n_pred){
  #   eta[i] ~ dnorm(mu_s[siteSurvival_pred[i],yearSurvival_pred[i]],tau_s[siteSurvival_pred[i],yearSurvival_pred[i]])
  #   b.pred[i] <- exp(eta_surv[i])
  #   # model without germination pieces
  #   mu_survival.pred[i] <- exp(-b.pred[i]*months_pred[i])
  #   mu.pred[i] <- theta_c[siteSurvival_pred[i],compIndex_pred[i]]*exp(-b.pred[i]*months_pred[i])
  #   y_prior.pred[i] ~ dbinom(mu.pred[i], seedStart_pred[i])
  # }

  for(i in 1:n2){

    alpha_g[i] ~ dnorm(mu_g[siteGermination[i],yearGermination[i],gIndex[i]],tau_g[siteGermination[i],yearGermination[i],gIndex[i]])
    logit(g[i]) <- alpha_g[i]

    ## try a couple of alternatives
    ## alpha_g[i] ~ dnorm(mu0_g[siteGermination[i],gIndex[i]],tau_g[siteGermination[i],yearGermination[i],gIndex[i]])
    ## alpha_g[i] <- mu_g[siteGermination[i],yearGermination[i],gIndex[i]]

    # # LIKELIHOOD
     seedlingJan[i] ~ dbinom(g[i], totalJan[i])

    # POSTERIOR PREDICTIVE
    # y_sim[i] ~ dbinom(g[i], totalJan[i])

    ## LOG-LIKELIHOOD
    # logLik_g[i] <- logdensity.bin(seedlingJan[i],g[i] ,totalJan[i])

    # PRIOR PREDICTIVE
    alpha_pred[i] ~ dnorm(mu_pred[siteGermination[i],yearGermination[i],gIndex[i]],tau_pred[siteGermination[i],yearGermination[i],gIndex[i]])
    logit(g_pred[i]) <- alpha_pred[i]
    #seedlingJan2[i] ~ dbinom(g_pred[i], totalJan[i])
    
    y_pred[i] ~ dbinom(g_pred[i], totalJan[i])


  }


  # derived quantities --------------------------------------------------------------



}
