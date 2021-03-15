model {

  # https://github.com/kionaogle/EcoApps_Random_Effects/blob/v1.0/Code%20A3.R
  # https://github.com/CCheCastaldo/CheCastaldo_etal_2019_EcolMongr/blob/master/ModelBuild/GlobalModel/i-009f9c75a9206a0ac/Implementation.R
  # https://github.com/colemonnahan/gradmcmc/tree/v1.0/models/wildflower
  # Priors --------------------------------------------------------------
  for(i in 1:n_siteSurvival){

    # Smits 2015 PNAS: assume the shape is the same but that
    # the scale varies by year
    # saying that a population has a characteristic pattern through time
    # but that the magnitude at which that pattern proceeds varies

    # shape parameter: population parameter
    # controls rate of change of age-specific mortality rate
    # change prior to be half-cauchy

    # seed survival ---------------------------------------------------------------
    a[i] ~ dgamma(2,2)

    # WEAKLY INFORMATIVE
    mu0_s[i] ~ dnorm(0,1)

    # WEAKLY INFORMATIVE
    sigma0_s[i] ~ dnorm(0,1) T(0,)
    tau0_s[i] <- 1/(sigma0_s[i]*sigma0_s[i])

    for(j in 1:n_yearSurvival){
      mu_s[i,j] ~ dnorm(mu0_s[i], tau0_s[i])

      # WEAKLY INFORMATIVE
      sigma_s[i,j] ~ dnorm(0,1) T(0,)
      tau_s[i,j] <- 1/(sigma_s[i,j]*sigma_s[i,j])
    }

    # seed germination ---------------------------------------------------------------
    # population*year parameter, hierarchical logit
    for(k in 1:3){
      mu0_g[i,k] ~ dnorm(0, 1)

      # Weakly informative prior (Rosenbaum et al. 2019)
      # slightly more informative prior because few groups https://statmodeling.stat.columbia.edu/2015/12/08/hierarchical-modeling-when-you-have-only-2-groups-i-still-think-its-a-good-idea-you-just-need-an-informative-prior-on-the-group-level-variation/
      sigma0_g[i,k] ~ dnorm(0,2) T(0,)
      tau0_g[i,k] <- 1/(sigma0_g[i,k]*sigma0_g[i,k])
    }

    # age 1 seeds observed in round 1-3
    for(j in c(1,4,6)){
      # i indexes site; j indexes year; k indexes germination variable
      mu_g[i,j] ~ dnorm(mu0_g[i,1], tau0_g[i,1])

      # Weakly informative prior (Rosenbaum et al. 2019)
      sigma_g[i,j] ~ dnorm(0,1) T(0,)
      tau_g[i,j] <- 1/(sigma_g[i,1]*sigma_g[i,1])
    }

    # age 2 seeds observed in round 1-2
    for(j in c(2,5)){
      # i indexes site; j indexes year; k indexes germination variable
      mu_g[i,j] ~ dnorm(mu0_g[i,2], tau0_g[i,2])
      
      # Weakly informative prior (Rosenbaum et al. 2019)
      sigma_g[i,j] ~ dnorm(0,1) T(0,)
      tau_g[i,j] <- 1/(sigma_g[i,2]*sigma_g[i,2])
    }
    
    # age 3 seeds observed in round 1
    for(j in 3){
      # i indexes site; j indexes year; k indexes germination variable
     # mu_g[i,yearRefGerm[j],indexRefGerm[j]] ~ dnorm(0, 1)
     mu_g[i,j] ~ dnorm(mu0_g[i,3], tau0_g[i,3])

      # Weakly informative prior (Rosenbaum et al. 2019)
      sigma_g[i,j] ~ dnorm(0,1) T(0,)
      tau_g[i,j] <- 1/(sigma_g[i,3]*sigma_g[i,3])
    }

    # unobserved s0 ---------------------------------------------------------------
    mu0_s0[i] ~ dnorm(0,1)

    # WEAKLY INFORMATIVE
    sigma0_s0[i] ~ dnorm(0,2) T(0,)
    tau0_s0[i] <- 1/(sigma0_s0[i]*sigma0_s0[i])

    for(j in 1:2){

      mu_s0[i,j] ~ dnorm(mu0_s0[i],tau0_s0[i])

      # WEAKLY INFORMATIVE
      sigma_s0[i,j] ~ dnorm(0,1) T(0,)
      tau_s0[i,j] <- 1/(sigma_s0[i,j]*sigma_s0[i,j])

    }

  }

  # Likelihood --------------------------------------------------------------

  # seed survival --------------------------------------------------------------
  for(i in 1:n1){

    # POPULATION*YEAR MARGINAL GERMINATION RATES FOR SURVIVAL
    # Round 1
    alpha_g1[i,1] ~ dnorm(mu_g[siteSurvival[i],1],tau_g[siteSurvival[i],1])
    logit(theta_g1[i,1]) <- alpha_g1[i,1]

    alpha_g1[i,2] ~ dnorm(mu_g[siteSurvival[i],2],tau_g[siteSurvival[i],2])
    logit(theta_g1[i,2]) <- alpha_g1[i,2]

    alpha_g1[i,3] ~ dnorm(mu_g[siteSurvival[i],3],tau_g[siteSurvival[i],3])
    logit(theta_g1[i,3]) <- alpha_g1[i,3]

    # Round 2
    alpha_g1[i,4] ~ dnorm(mu_g[siteSurvival[i],4],tau_g[siteSurvival[i],4])
    logit(theta_g1[i,4]) <- alpha_g1[i,4]

    alpha_g1[i,5] ~ dnorm(mu_g[siteSurvival[i],5],tau_g[siteSurvival[i],5])
    logit(theta_g1[i,5]) <- alpha_g1[i,5]

    # Round 3
    alpha_g1[i,6] ~ dnorm(mu_g[siteSurvival[i],6],tau_g[siteSurvival[i],6])
    logit(theta_g1[i,6]) <- alpha_g1[i,6]

    # COMPOSITE EVENT HISTORIES
    theta_c[i,1] = 1                   # jan - year/round 1 - age 1
    theta_c[i,2] = (1-theta_g1[i,1]) # oct - year/round 1 - age 1
    theta_c[i,3] = (1-theta_g1[i,1]) # jan - year/round 1 - age 2
    theta_c[i,4] = (1-theta_g1[i,1])*(1-theta_g1[i,2]) # oct - year/round 1 - age 2
    theta_c[i,5] = (1-theta_g1[i,1])*(1-theta_g1[i,2]) # jan - year/round 1 - age 3
    theta_c[i,6] = (1-theta_g1[i,1])*(1-theta_g1[i,2])*(1-theta_g1[i,3]) # oct - year/round 1 - age 3
    theta_c[i,7] =  1                   # jan - year/round 2 - age 1
    theta_c[i,8] = (1-theta_g1[i,4]) # oct - year/round 2 - age 1
    theta_c[i,9] = (1-theta_g1[i,4]) # jan - year/round 2 - age 2
    theta_c[i,10] = (1-theta_g1[i,4])*(1-theta_g1[i,5]) # oct - year/round 2 - age 2
    theta_c[i,11] =  1                   # jan - year/round 3 - age 1
    theta_c[i,12] = (1-theta_g1[i,6]) # oct - year/round 3 - age 1

    # DETERMINISTIC WEIBULL SURVIAL MODEL
    eta_surv[i] ~ dnorm(mu_s[siteSurvival[i],yearSurvival[i]],tau_s[siteSurvival[i],yearSurvival[i]])
    inv.b[i] <- exp(-(eta_surv[i])/(a[siteSurvival[i]]))
    mu[i] <- theta_c[siteSurvival[i],compIndex[i]]*exp(-(months[i]/inv.b[i])^a[siteSurvival[i]])

    ## LIKELIHOOD
    y[i] ~ dbinom(mu[i], seedStart[i])

    ## POSTERIOR PREDICTIVE
    y_sim[i] ~ dbinom(mu[i], seedStart[i])

    # Chi-squared
    chi2.yobs[i] <- pow((y[i]- mu[i]*seedStart[i]),2) / (mu[i]*seedStart[i]+.001)
    chi2.ysim[i] <- pow((y_sim[i]- mu[i]*seedStart[i]),2) / (mu[i]*seedStart[i]+.001)

    # LOG-LIKELIHOOD
    # logLik_y[i] <- logdensity.bin(y[i],mu[i],seedStart[i])

  }

  # seed germination --------------------------------------------------------------
  for(i in 1:n2){

    alpha_g[i] ~ dnorm(mu_g[siteGermination[i],germinationIndex[i]],tau_g[siteGermination[i],germinationIndex[i]])
    logit(g[i]) <- alpha_g[i]

    # LIKELIHOOD
    seedlingJan[i] ~ dbinom(g[i], totalJan[i])

    # POSTERIOR PREDICTIVE
    seedlingJan_sim[i] ~ dbinom(g[i], totalJan[i])

    # Chi-squared
    chi2.obs[i] <- pow((seedlingJan[i]- g[i]*totalJan[i]),2) / (g[i]*totalJan[i]+.001)
    chi2.sim[i] <- pow((seedlingJan_sim[i]- g[i]*totalJan[i]),2) / (g[i]*totalJan[i]+.001)

    # LOG-LIKELIHOOD
    # logLik_g[i] <- logdensity.bin(seedlingJan[i],g[i] ,totalJan[i])


  }

  # unobserved s0 -------------------------------------------------------------
  for(i in 1:n3){

    alpha_s0[i] ~ dnorm(mu_s0[sitePlot[i],yearPlot[i]],tau_s0[sitePlot[i],yearPlot[i]])
    logit(s0[i]) <- alpha_s0[i]

    # note the addition of +1 to the year index; data from year 1 and 2 for the aboveground
    # corresponds to years 2 and 3 for the seed bag data but indexing has to start at 1
    g_marg[i] ~ dnorm(mu_g[sitePlot[i],fecIndex[i]],tau_g[sitePlot[i],fecIndex[i]])
    logit(g_plot[i]) <- g_marg[i]

    # DETERMINISTIC WEIBULL SURVIAL MODEL
    eta_plot[i] ~ dnorm(mu_s[sitePlot[i],yearPlot[i]+1],tau_s[sitePlot[i],yearPlot[i]+1])
    inv.b_plot[i] <- exp(-(eta_plot[i])/(a[sitePlot[i]]+.01))
    mu_plot[i] <- g_plot[i]*exp(-((3/36)/inv.b_plot[i])^a[sitePlot[i]])*s0[i]

    # LIKELIHOOD
    plotSeedlings[i] ~ dbinom(mu_plot[i], fec[i])

    # POSTERIOR PREDICTIVE
    plotSeedlings_sim[i] ~ dbinom(mu_plot[i], fec[i])

    # Chi-squared
    chi2.plot.obs[i] <- pow((plotSeedlings[i]- mu_plot[i]*totalJan[i]),2) / (mu_plot[i]*fec[i]+.001)         # obs.
    chi2.plot.sim[i] <- pow((plotSeedlings_sim[i]- mu_plot[i]*totalJan[i]),2) / (mu_plot[i]*fec[i]+.001)         # obs.

    # LOG-LIKELIHOOD
    # logLik_yplot[i] <- logdensity.bin(plotSeedlings[i],mu_plot[i] , fec[i])

  }


}
