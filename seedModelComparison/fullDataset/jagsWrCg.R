model {

  # seed bags ---------------------------------------------------------------

  for(i in 1:n_siteSurvival){
    for(j in 1:6){
      beta[i,j] ~ dgamma(1,0.001)
    }
  }

  for(i in 1:n_siteSurvival){
    a[i] ~ dgamma(0.001,0.001)
    b[i] ~ dgamma(0.001,0.001)
    g[i] ~ dbeta(1,1)
  }
  # likelihood (seed bags) --------------------------------------------------------------

  for(i in 1:n1){

    mu[i] <- (((1-g[siteSurvival[i]])^(g_1[i]))*((1-g[siteSurvival[i]])^(g_2[i]))*((1-g[siteSurvival[i]])^(g_3[i])))*exp(-(months[i]/b[siteSurvival[i]])^a[siteSurvival[i]])
    beta_mod[i] <- beta[siteSurvival[i],betaIndex[i]]+.01
    alpha[i] <- (mu[i]*beta_mod[i])/(1-mu[i])
    theta[i] ~ dbeta(alpha[i]+.01,beta_mod[i]) T(0.001,0.999)

    # likelihood
    y[i] ~ dbinom(theta[i], seedStart[i])

    # log-likelihood
    logLik_y[i] <- logdensity.bin(y[i],theta[i],seedStart[i])
  }

  for(i in 1:n2){

    # likelihood
    seedlingJan[i] ~ dbinom(g[siteGermination[i]], totalJan[i])

    # log-likelihood
    logLik_g[i] <- logdensity.bin(seedlingJan[i],g[siteGermination[i]] ,totalJan[i])

  }

  # derived quantities --------------------------------------------------------------



}
