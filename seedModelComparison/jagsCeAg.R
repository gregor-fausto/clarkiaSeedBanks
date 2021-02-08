model {

  # seed bags ---------------------------------------------------------------

  a ~ dgamma(0.001,0.001)
  b ~ dgamma(0.001,0.001)

  for(i in 1:6){
    beta[i] ~ dgamma(1,0.001)
  }

  for(i in 1:3){
    g[i] ~ dbeta(1,1)
  }

  # likelihood (seed bags) --------------------------------------------------------------

  for(i in 1:n1){

    mu[i] <- (((1-g[germIndex[i]])^(g_1[i]))*((1-g[germIndex[i]])^(g_2[i]))*((1-g[germIndex[i]])^(g_3[i])))*(1/((1+b*months[i])^a))
    beta_mod[i] <- beta[betaIndex[i]]+.01
    alpha[i] <- (mu[i]*beta_mod[i])/(1-mu[i])
    theta[i] ~ dbeta(alpha[i]+.01,beta_mod[i])

    # likelihood
    y[i] ~ dbinom(theta[i], seedStart[i])

    # log-likelihood
    logLik_y[i] <- logdensity.bin(y[i],theta[i],seedStart[i])
  }

  for(i in 1:n2){
    # likelihood
    seedlingJan[i] ~ dbinom(g[gIndex[i]], totalJan[i])

    # log-likelihood
    logLik_g[i] <- logdensity.bin(seedlingJan[i],g[gIndex[i]] ,totalJan[i])
  }


  # derived quantities --------------------------------------------------------------



}
