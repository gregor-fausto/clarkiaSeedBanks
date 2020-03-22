model { 

  # priors

    for(i in 1:n_site){
    mu0[i] ~ dnorm(0, .001)
    sigma0[i] ~ dunif(0,1.5)
    tau0[i] <- 1/(sigma0[i]*sigma0[i])
   
    for(k in 1:n_year){
      #  alpha_mu.std[i,k] ~ dnorm(0,1)
       # alpha_mu[i,k] <- mu0[i] + sigma0[i]*alpha_mu.std[i,k]
        mu[i,k] ~ dnorm(mu0[i], tau0[i])
        sigma[i,k] ~ dunif(0,1.5)
    }
    }
  
  # Likelihood
  for(i in 1:n){
    alpha.std[i] ~ dnorm(0,1)
    alpha[i]<- mu[site[i],year[i]] + sigma[site[i],year[i]]*alpha.std[i]
    logit(theta[i]) <- alpha[i]
    
    totalJan[i] ~ dbinom(theta[i], seedStart[i])
  }

                                        # derived quantities
    for(i in 1:n_site){
   p_pop.i[i] ~ dnorm(mu0[i], tau0[i])
   logit(p_pop[i]) <- p_pop.i[i]
   }
  
}
