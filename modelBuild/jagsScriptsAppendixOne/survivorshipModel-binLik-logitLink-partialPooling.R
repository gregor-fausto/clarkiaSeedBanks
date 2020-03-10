model { 
  ## Priors
  ##hyperprior for intercept alpha
  
  for(i in 1:n_site){
    sigma.site[i] ~ dunif(0,100)
    tau.site[i] <- 1/(sigma.site[i]*sigma.site[i])
    mu.alpha[i] ~ dnorm(0, 0.0001)
  }
  
  for(i in 1:n_site){
    for(j in 1:n_year){
      #site intercepts
      alpha.i[i,j] ~ dnorm(mu.alpha[site[i]], tau.site[site[i]])
    }
  }
  
  ## Likelihood
  for(i in 1:n){
    theta[i] <- ilogit(alpha.i[site[i],year[i]] )
    fruitingPlantNumber[i] ~ dbin(theta[i], seedlingNumber[i])
  }
  
}