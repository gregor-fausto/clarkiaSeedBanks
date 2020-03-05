model {
  
  mu0 ~  dnorm(0, 0.001)
  sigma0 ~ dunif(0,100)
  tau0 <- 1/(sigma0*sigma0)
  
  for(k in 1:n_site){
    for(i in 1:n_year){
      mu[k,i] ~ dnorm(mu0, tau0)
      sigma[k,i] ~ dunif(0,100)
      tau[k,i] <- 1/(sigma[k,i]*sigma[k,i])
      
    }
  }
  
  for(k in 1:n_site) {
  for(i in 1:n_year){   
    for ( j in 1:n_plot ) {
        alpha[k,i,j] ~ dnorm(mu[k,i],tau[k,i])
      }
  }
  }
  
  ## Likelihood
  
  for ( i in 1:n ) {
    theta[i] <- ilogit(alpha[site[i],year[i],plot[i]] )
    
    fruitingPlantNumber[i] ~ dbin( theta[i], seedlingNumber[i] )
    fruitingPlantNumberSim[i] ~ dbin( theta[i], seedlingNumber[i])
    
    }
  
  # for(i in 1:n_year){   
  #     alpha.year[i] ~ dnorm(mu[i], tau[i])
  #     p[i] <- ilogit(alpha.year[i])
  # }
}
