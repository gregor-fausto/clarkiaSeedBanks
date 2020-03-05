model {
  
  for(i in 1:n_year){
    mu[i] ~ dnorm(0, 0.001)
    sigma[i] ~ dunif(0,100)
    tau[i] <- 1/(sigma[i]*sigma[i])
  }
  
  for(i in 1:n_year){   
    for ( j in 1:n_plot ) {
        alpha[i,j] ~ dnorm(mu[i], tau[i])
    }
  }
  
  ## Likelihood
  
  for ( i in 1:n ) {
    theta[i] <- ilogit(alpha[year[i],plot[i]] )
    
    fruitingPlantNumber[i] ~ dbin( theta[i], seedlingNumber[i] )
    fruitingPlantNumberSim[i] ~ dbin( theta[i], seedlingNumber[i])
    
    }
  
  for(i in 1:n_year){   
      alpha.year[i] ~ dnorm(mu[i], tau[i])
      p[i] <- ilogit(alpha.year[i])
  }
}
