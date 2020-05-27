
model { 
 
  # undamaged
  for(i in 1:n_site){
    
    # theta 1
    mu0[i] ~  dnorm(0, 0.001)
    sigma0[i] ~ dunif(0,1.5)
    tau0[i] <- 1/(sigma0[i]*sigma0[i])
    
    for(j in 1:n_year){
      
      # theta 1
      gamma[i,j] ~ dnorm(mu0[i], tau0[i])
      r[i,j] ~ dgamma(.001,.001)
      
    }
  }
  
  # damaged
  for(i in 1:n_site){
    
    # theta 1
    mu0_dam[i] ~  dnorm(0, 0.001)
    sigma0_dam[i] ~ dunif(0,1.5)
    tau0_dam[i] <- 1/(sigma0_dam[i]*sigma0_dam[i])
    
    for(j in 1:n_year){
      
      # theta 1
      gamma_dam[i,j] ~ dnorm(mu0_dam[i], tau0_dam[i])
      r_dam[i,j] ~ dgamma(.001,.001)
      
    }
  }
  
  # likelihoods
  for (i in 1:n){
    lambda[i] = exp(gamma[site[i],year[i]])
    sdno[i] ~ dnegbin(r[site[i],year[i]]/(r[site[i],year[i]]+lambda[i]),r[site[i],year[i]])
  }
  
  for (i in 1:n){
    lambda_dam[i] = exp(gamma_dam[site[i],year[i]])
    sdno_dam[i] ~ dnegbin(r_dam[site[i],year[i]]/(r_dam[site[i],year[i]]+lambda_dam[i]),r_dam[site[i],year[i]])
  }
  
  # for other use
  exp(gamma_dam[site[i],year[i]])/exp(gamma[site[i],year[i]])
  
  ratio[j,k] exp(gamma_dam[site[i],year[i]] - gamma[site[i],year[i]])
  
  for (i in 1:n){
    lambda[i] = exp(gamma[site[i],year[i]])
    
    for (year[i] in 2006:2012) {
      countFruitsPerPlant[i] ~ dnegbin(r[site[i],year[i]]/(r[site[i],year[i]]+lambda[i]),r[site[i],year[i]])
    }
    for (year[i] in 2013:2018) {
      ratio[site[i],year[i]] = exp(gamma_dam[site[i],year[i]] - gamma[site[i],year[i]])
      countTFE[i] = n_undamaged[i] + round( ratio[j,k] * n_damaged[i] )
      
      countTFE[i] ~ dnegbin(r[site[i],year[i]]/(r[site[i],year[i]]+lambda[i]),r[site[i],year[i]])
      
    }
    
    
    lambda[i] = exp(gamma[site[i],year[i]])
    
    ratio[site[i],year[i]] exp(gamma_dam[site[i],year[i]] - gamma[site[i],year[i]])
    
    n_undamaged + round( ratio[j,k] * n_damaged ) ~ dnegbin(r[site[i],year[i]]/(r[site[i],year[i]]+lambda[i]),r[site[i],year[i]])
    countFruitsPerPlant[i] ~ dnegbin(r[site[i],year[i]]/(r[site[i],year[i]]+lambda[i]),r[site[i],year[i]])
    
  
}
