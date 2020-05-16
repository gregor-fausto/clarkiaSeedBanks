
model { 
 
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
   
  # for(j in 1:n_site){
  #   
  #   alpha[j] ~ dnorm(0, .001)
  #   
  #   for(k in 1:n_year){
  #     
  #     gamma[j,k] ~ dnorm(0, 0.001)
  #     r[j,k] ~ dgamma(.001,.001)
  # 
  #   }
  # }
  
  # likelihoods
  for (i in 1:n){
    lambda[i] = exp(gamma[site[i],year[i]])
    countFruitsPerPlant[i] ~ dnegbin(r[site[i],year[i]]/(r[site[i],year[i]]+lambda[i]),r[site[i],year[i]])

  # simulated data for posterior predictive checks
  # y1.sim[i] ~ dnegbin(rF[siteFec[i],yearFec[i]]/(rF[siteFec[i],yearFec[i]]+lambdaF[i]),rF[siteFec[i],yearFec[i]])
  }

   # derived quantity
   for(i in 1:n_site){

         p0.i[i] ~ dnorm(mu0[i],tau0[i])
       p0[i] <- exp(p0.i[i])
        
       }
  
  
}
