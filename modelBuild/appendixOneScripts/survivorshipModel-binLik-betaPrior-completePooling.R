model { 
  
  # Priors
  for(j in 1:n_site){
    for(k in 1:n_year){
      p[j,k] ~ dbeta(1,1)
    }
  }
  
  # Likelihood
  for(i in 1:n){
    fruitingPlantNumber[i] ~ dbin(p[site[i],year[i]], seedlingNumber[i])
  }
  
}