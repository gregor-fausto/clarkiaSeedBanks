model { 
  
  # Priors
  for(j in 1:n_year){
      p[j] ~ dbeta(1,1)
  }
  
  # Likelihood
  for(i in 1:n){
    fruitingPlantNumber[i] ~ dbin(p[year[i]], seedlingNumber[i])
    fruitingPlantNumberSim[i] ~ dbin(p[year[i]], seedlingNumber[i])
  }
  
}