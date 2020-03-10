model { 
  
  # Likelihood
  for(i in 1:n){
    fruitingPlantNumber[i] ~ dbin(p[site[i],year[i]], seedlingNumber[i])
  }
  
}