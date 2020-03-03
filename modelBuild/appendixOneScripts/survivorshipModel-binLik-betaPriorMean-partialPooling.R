model {
  
  for(j in 1:n_site){
      phi[j] ~ dunif(0,1)
      kappa[j] ~ dpar(1.5,1) 
  }
  
  for(j in 1:n_site){   
    for ( i in 1:n_year ) {
        theta[j,i] ~ dbeta(kappa[j]*phi[j], kappa[j]*(1-phi[j]))
    }
  }
  
  for ( i in 1:n ) {
    fruitingPlantNumber[i] ~ dbin( theta[site[i],year[i]], seedlingNumber[i] )
  }
}
