model {
  
  for(i in 1:n_year){
      phi[i] ~ dunif(0,1)
      kappa[i] ~ dpar(1.5,1) 
  }
  
  for(i in 1:n_year){   
    for ( j in 1:n_plot ) {
        theta[i,j] ~ dbeta(kappa[i]*phi[i], kappa[i]*(1-phi[i]))
    }
  }
  
  for ( i in 1:n ) {
    fruitingPlantNumber[i] ~ dbin( theta[year[i],plot[i]], seedlingNumber[i] )
    fruitingPlantNumberSim[i] ~ dbin( theta[year[i],plot[i]], seedlingNumber[i])
    
    }
  
  for(i in 1:n_year){   
      p[i] ~ dbeta(kappa[i]*phi[i], kappa[i]*(1-phi[i]))
  }
}
