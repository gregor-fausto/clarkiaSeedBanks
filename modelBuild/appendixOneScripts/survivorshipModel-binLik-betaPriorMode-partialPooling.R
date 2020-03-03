model {
  
  for(j in 1:n_site){
    omega[j] ~ dbeta( 1 , 1 ) # broad uniform
    kappaMinusTwo[j]  ~ dgamma( .01 , .01 )  # mode=1 , sd=10 
    kappa[j]  <- kappaMinusTwo[j]  + 2
  }
  
  for(j in 1:n_site){   
    for ( i in 1:n_year ) {
      theta[j,i] ~ dbeta( omega[j]*(kappa[j]-2)+1 , (1-omega[j])*(kappa[j]-2)+1 ) 
    }
  }
  
  for ( i in 1:n ) {
    fruitingPlantNumber[i] ~ dbin( theta[site[i],year[i]], seedlingNumber[i] )
  }
}