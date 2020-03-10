model {
  
  omega ~ dbeta( 1 , 1 ) # broad uniform
  kappaMinusTwo ~ dgamma( .01 , .01 )  # mode=1 , sd=10 
  kappa <- kappaMinusTwo + 2
  
  for ( i in 1:n_year ) {
    theta[i] ~ dbeta( omega*(kappa-2)+1 , (1-omega)*(kappa-2)+1 ) 
  }
  
  for ( i in 1:n ) {
    fruitingPlantNumber[i] ~ dbin( theta[year[i]], seedlingNumber[i] )
  }
}