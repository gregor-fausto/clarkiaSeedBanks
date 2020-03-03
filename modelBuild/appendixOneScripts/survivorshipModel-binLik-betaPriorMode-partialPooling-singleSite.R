model {
  
  for(i in 1:n_year){
    omega[i] ~ dbeta( 1 , 1 ) # broad uniform
    kappaMinusTwo[i]  ~ dgamma( .01 , .01 )  # mode=1 , sd=10 
    kappa[i]  <- kappaMinusTwo[i]  + 2
  }
  
  for(i in 1:n_year){   
    for ( j in 1:n_plot ) {
      theta[i,j] ~ dbeta( omega[i]*(kappa[i]-2)+1 , (1-omega[i])*(kappa[i]-2)+1 ) 
    }
  }
  
  for ( i in 1:n ) {
    fruitingPlantNumber[i] ~ dbin( theta[year[i],plot[i]], seedlingNumber[i] )
  }
  
  for(i in 1:n_year){   
      p[i] ~ dbeta( omega[i]*(kappa[i]-2)+1 , (1-omega[i])*(kappa[i]-2)+1 ) 
    }

}