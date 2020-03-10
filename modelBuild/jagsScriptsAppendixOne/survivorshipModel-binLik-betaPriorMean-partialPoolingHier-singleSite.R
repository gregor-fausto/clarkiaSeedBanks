model {

    phi0 ~ dunif(0,1)
  # check dpar
    kappa0 ~ dpar(1.5,1) 
  
  for(k in 1:n_site){
    for(i in 1:n_year){
      phi[k,i] ~ dbeta(kappa0*phi0, kappa0*(1-phi0)) #T(0,0.95)
      kappa[k,i] ~ dpar(1.5,1) 
  }
  }

  for(k in 1:n_site)  {  
    for(i in 1:n_year) {   
      for ( j in 1:n_plot ) {
        theta[k,i,j] ~ dbeta(kappa[k,i]*phi[k,i], kappa[k,i]*(1-phi[k,i]))
      }
    }
  }
  
  for ( i in 1:n ) {
    fruitingPlantNumber[i] ~ dbin( theta[site[i],year[i],plot[i]], seedlingNumber[i] )
    fruitingPlantNumberSim[i] ~ dbin( theta[site[i],year[i],plot[i]], seedlingNumber[i])
    
    }
  
  # marginalize
  for(k in 1:n_site){
    p_site[k] ~ dbeta(kappa0[k]*phi0[k], kappa0[k]*(1-phi0[k]))

  for(i in 1:n_year){
      p[k,i] ~ dbeta(kappa[k,i]*phi[k,i], kappa[k,i]*(1-phi[k,i]))
  }
  }
  
    # create output for ppc checks
}
