model { 
 
  # hyperpriors
  
   theta.g ~ dunif(0,1)
  # kappa ~ dpar(alpha=shape,c=scale)
    kappa.g ~ dpar(1.5,1)

    theta.v ~ dunif(0,1)
    kappa.v ~ dpar(1.5,1)
  
  # priors
  
  for(j in 1:nbags){
      pg[j] ~ dbeta(kappa.g*theta.g, kappa.g*(1-theta.g))
      pv[j] ~ dbeta(kappa.v*theta.v, kappa.v*(1-theta.v))

  }
  
  for(i in 1:N){
      # g germination
      yg[i] ~ dbinom( pg[bag[i]] , ng[i] )
      
      # v viability
      yv[i] ~ dbinom( pv[bag[i]] , nv[i] )

      yg.sim[i] ~ dbinom( pg[bag[i]], ng[i] )
      yv.sim[i] ~ dbinom( pv[bag[i]], nv[i] )
    
  }
    
      # derived quantity: viability
    for ( j in 1:nbags){
        viability[j] = pg[j] + pv[j]*(1-pg[j])
        }
  
 
  
}
