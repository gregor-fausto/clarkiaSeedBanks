model { 
  
  # hyperpriors

  mu.g ~ dnorm(0,.0001)    
  sigma.g ~ dunif(0,100)
  tau.g <- 1/(sigma.g * sigma.g)
    
  mu.v ~ dnorm(0,.0001)      
  sigma.v ~ dunif(0,100)
  tau.v <- 1/(sigma.v * sigma.v)

  # priors
  
  for(j in 1:nbags){
      alpha.g[j] ~ dnorm(mu.g, tau.g)
      alpha.v[j] ~ dnorm(mu.v, tau.v)
  }
  
  #l ikelihood
  
    for(i in 1:N){
        
    # g germination
    pg[i] <- ilogit(alpha.g[bag[i]])
    yg[i] ~ dbinom( pg[i] , ng[i] )
    
    # v viability
    pv[i] <- ilogit(alpha.v[bag[i]])
    yv[i] ~ dbinom( pv[i] , nv[i] )

        yg.sim[i] ~ dbinom( pg[i], ng[i] )
        yv.sim[i] ~ dbinom( pv[i], nv[i] )
    
  }
  
  # derived quantity: viability
    for ( j in 1:nbags){
        t.g[j] = ilogit(alpha.g[j])
        t.v[j] = ilogit(alpha.v[j])
        viability[j] = t.g[j] + t.v[j]*(1-t.g[j])
        }
  
}
