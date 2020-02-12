
model { 
    
  for(i in 1:N){
    
    # bag prior
    pg[i] ~ dunif(0, 1)
    pv[i] ~ dunif(0, 1)

    # g germination
      yg[i] ~ dbinom( pg[i], ng[i] )
      
    # v viability
      yv[i] ~ dbinom( pv[i] , nv[i] )

    yg.sim[i] ~ dbinom( pg[i], ng[i] )
    yv.sim[i] ~ dbinom( pv[i], nv[i] )
    
  }
  
    for(i in 1:N){
        viability[i] = pg[i] + pv[i]*(1-pg[i])
        }
  
}

