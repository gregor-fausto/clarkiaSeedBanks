model { 
  
  for(j in 1:nbags){
      pg[j] ~ dbeta(1, 1)
      pv[j] ~ dbeta(1, 1)
  }
  
  for(i in 1:N){

      # g germination
      yg[i] ~ dbinom( pg[bag[i]] , ng[i] )
      
    # v viability
    yv[i] ~ dbinom( pv[bag[i]] , nv[i] )

      yg.sim[i] ~ dbinom( pg[bag[i]], ng[i] )
      yv.sim[i] ~ dbinom( pv[bag[i]], nv[i] )
    
  }
  
                                        # derived quantity block: viability
    for(j in 1:nbags){
    viability[j] = pg[j]+pv[j]*(1-pg[j])
  }
}
