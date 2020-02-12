
model { 
    
  for(i in 1:N){
    
    # bag prior
    p[i] ~ dunif(0, 1)
    
    # v viability
    yv[i] ~ dbinom( p[i] , nv[i] )
    
    yv.sim[i] ~ dbinom( p[i], nv[i] )
    
    # # code for deviance from Lunn 2013
    prop[i] <- yv[i]/nv[i]
    
    Ds[i] <- 2*nv[i]*(prop[i])*log((prop[i]+0.00001)/p[i]) + (1-prop[i])*log((1-prop[i]+0.00001)/(1-p[i]))
    
  }
  
  # calculate saturated deviance
  dev.sat <- sum(Ds[])
  
  mean.data <- mean(yv)
  mean.sim <- mean(yv.sim)
  p.mean <- step(mean.sim - mean.data)
  
}

