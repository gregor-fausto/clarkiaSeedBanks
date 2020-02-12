model { 
 
  # hyperpriors
  
   theta ~ dunif(0,1)
  # kappa ~ dpar(alpha=shape,c=scale)
   kappa ~ dpar(1.5,1)
  
  # priors
  
  for(j in 1:nbags){
    p[j] ~ dbeta(kappa*theta, kappa*(1-theta))
  }
  
  # for(j in 1:nbags){
  #   p[j] ~ dbeta(1, 1)
  # }
  
  for(i in 1:N){
    
    # v viability
    yv[i] ~ dbinom( p[bag[i]] , nv[i] )

    yv.sim[i] ~ dbinom( p[bag[i]], nv[i] )

    # # code for deviance from Lunn 2013
    prop[i] <- yv[i]/nv[i]

    Ds[i] <- 2*nv[i]*(prop[i])*log((prop[i]+0.00001)/p[bag[i]]) + (1-prop[i])*log((1-prop[i]+0.00001)/(1-p[bag[i]]))
    
  }
  
  # calculate saturated deviance
  dev.sat <- sum(Ds[])
  
  mean.data <- mean(yv)
  mean.sim <- mean(yv.sim)
  p.mean <- step(mean.sim - mean.data)
  
}