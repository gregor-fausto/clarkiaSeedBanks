model { 
## Priors
##hyperprior for intercept alpha
  
  mu.alpha ~ dnorm(0, 0.0001)
  
  sigma.plot~dunif(0,1000)
  tau.plot<-1/(sigma.plot*sigma.plot)

  for(i in 1:N.plots){      
    eps.plot[i]~dnorm(0,tau.plot)
  }
  
## Likelihood
for(i in 1:N){
    p[i] <- ilogit(mu.alpha+eps.plot[plot[i]])
  y[i] ~ dbin(p[i], n[i])

}
}


