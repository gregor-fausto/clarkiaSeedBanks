model { 
## Priors
##hyperprior for intercept alpha
  
  mu.alpha ~ dnorm(0, 0.0001)

  sigma.plot ~ dunif(0,100) #notated as varsigma in model documentation
  tau.plot <- 1 / sigma.plot^2
  
   for(i in 1:N.plots){     
     eps.plot[i]~dnorm(mu.alpha,tau.plot)
   }
  
## Likelihood
for(i in 1:N){
    p[i] <- ilogit(eps.plot[plot[i]])
  y[i] ~ dbin(p[i], n[i])

}
}


