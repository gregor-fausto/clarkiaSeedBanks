model { 
## Priors
##hyperprior for intercept alpha
  
  mu.alpha ~ dnorm(0, 0.0001)

  sigma.plot ~ dunif(0,100) #notated as varsigma in model documentation
  tau.plot <- 1 / sigma.plot^2
  
  # sigma.block ~ dunif(0,100) #notated as varsigma in model documentation
  # tau.block <- 1 / sigma.block^2
  
   for(i in 1:N.plots){     
     eps.plot[i]~dnorm(0,tau.plot)
   }
  
  # for(i in 1:N.blocks){     
  #   eps.block[i]~dnorm(0,tau.block)
  # }
  
## Likelihood
for(i in 1:N){
    logit(p[i]) <- mu.alpha + eps.plot[plot[i]]#+eps.block[block[i]]
    y[i] ~ dbin(p[i], n[i])

}
}


