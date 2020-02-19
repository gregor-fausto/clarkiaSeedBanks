model{
  
  ## Priors
  ##hyperprior for intercept alpha
  surv.mu~dnorm(0,0.001) 
  
  ## plot variances
  surv.sigma.plot~dunif(0,1000)
  surv.tau.plot<-1/(surv.sigma.plot*surv.sigma.plot)
  
  ## year variances
  surv.sigma.year~dunif(0,1000)
  surv.tau.year<-1/(surv.sigma.year*surv.sigma.year)
  
  ## prior for size slope - also do SVS
  # surv.bsize~dnorm(0,surv.precsize) 
  # surv.precsize <- 1/surv.Vsize
  # surv.Vsize <- (1-surv.zsize)*0.001 + surv.zsize*10
  # surv.zsize ~ dbern(0.5)
  
  ##plot and year deviates on the intercept
  for(i in 1:N.plots){      
    surv.eps.plot[i]~dnorm(0,surv.tau.plot)
  }
  for(i in 1:N.years){      
    surv.eps.year[i]~dnorm(0,surv.tau.year)
  }
  
  ## Likelihoods

  ## Survival (same model as flowering)
  for(i in 1:surv.N.obs){
    logit(surv.p[i]) <- surv.mu + surv.eps.plot[surv.plot[i]] + surv.eps.year[surv.year[i]]# + surv.bsize*surv.size[i]
    
    surv.y[i]~dbern(surv.p[i])
    
    surv.Presi[i] <- abs(surv.y[i]-surv.p[i])
    surv.y.new[i] ~ dbern(surv.p[i])
    surv.Presi.new[i] <- abs(surv.y.new[i] - surv.p[i])
  }
  
  
  ## sum up posterior predictive checks
  surv.fit <- sum(surv.Presi[]) 
  surv.fit.new <- sum(surv.Presi.new[]) 

  
}