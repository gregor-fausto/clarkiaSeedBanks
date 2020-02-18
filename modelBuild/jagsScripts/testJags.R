model { 
## Priors
##hyperprior for intercept alpha
  for(i in 1:nsites){
mu.alpha[i]~dnorm(0,0.001)   
  
sigma.year[i]~dunif(0,100)
tau.year[i]<-1/(sigma.year[i]*sigma.year[i])
  }
  
##plot and year deviates on the intercept
for(i in 1:nsiteyears){      
  eps.year[i]~dnorm(0,tau.year[site_year[i]])
}

## Likelihood
for(i in 1:N_burial){
  pi[i] <- ilogit(mu.alpha[site[i]] + eps.year[siteyear[i]])
  
  y_total[i] ~ dbin(pi[i], n_buried[i])

}


}
