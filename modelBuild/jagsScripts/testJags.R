model { 
## Priors
##hyperprior for intercept alpha
mu.alpha~dnorm(0,0.001)   

sigma.year~dunif(0,1000)
tau.year<-1/(sigma.year*sigma.year)

##plot and year deviates on the intercept
for(i in 1:nsiteyears){      
  eps.year[i]~dnorm(0,tau.year)
}

## Likelihood
for(i in 1:N_burial){
  logit(pi[i]) <- mu.alpha + eps.year[siteyear[i]]
  
  y_total[i] ~ dbin(pi[i], n_buried[i])

}


}
