model { 
## Priors
##hyperprior for intercept alpha
  
  mu.alpha ~ dnorm(0, 0.0001)
  
    sigma.year~dunif(0,100)
    tau.year<-1/(sigma.year*sigma.year)

  for(j in 1:nsiteyears){
    beta.i[j] ~ dnorm(0, tau.year)
  }

## Likelihood
for(i in 1:N_burial){
    pi[i] <- ilogit(mu.alpha + beta.i[siteyear[i]])
  y_total[i] ~ dbin(pi[i], n_buried[i])

}
}


