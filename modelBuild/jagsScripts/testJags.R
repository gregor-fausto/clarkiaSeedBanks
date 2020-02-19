model { 
## Priors
##hyperprior for intercept alpha
  
  mu.alpha ~ dnorm(0, 0.0001)

  for(j in 1:nsiteyears){
    beta.i[j] ~ dnorm(0, 0.0001)
  }
  
  # sum to zero constraint
  alpha.i[1] <- 0 - sum( beta.i[2:nsiteyears])
  alpha.i[2:nsiteyears] <- beta.i[2:nsiteyears]

## Likelihood
for(i in 1:N_burial){
    pi[i] <- ilogit(mu.alpha + alpha.i[siteyear[i]])
  y_total[i] ~ dbin(pi[i], n_buried[i])

}
}

