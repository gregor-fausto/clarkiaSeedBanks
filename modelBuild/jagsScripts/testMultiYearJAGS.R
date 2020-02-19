model { 
  ## Priors
  ##hyperprior for intercept alpha
  
  for(i in 1:nsites){
    mu.alpha[i] ~ dnorm(0, 0.0001)
    for(k in 1:nyears){
      beta.i[i,k] ~ dnorm(0, 0.0001)
      }
  }
  
  # sum to zero constraint
  for(i in 1:nsites){
    alpha.i[i,1] <- 0 - sum( beta.i[i,2:nyears] )
    alpha.i[i,2:nyears] <- beta.i[i,2:nyears]
  }
  
  ## Likelihood
  for(i in 1:N_burial){
    pi[i] <- ilogit(mu.alpha[site[i]] + alpha.i[site[i],year[i]])
    y_total[i] ~ dbin(pi[i], n_buried[i])
    
  }
  
  for(i in 1:nsites){
    for(j in 1:nyears){
      pred[i,j] <- mu.alpha[i] + alpha.i[i,j]
    }
  }
}

