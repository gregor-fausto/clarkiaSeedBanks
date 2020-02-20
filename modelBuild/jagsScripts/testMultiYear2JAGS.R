model { 
  ## Priors
  ##hyperprior for intercept alpha
  
  for(i in 1:nsites){
    sigma.site[i] ~ dunif(0,100)
    tau.site[i] <- 1/(sigma.site[i]*sigma.site[i])
    mu.alpha[i] ~ dnorm(0, tau.site[i])
    for(k in 1:nyears){
      beta.i[i,k] ~ dnorm(0, 0.0001)
    }
  }
  
  # sum to zero constraint
  for(i in 1:nsites){
    mu_adj[i] = mu.alpha[i] + mean(beta.i[i,])
    for(k in 1:nyears){
      beta.i_adj[i,k] = beta.i[i,k] - mean(beta.i[i,])
      }
  }
  
  ## Likelihood
  for(i in 1:N_burial){
    pi[i] <- ilogit(mu.alpha[site[i]] + beta.i[site[i],year[i]])
    y_total[i] ~ dbin(pi[i], n_buried[i])
    yTotalSim[i] ~ dbin(pi[i], n_buried[i])
  }
  
  for(i in 1:nsites){
    for(j in 1:nyears){
      pred[i,j] <- mu.alpha[i] + beta.i[i,j]
      pred_adj[i,j] <- mu_adj[i] + beta.i_adj[i,j]
    }
  }
}

