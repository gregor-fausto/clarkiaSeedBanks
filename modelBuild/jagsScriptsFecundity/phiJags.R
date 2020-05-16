
model { 
  # hyperpriors
  for(i in 1:nsites){
    alphaP[i] ~ dnorm(0, .001)
  }
  
  for(i in 1:nyears){
    betaP[i] ~ dnorm(0, .001)
  }
  
  for(j in 1:nsites){
    for(k in 1:nyears){
      gammaP[j,k] ~ dnorm(0, 0.001)
      rP[j,k] ~ dgamma(.001,.001)
    }
  }
  
  
  # likelihoods
  for (i in 1:N2){
  log(lambdaP[i]) = alphaP[sitePhi[i]] + betaP[yearPhi[i]] + gammaP[sitePhi[i],yearPhi[i]]
  y2[i] ~ dnegbin(rP[sitePhi[i],yearPhi[i]]/(rP[sitePhi[i],yearPhi[i]]+lambdaP[i]),rP[sitePhi[i],yearPhi[i]])

  # simulated data for posterior predictive checks
  y2.sim[i] ~ dnegbin(rP[sitePhi[i],yearPhi[i]]/(rP[sitePhi[i],yearPhi[i]]+lambdaP[i]),rP[sitePhi[i],yearPhi[i]])
  }

}

