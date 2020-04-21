

model{
  
    # hyperpriors
    alpha ~ dgamma(0.001, 0.001)
    upsilon ~ dgamma(0.001, 0.001)
    
    # priors
    # place a (0,.5) prior on pi
    # assuming that > half of seedlings are counted
    # as parameter decreases, true number of seedlings = observed number
    # as parameter increases, true number of seedlings > observed number (max p = 0)
    # as parameter = .5 half of seeldings are observed
    for(i in 1:n_transect){ pi[i] ~ dbeta(1,1) }
    phi ~ dbeta(1, 1) 
    for(i in 1:n_transect){ lambda[i] ~ dgamma(alpha, upsilon) } 
    
    # likelihood
    for (i in 1:n){ 
    # process model for true number of seedlings drawn from lambda
    z[i] ~ dpois(lambda[transect[i]])
    
    # data model
    # counts of seedling number related to true number of seedlings (z[i])
    # via detection probability (1-pi)
    seedlingNumber[i] ~ dbinom(1-pi[transect[i]], z[i])
    
    # data model 
    # counts of fruiting plants relating true number of seedlings (z[i])
    # via survival probability (phi)
    fruitplNumber[i] ~ dbinom(phi, z[i])
    } 
    
  # end of data model



} #end of model


