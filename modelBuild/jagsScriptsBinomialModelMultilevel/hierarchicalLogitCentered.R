model { 

  # priors

    for(k in 1:n_site){
          mu0[k] ~  dnorm(0, 0.001)
          sigma0[k] ~ dunif(0,100)
          tau0[k] <- 1/(sigma0[k]*sigma0[k])
    for(i in 1:n_year){
      mu[k,i] ~ dnorm(mu0[k], tau0[k])
      sigma[k,i] ~ dunif(0,100)
      tau[k,i] <- 1/(sigma[k,i]*sigma[k,i])
      
    }
  }

    
  # Likelihood
        for(i in 1:n){
            
            alpha[i] ~ dnorm(mu[site[i],year[i]],tau[site[i],year[i]])
            logit(theta[i]) <- alpha[i]
            
            totalJan[i] ~ dbinom(theta[i], seedStart[i])
    }

   for(i in 1:n_site){
    p.i[i] ~ dnorm(mu0[i],tau0[i])
    logit(p[i]) <- p.i[i]
    }
    
}

