model { 
  
  # priors
  for(k in 1:n_site){
    
    # hyperpriors -------------------------------------------------------------
    mu0[k] ~  dnorm(0, 1)
    sigma0[k] ~ dnorm(0, 1) T(0,)
    tau0[k] <- 1/(sigma0[k]*sigma0[k])
    
    for(i in 1:n_year){
      
      # priors ------------------------------------------------------------------
      mu[k,i] ~ dnorm(mu0[k], tau0[k])
      # weakly informative on sigma
      sigma[k,i] ~ dnorm(0, 1) T(0,)
      tau[k,i] <- 1/(sigma[k,i]*sigma[k,i])
      
      # pooling factor ------------------------------------------------------------------
      # e.mu[k,i] <- mu[k,i]-mu0[k]
    }
  }
  
  # Likelihood ------------------------------------------------------------------
  
  for(i in 1:n){
    
    # alpha
    alpha[i] ~ dnorm(mu[site[i],year[i]],tau[site[i],year[i]])
    
    # logit 
    logit(theta[i]) <- alpha[i]
    
    # Prior predictive
    #  y_prior[i] ~ dbinom(theta[i], seedlingNumber[i])
    
    # LIKELIHOOD
    fruitplNumber[i] ~ dbinom(theta[i], seedlingNumber[i])
    
    # POSTERIOR PREDICTIVE
    fruitplNumber_sim[i] ~ dbinom(theta[i], seedlingNumber[i])
    
    chi2.obs[i] <- pow((fruitplNumber[i]- theta[i]*seedlingNumber[i]),2) / (theta[i]*seedlingNumber[i]+.001)         # obs.
    chi2.sim[i] <- pow((fruitplNumber_sim[i]- theta[i]*seedlingNumber[i]),2) / (theta[i]*seedlingNumber[i]+.001)         # obs.
    
    ## LOG-LIKELIHOOD
    logLik[i] <- logdensity.bin(fruitplNumber[i],theta[i],seedlingNumber[i])
    
    # https://people.ucsc.edu/~abrsvn/bayes_winbugs_jags_5.r
    #Pearson.res[i] <- (fruitplNumber[i] - (seedlingNumber[i]*theta_1[i]))/(sqrt(seedlingNumber[i]*theta_1[i]*(1-theta_1[i])))
    #Pearson.res.sim[i] <- (fruitplNumber.sim[i] - (seedlingNumber[i]*theta_1[i]))/(sqrt(seedlingNumber[i]*theta_1[i]*(1-theta_1[i])))
    
  }
  
}
