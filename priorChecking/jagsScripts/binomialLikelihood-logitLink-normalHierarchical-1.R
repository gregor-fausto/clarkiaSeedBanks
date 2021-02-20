model { 
  
  # priors
  for(k in 1:n_site){
    
    # theta 1
    mu0[k] ~  dnorm(0, 0.001)
    sigma0[k] ~ dnorm(0, 0.3) T(0,)
    tau0[k] <- 1/(sigma0[k]*sigma0[k])
    
    for(i in 1:n_year){
      
      # theta 1
      mu[k,i] ~ dnorm(mu0[k], tau0[k])
      sigma[k,i] ~ dnorm(0, 0.3) T(0,)
      tau[k,i] <- 1/(sigma[k,i]*sigma[k,i])
      
    }
  }
  
  
  # Likelihood
  
  for(i in 1:n){
    
    # alpha
    alpha[i] ~ dnorm(mu[site[i],year[i]],tau[site[i],year[i]])
    
    # logit 
    logit(theta[i]) <- alpha[i]
    
    # Prior predictive
    y_prior[i] ~ dbinom(theta[i], seedlingNumber[i])
    
    # likelihood
    # fruitplNumber[i] ~ dbinom(theta[i], seedlingNumber[i])
    
    # posterior predictive
    fruitplNumber.sim[i] ~ dbinom(theta[i], seedlingNumber[i])
    
    
    # https://people.ucsc.edu/~abrsvn/bayes_winbugs_jags_5.r
    #Pearson.res[i] <- (fruitplNumber[i] - (seedlingNumber[i]*theta_1[i]))/(sqrt(seedlingNumber[i]*theta_1[i]*(1-theta_1[i])))
    #Pearson.res.sim[i] <- (fruitplNumber.sim[i] - (seedlingNumber[i]*theta_1[i]))/(sqrt(seedlingNumber[i]*theta_1[i]*(1-theta_1[i])))
    
  }
  
  # for(i in 1:n_site){
  #   
  #   # mean.y <- mean(fruitplNumber[i])
  #   # mean.y.sim <- mean(fruitplNumber.sim[i])
  #   # p.value.mean <- step(mean.y.sim-mean.y)
  #   # 
  #   # sd.y <- sd(fruitplNumber[i])
  #   # sd.y.sim <- sd(fruitplNumber.sim[i])
  #   # p.value.sd <- step(sd.y.sim-sd.y)
  #   
  #   p.i0_1[i] ~ dnorm(mu0_1[i],tau0_1[i])
  #   logit(p0_1[i]) <- p.i0_1[i]
  #   
  #   for(j in 1:n_year){
  #     
  #     p.i_1[i,j] ~ dnorm(mu_1[i,j],tau_1[i,j])
  #     logit(p_1[i,j]) <- p.i_1[i,j]
  #     
  #   }
  # }
  
}
