model { 

    # priors
    for(k in 1:n_site){

        # theta 1
          mu0_1[k] ~  dnorm(0, 0.001)
          sigma0_1[k] ~ dnorm(0, 0.3) T(0,)
        tau0_1[k] <- 1/(sigma0_1[k]*sigma0_1[k])
                 
        for(i in 1:n_year){

            # theta 1
      mu_1[k,i] ~ dnorm(mu0_1[k], tau0_1[k])
      sigma_1[k,i] ~ dnorm(0, 0.3) T(0,)
            tau_1[k,i] <- 1/(sigma_1[k,i]*sigma_1[k,i])

    }
  }

    
    # Likelihood
      
    for(i in 1:n){
        
      # alpha
      alpha_1[i] ~ dnorm(mu_1[site[i],year[i]],tau_1[site[i],year[i]])
      
      # logit 
      logit(theta_1[i]) <- alpha_1[i]
      
       # likelihood
        fruitplNumber[i] ~ dbinom(theta_1[i], seedlingNumber[i])

        # posterior predictive
        fruitplNumber.sim[i] ~ dbinom(theta_1[i], seedlingNumber[i])
        
        
        # https://people.ucsc.edu/~abrsvn/bayes_winbugs_jags_5.r
        #Pearson.res[i] <- (fruitplNumber[i] - (seedlingNumber[i]*theta_1[i]))/(sqrt(seedlingNumber[i]*theta_1[i]*(1-theta_1[i])))
        #Pearson.res.sim[i] <- (fruitplNumber.sim[i] - (seedlingNumber[i]*theta_1[i]))/(sqrt(seedlingNumber[i]*theta_1[i]*(1-theta_1[i])))

    }

    for(i in 1:n_site){
      
      # mean.y <- mean(fruitplNumber[i])
      # mean.y.sim <- mean(fruitplNumber.sim[i])
      # p.value.mean <- step(mean.y.sim-mean.y)
      # 
      # sd.y <- sd(fruitplNumber[i])
      # sd.y.sim <- sd(fruitplNumber.sim[i])
      # p.value.sd <- step(sd.y.sim-sd.y)
      
         p.i0_1[i] ~ dnorm(mu0_1[i],tau0_1[i])
    logit(p0_1[i]) <- p.i0_1[i]
        
        for(j in 1:n_year){
            
       p.i_1[i,j] ~ dnorm(mu_1[i,j],tau_1[i,j])
    logit(p_1[i,j]) <- p.i_1[i,j]

        }
       }
    
}
