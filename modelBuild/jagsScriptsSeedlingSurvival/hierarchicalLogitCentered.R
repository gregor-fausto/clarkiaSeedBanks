model { 

    # priors
    for(k in 1:n_site){

        # theta 1
          mu0_1[k] ~  dnorm(0, 0.001)
          sigma0_1[k] ~ dunif(0,1.5)
        tau0_1[k] <- 1/(sigma0_1[k]*sigma0_1[k])

                 
        for(i in 1:n_year){

            # theta 1
      mu_1[k,i] ~ dnorm(mu0_1[k], tau0_1[k])
      sigma_1[k,i] ~ dunif(0,1.5)
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
        
    }

    for(i in 1:n_site){

         p.i0_1[i] ~ dnorm(mu0_1[i],tau0_1[i])
    logit(p0_1[i]) <- p.i0_1[i]
        
        for(j in 1:n_year){
            
       p.i_1[i,j] ~ dnorm(mu_1[i,j],tau_1[i,j])
    logit(p_1[i,j]) <- p.i_1[i,j]

        }
       }
    
}
