model { 

 # priors
    for(k in 1:n_siteViab){

        # theta 1
          mu0_g[k] ~  dnorm(0, 0.001)
          sigma0_g[k] ~ dnorm(0,1.0E-3) T(0,)
        tau0_g[k] <- 1/(sigma0_g[k]*sigma0_g[k])

                                        # theta 2
          mu0_v[k] ~  dnorm(0, 0.001)
          sigma0_v[k] ~ dnorm(0,1.0E-3) T(0,)
        tau0_v[k] <- 1/(sigma0_v[k]*sigma0_v[k])
        
        for(i in 1:n_yearViab){

            # theta 1
      mu_g[k,i] ~ dnorm(mu0_g[k], tau0_g[k])
      sigma_g[k,i] ~ dnorm(0,1.0E-3) T(0,)
            tau_g[k,i] <- 1/(sigma_g[k,i]*sigma_g[k,i])

            # theta 2
            mu_v[k,i] ~ dnorm(mu0_v[k], tau0_v[k])
      sigma_v[k,i] ~ dnorm(0,1.0E-3) T(0,)
            tau_v[k,i] <- 1/(sigma_v[k,i]*sigma_v[k,i])
    }
  }

                                        # priors
    for(k in 1:n_siteBags){

        # theta 1
          mu0_1[k] ~  dnorm(0, 0.001)
          sigma0_1[k] ~ dnorm(0,1.0E-3) T(0,)
        tau0_1[k] <- 1/(sigma0_1[k]*sigma0_1[k])

                                        # theta 2
          mu0_2[k] ~  dnorm(0, 0.001)
          sigma0_2[k] ~ dnorm(0,1.0E-3) T(0,)
        tau0_2[k] <- 1/(sigma0_2[k]*sigma0_2[k])
        
        # theta 3
        mu0_3[k] ~  dnorm(0, 0.001)
        sigma0_3[k] ~ dnorm(0,1.0E-3) T(0,)
        tau0_3[k] <- 1/(sigma0_3[k]*sigma0_3[k])
        
        for(i in 1:n_yearBags){

            # theta 1
      mu_1[k,i] ~ dnorm(mu0_1[k], tau0_1[k])
      sigma_1[k,i] ~ dnorm(0,1.0E-3) T(0,)
            tau_1[k,i] <- 1/(sigma_1[k,i]*sigma_1[k,i])

            # theta 2
            mu_2[k,i] ~ dnorm(mu0_2[k], tau0_2[k])
      sigma_2[k,i] ~ dnorm(0,1.0E-3) T(0,)
            tau_2[k,i] <- 1/(sigma_2[k,i]*sigma_2[k,i])
            
            # theta 3
            mu_3[k,i] ~ dnorm(mu0_3[k], tau0_3[k])
            sigma_3[k,i] ~ dnorm(0,1.0E-3) T(0,)
            tau_3[k,i] <- 1/(sigma_3[k,i]*sigma_3[k,i])
    }
  }

    
                                        # Likelihood
    for(i in 1:n_bag){

                
      # alpha
      alpha_g[i] ~ dnorm(mu_g[siteViab[i],yearViab[i]],tau_g[siteViab[i],yearViab[i]])
      alpha_v[i] ~ dnorm(mu_v[siteViab[i],yearViab[i]],tau_v[siteViab[i],yearViab[i]])

      # logit 
      logit(theta_g[i]) <- alpha_g[i]
        logit(theta_v[i]) <- alpha_v[i]
        
                germCount[i] ~ dbinom( theta_g[i] , germStart[i] )
                viabStain[i] ~ dbinom( theta_v[i] , viabStart[i] )

                }
    
    for(i in 1:n){
        
      # alpha
      alpha_1[i] ~ dnorm(mu_1[siteBags[i],yearBags[i]],tau_1[siteBags[i],yearBags[i]])
      alpha_2[i] ~ dnorm(mu_2[siteBags[i],yearBags[i]],tau_2[siteBags[i],yearBags[i]])
      alpha_3[i] ~ dnorm(mu_3[siteBags[i],yearBags[i]],tau_3[siteBags[i],yearBags[i]])
      
      # logit 
      logit(theta_1[i]) <- alpha_1[i]
      logit(theta_2[i]) <- alpha_2[i]
      logit(theta_3[i]) <- alpha_3[i]
      
       # likelihood
       totalJan[i] ~ dbinom(theta_1[i], seedStart[i])
        seedlingJan[i] ~ dbinom(theta_2[i], totalJan[i])
        intactJan[i] = totalJan[i]-seedlingJan[i]
        intactOct[i] ~ dbinom(theta_3[i], intactJan[i])
        
    }

   for(i in 1:n_siteBags){

    p.i_g[i] ~ dnorm(mu0_g[i],tau0_g[i])
    logit(p_g[i]) <- p.i_g[i]

    p.i_v[i] ~ dnorm(mu0_v[i],tau0_v[i])
       logit(p_v[i]) <- p.i_v[i]
       
    nu_1[i] = p_g[i] + p_v[i]*(1-p_g[i])

       p.i_1[i] ~ dnorm(mu0_1[i],tau0_1[i])
    logit(p_1[i]) <- p.i_1[i]

    p.i_2[i] ~ dnorm(mu0_2[i],tau0_2[i])
    logit(p_2[i]) <- p.i_2[i]
    
    p.i_3[i] ~ dnorm(mu0_3[i],tau0_3[i])
    logit(p_3[i]) <- p.i_3[i]

       s1[i] = p_1[i]*(p_2[i] + (1-p_2[i])*(nu_1[i])^(1/3))
       g1[i] = p_2[i]/(1-(1-(nu_1[i]^(1/3)))*(1-p_2[i]))
       s2[i] = p_3[i]*(nu_1[i]^(2/3))
   
  
  for(k in 1:n_yearBags){
    
    p.i0_g[i,k] ~ dnorm(mu_g[i,k],tau_g[i,k])
    logit(p0_g[i,k]) <- p.i0_g[i,k]
    
    p.i0_v[i,k] ~ dnorm(mu_v[i,k],tau_v[i,k])
    logit(p0_v[i,k]) <- p.i0_v[i,k]
    
    nu0_1[i,k] = p0_g[i,k] + p0_v[i,k]*(1-p0_g[i,k])
    
    p.i0_1[i,k] ~ dnorm(mu_1[i,k],tau_1[i,k])
    logit(p0_1[i,k]) <- p.i0_1[i,k]
    
    p.i0_2[i,k] ~ dnorm(mu_2[i,k],tau_2[i,k])
    logit(p0_2[i,k]) <- p.i0_2[i,k]
  
    p.i0_3[i,k] ~ dnorm(mu_3[i,k],tau_3[i,k])
    logit(p0_3[i,k]) <- p.i0_3[i,k]
    
    s1.0[i,k] = p0_1[i,k]*(p0_2[i,k] + (1-p0_2[i,k])*(nu0_1[i,k])^(1/3))
    g1.0[i,k] = p0_2[i,k]/(1-(1-(nu0_1[i,k]^(1/3)))*(1-p0_2[i,k]))
    s2.0[i,k] = p0_3[i,k]*(nu0_1[i,k]^(2/3))
    
    
  }
}
    
} 

