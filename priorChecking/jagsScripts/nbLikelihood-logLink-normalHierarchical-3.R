
model {

  # PRIORS -----------------------------------------------
  
  for(i in 1:n_site){

# TOTAL FRUIT EQUIVALENTS -----------------------------------------

    mu0_tfe[i] ~  dnorm(1, 1)
    
   # sigma0_tfe[i] ~ dexp(2)
    sigma0_tfe[i] ~ dnorm(0, 1) T(0,)
    tau0_tfe[i] <- 1/(sigma0_tfe[i]*sigma0_tfe[i])   
    
    #tau0_tfe[i] ~ dgamma(0.001,0.001) # This is inverse gamma
  #  sigma0_tfe[i] <- 1/sqrt( tau0_tfe[i]); # sd is treated as derived parameter
    
    for(j in 1:n_year){

      # theta 1
      mu_tfe[i,j] ~ dnorm(mu0_tfe[i], tau0_tfe[i])
      
     # tau_tfe[i,j] ~ dgamma(0.001,0.001); # This is inverse gamma
     #sigma_tfe[i,j] <- 1/sqrt( tau_tfe[i,j]); # sd is treated as derived parameter
      
     # sigma_tfe[i,j] ~ dexp(2)
      sigma_tfe[i,j] ~ dnorm(0, 1) T(0,)
      tau_tfe[i,j] <- 1/(sigma_tfe[i,j]*sigma_tfe[i,j])   
      
      #kappa_tfe[i,j] ~ dt(0,1/5^2,1) I(0,)
    #  kappa_tfe[i,j] ~ dnorm(0,1) T(0,)
    #  kappa_tfe[i,j] ~ dgamma(2,.1)
      
     #  mu_tfe[i,j] ~ dlnorm(log(g_tfe[i]), tau0_tfe[i])
     # # tau_tfe[i,j] ~ dgamma(0.001,0.001)  # Inverse gamma is not a built in distribution
     #  sigma_tfe[i,j] ~ dnorm(0, 1) T(0,)
     #  tau_tfe[i,j] <- 1/(sigma_tfe[i,j]*sigma_tfe[i,j])
    }
  }

# LIKELIHOODS -------------------------------------------------------------

  for (i in 1:n){
 alpha_tfe[i] ~ dnorm(mu_tfe[site[i],year[i]],tau_tfe[site[i],year[i]])
 log(lambda_tfe[i]) = alpha_tfe[i]

 y_prior[i] ~ dpois(lambda_tfe[i])
      }
  


}
