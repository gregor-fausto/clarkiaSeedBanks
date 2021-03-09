
model {
  
  # PRIORS -----------------------------------------------
  
  # for(i in 1:n_site){
  #   
  #   # TOTAL FRUIT EQUIVALENTS -----------------------------------------
  #   
  #   mu0_tfe[i] ~  dnorm(1, 1)
  #   
  #   # sigma0_tfe[i] ~ dexp(2)
  #   sigma0_tfe[i] ~ dnorm(0, 5) T(0,)
  #   tau0_tfe[i] <- 1/(sigma0_tfe[i]*sigma0_tfe[i])   
  #   
  #   #tau0_tfe[i] ~ dgamma(0.001,0.001) # This is inverse gamma
  #   #  sigma0_tfe[i] <- 1/sqrt( tau0_tfe[i]); # sd is treated as derived parameter
  #   
  #   for(j in 1:n_year){
  #     
  #     # theta 1
  #     mu_tfe[i,j] ~ dnorm(mu0_tfe[i], tau0_tfe[i])
  #     
  #     # tau_tfe[i,j] ~ dgamma(0.001,0.001); # This is inverse gamma
  #     #sigma_tfe[i,j] <- 1/sqrt( tau_tfe[i,j]); # sd is treated as derived parameter
  #     
  #     # sigma_tfe[i,j] ~ dexp(2)
  #     sigma_tfe[i,j] ~ dnorm(0, 5) T(0,)
  #     tau_tfe[i,j] <- 1/(sigma_tfe[i,j]*sigma_tfe[i,j])   
  #     
  #     #kappa_tfe[i,j] ~ dt(0,1/5^2,1) I(0,)
  #     #  kappa_tfe[i,j] ~ dnorm(0,1) T(0,)
  #     #  kappa_tfe[i,j] ~ dgamma(2,.1)
  #     
  #     #  mu_tfe[i,j] ~ dlnorm(log(g_tfe[i]), tau0_tfe[i])
  #     # # tau_tfe[i,j] ~ dgamma(0.001,0.001)  # Inverse gamma is not a built in distribution
  #     #  sigma_tfe[i,j] ~ dnorm(0, 1) T(0,)
  #     #  tau_tfe[i,j] <- 1/(sigma_tfe[i,j]*sigma_tfe[i,j])
  #   }
  # }
  
  # for(i in 1:n_site){
  #   
  #   # TOTAL FRUIT EQUIVALENTS -----------------------------------------
  #   nu_tfe_prior[i] ~ dgamma(2, 2)
  #   nu_tfe[i] ~ dgamma(2, 2)
  #   g_tfe[i] = exp(nu_tfe[i])
  #   
  #   # page 76 https://nwfsc-timeseries.github.io/atsa-2017/Labs/Week%203%20intro%20to%20jags/intro-to-jags.pdf
  #  # tau0_tfe[i] ~ dgamma(0.001,0.001)  # Inverse gamma is not a built in distribution
  #   sigma0_tfe_prior[i] ~ dnorm(0, 1) T(0,)
  #     sigma0_tfe[i] ~ dnorm(0, 1) T(0,)
  #     tau0_tfe[i] <- 1/(sigma0_tfe[i]*sigma0_tfe[i])
  #     
  #   for(j in 1:n_year){
  #     
  #     mu_tfe[i,j] ~ dlnorm(log(g_tfe[i]), tau0_tfe[i])
  #   #  tau_tfe[i,j] ~ dgamma(0.001,0.001)  # Inverse gamma is not a built in distribution
  #     sigma_tfe_prior[i,j] ~ dnorm(0, 1) T(0,)
  #     sigma_tfe[i,j] ~ dnorm(0, 1) T(0,)
  #       tau_tfe[i,j] <- 1/(sigma_tfe[i,j]*sigma_tfe[i,j])
  #     
  #   }
  # }
  
  for(i in 1:n_site){
    
    # DAMAGED FRUITS PER PLANT  -----------------------------------------
    
    nu_dam[i] ~ dgamma(1, 1)
    g_dam[i] = exp(nu_dam[i])
    
    # page 76 https://nwfsc-timeseries.github.io/atsa-2017/Labs/Week%203%20intro%20to%20jags/intro-to-jags.pdf
    #  tau0_dam[i] ~ dgamma(0.001,0.001)  # Inverse gamma is not a built in distribution
    sigma0_dam[i] ~ dnorm(0, 1) T(0,)
    tau0_dam[i] <- 1/(sigma0_dam[i]*sigma0_dam[i])
    
    for(j in 1:n_year2){
      
      mu_dam[i,j] ~ dlnorm(log(g_dam[i]), tau0_dam[i])
      #      tau_dam[i,j] ~ dgamma(0.001,0.001)  # Inverse gamma is not a built in distribution
      sigma_dam[i,j] ~ dnorm(0, 1) T(0,)
      tau_dam[i,j] <- 1/(sigma_dam[i,j]*sigma_dam[i,j])
    }
  }
  
  
#   for(i in 1:n_site){
# 
#     # SEEDS PER DAMAGED FRUITS  -----------------------------------------
# 
#     nu_dam_seeds[i] ~ dgamma(2, 2)
#     g_dam_seeds[i] = exp(nu_dam_seeds[i])
# 
#     # page 76 https://nwfsc-timeseries.github.io/atsa-2017/Labs/Week%203%20intro%20to%20jags/intro-to-jags.pdf
#   #  tau0_dam_seeds[i] ~ dgamma(0.001,0.001)  # Inverse gamma is not a built in distribution
#         sigma0_dam_seeds[i] ~ dnorm(0, 1) T(0,)
#         tau0_dam_seeds[i] <- 1/(sigma0_dam_seeds[i]*sigma0_dam_seeds[i])
#         
#     for(j in 1:n_year4){
# 
#       mu_dam_seeds[i,j] ~ dlnorm(log(g_dam_seeds[i]), tau0_dam_seeds[i])
# #      tau_dam_seeds[i,j] ~ dgamma(0.001,0.001)  # Inverse gamma is not a built in distribution
#       sigma_dam_seeds[i,j] ~ dnorm(0, 1) T(0,)
#       tau_dam_seeds[i,j] <- 1/(sigma_dam_seeds[i,j]*sigma_dam_seeds[i,j])
#     }
#   }

  
  # LIKELIHOODS -------------------------------------------------------------
  
  # for (i in 1:n){
  #   # alpha_tfe[i] ~ dnorm(mu_tfe[site[i],year[i]],tau_tfe[site[i],year[i]])
  #   # log(lambda_tfe[i]) = alpha_tfe[i]
  #   # 
  #   
  #   lambda_tfe[i] ~ dlnorm(log(mu_tfe[site[i],year[i]]),tau_tfe[site[i],year[i]])
  #   # y_tfe[i] ~ dpois(lambda_tfe[i])
  #    y_prior[i] ~ dpois(lambda_tfe[i])
  #   
  # }
  
  for (i in 1:n2){
    #  DAMAGED FRUITS -------------------------------------------------
    z_dam[i] ~ dlnorm(log(mu_dam[site2[i],year2[i]]),tau_dam[site2[i],year2[i]])
    y_dam[i] ~ dpois(z_dam[i])
  }
  
  # SEEDS PER DAMAGED FRUITS -------------------------------------------------
  # for (i in 1:n4){
  #   z_sd_dam[i] ~ dlnorm(log(mu_dam_seeds[site4[i],year4[i]]),tau_dam_seeds[site4[i],year4[i]])
  #   sdno_dam[i] ~ dpois(z_sd_dam[i])
  #   }


  
  
}
