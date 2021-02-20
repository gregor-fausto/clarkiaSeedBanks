
model {

  # PRIORS ------------------------------------------------

  
  # for(i in 1:n_site){
  #   
  #   # UNDAMAGED FRUITS
  #   # theta 1
  #   mu0_und[i] ~  dnorm(0, 1)
  #   sigma0_und[i] ~ dnorm(0,1) T(0,)
  #   tau0_und[i] <- 1/(sigma0_und[i]*sigma0_und[i])
  #   
  #   # DAMAGED FRUITS
  #   # theta 1
  #   # mu0_dam[i] ~  dnorm(0, 0.001)
  #   # sigma0_dam[i] ~ dunif(0,10)
  #   # tau0_dam[i] <- 1/(sigma0_dam[i]*sigma0_dam[i])
  #   
  #   for(j in 1:n_year2){
  #     
  #     # theta 1
  #     gamma_und[i,j] ~ dnorm(mu0_und[i], tau0_und[i])
  #     #r_und[i,j] ~ dgamma(.001,.001)
  #     r_und.inv[i,j] ~ dnorm(0,1) T(0,)
  #     r_und[i,j] = 1/r_und.inv[i,j]
  #     
  #     # theta 1
  #     # gamma_dam[i,j] ~ dnorm(mu0_dam[i], tau0_dam[i])
  #     # r_dam[i,j] ~ dgamma(.001,.001)
  #     
  #   }
  # }
  
  
  for(i in 1:n_site){

# TOTAL FRUIT EQUIVALENTS -----------------------------------------

    mu0_tfe[i] ~  dnorm(0, 1)
  #  mu0_tfe[i] ~ dt(0,1/2.5^2,1) I(0,)
    sigma0_tfe[i] ~ dnorm(0,2) T(0,)
    #sigma0_tfe[i] ~ dt(0,1,1) I(0,)
   # sigma0_tfe[i] ~ dexp(1)
    tau0_tfe[i] <- 1/(sigma0_tfe[i]*sigma0_tfe[i])
    
    # nu_tfe[i] ~ dgamma(1, 1)
    # g_tfe[i] = exp(nu_tfe[i])
    # 
    # # page 76 https://nwfsc-timeseries.github.io/atsa-2017/Labs/Week%203%20intro%20to%20jags/intro-to-jags.pdf
    # # tau0_tfe[i] ~ dgamma(0.001,0.001)  # Inverse gamma is not a built in distribution
    # sigma0_tfe[i] ~ dnorm(0, 1) T(0,)
    # tau0_tfe[i] <- 1/(sigma0_tfe[i]*sigma0_tfe[i])
    
    for(j in 1:n_year){

      # theta 1
      mu_tfe[i,j] ~ dnorm(mu0_tfe[i], tau0_tfe[i])
     # sigma_tfe[i,j] ~ dt(0,1,1) I(0,)
      sigma_tfe[i,j] ~ dnorm(0,2) T(0,)
     # sigma_tfe[i,j] ~ dexp(1)
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


  # for(i in 1:n_site2){
  # 
  # # UNDAMAGED/DAMAGED FRUITS ------------------------------------------------
  # 
  #   ## UNDAMAGED
  #   nu_und[i] ~ dgamma(1, 1)
  #   g_und[i] = exp(nu_und[i])
  # 
  #   # page 76 https://nwfsc-timeseries.github.io/atsa-2017/Labs/Week%203%20intro%20to%20jags/intro-to-jags.pdf
  #   tau0_und[i] ~ dgamma(0.001,0.001)  # Inverse gamma is not a built in distribution
  # 
  #   ## DAMAGED
  #   nu_dam[i] ~ dgamma(1, 1)
  #   g_dam[i] = exp(nu_dam[i])
  # 
  #   # page 76 https://nwfsc-timeseries.github.io/atsa-2017/Labs/Week%203%20intro%20to%20jags/intro-to-jags.pdf
  #   tau0_dam[i] ~ dgamma(0.001,0.001)  # Inverse gamma is not a built in distribution
  # 
  #   for(j in 1:n_year2){
  # 
  #     mu_und[i,j] ~ dlnorm(log(g_und[i]), tau0_und[i])
  #     tau_und[i,j] ~ dgamma(0.001,0.001)  # Inverse gamma is not a built in distribution
  # 
  #     mu_dam[i,j] ~ dlnorm(log(g_dam[i]), tau0_dam[i])
  #     tau_dam[i,j] ~ dgamma(0.001,0.001)  # Inverse gamma is not a built in distribution
  # 
  #   }
  # }
  # 
  # for(i in 1:n_site){
  # 
  #   # SEEDS PER UNDAMAGED FRUITS  -----------------------------------------
  # 
  #   nu_seeds[i] ~ dgamma(1, 1)
  #   g_seeds[i] = exp(nu_seeds[i])
  # 
  #   # page 76 https://nwfsc-timeseries.github.io/atsa-2017/Labs/Week%203%20intro%20to%20jags/intro-to-jags.pdf
  #   tau0_seeds[i] ~ dgamma(0.001,0.001)  # Inverse gamma is not a built in distribution
  # 
  #   for(j in 1:n_year3){
  # 
  #     mu_seeds[i,j] ~ dlnorm(log(g_seeds[i]), tau0_seeds[i])
  #     tau_seeds[i,j] ~ dgamma(0.001,0.001)  # Inverse gamma is not a built in distribution
  # 
  #   }
  # }
  # 
  # for(i in 1:n_site){
  # 
  #   # SEEDS PER DAMAGED FRUITS  -----------------------------------------
  # 
  #   nu_dam_seeds[i] ~ dgamma(1, 1)
  #   g_dam_seeds[i] = exp(nu_dam_seeds[i])
  # 
  #   # page 76 https://nwfsc-timeseries.github.io/atsa-2017/Labs/Week%203%20intro%20to%20jags/intro-to-jags.pdf
  #   tau0_dam_seeds[i] ~ dgamma(0.001,0.001)  # Inverse gamma is not a built in distribution
  # 
  #   for(j in 1:n_year4){
  # 
  #     mu_dam_seeds[i,j] ~ dlnorm(log(g_dam_seeds[i]), tau0_dam_seeds[i])
  #     tau_dam_seeds[i,j] ~ dgamma(0.001,0.001)  # Inverse gamma is not a built in distribution
  # 
  #   }
  # }

# LIKELIHOODS -------------------------------------------------------------

 # for (i in 1:n){
  #  TOTAL FRUIT EQUIVALENTS -------------------------------------------------
    # z_tfe[i] ~ dlnorm(log(mu_tfe[site[i],year[i]]),tau_tfe[site[i],year[i]])
    # 
    # # Prior predictive
    # y_tfe_prior[i] ~ dpois(z_tfe[i])
    # 
    # # Likelihood
    # #y_tfe[i] ~ dpois(z_tfe[i])
 # }
  
  for (i in 1:n){
    alpha_tfe[i] ~ dnorm(mu_tfe[site[i],year[i]],tau_tfe[site[i],year[i]])
    log(lambda_tfe[i]) = alpha_tfe[i]
   # log(lambda_tfe[i]) = mu_tfe[site[i],year[i]]
    
    #y_tfe[i] ~ dnegbin(kappa_tfe[site[i],year[i]]/(kappa_tfe[site[i],year[i]]+lambda_tfe[i]),kappa_tfe[site[i],year[i]])
    #y_prior[i] ~ dnegbin(kappa_tfe[site[i],year[i]]/(kappa_tfe[site[i],year[i]]+lambda_tfe[i]),kappa_tfe[site[i],year[i]])
    #y_prior[i] ~ dnegbin(kappa_tfe[site[i],year[i]]/(kappa_tfe[site[i],year[i]]+lambda_tfe[i]),kappa_tfe[site[i],year[i]])
    y_prior[i] ~ dpois(lambda_tfe[i])
      }
  

  # for (i in 1:n2){
  # #  UNDAMAGED FRUITS -------------------------------------------------
  #   z_und[i] ~ dlnorm(log(mu_und[site2[i],year2[i]]),tau_und[site2[i],year2[i]])
  #   y_und[i] ~ dpois(z_und[i])
  # }
  # 
  # for (i in 1:n2){
  #   #  DAMAGED FRUITS -------------------------------------------------
  #   z_dam[i] ~ dlnorm(log(mu_dam[site2[i],year2[i]]),tau_dam[site2[i],year2[i]])
  #   y_dam[i] ~ dpois(z_dam[i])
  # }
  # 
  # # SEEDS PER UNDAMAGED FRUITS -------------------------------------------------
  # for (i in 1:n3){
  #   z_sd[i] ~ dlnorm(log(mu_seeds[site3[i],year3[i]]),tau_seeds[site3[i],year3[i]])
  #   sdno[i] ~ dpois(z_sd[i])
  #   }
  # 
  # # SEEDS PER DAMAGED FRUITS -------------------------------------------------
  # for (i in 1:n4){
  #   z_sd_dam[i] ~ dlnorm(log(mu_dam_seeds[site4[i],year4[i]]),tau_dam_seeds[site4[i],year4[i]])
  #   sdno_dam[i] ~ dpois(z_sd_dam[i])
  #   }


# DERIVED QUANTITIES ------------------------------------------------------
# 
# 
# for(i in 1:n_site){
# 
#   # TOTAL FRUIT EQUIVALENTS : POPULATION ------------------------------------------------------
#   mu_p_tfe[i] ~ dlnorm(log(g_tfe[i]), tau0_tfe[i])
# 
#   # POPULATION AND YEAR------------------------------------------------------
#   for(j in sparse_tfe[id_tfe[i]:(id_tfe[i+1]-1)]){
#     mu_py_tfe[i,j] ~ dlnorm(log(mu_tfe[i,j]),tau_tfe[i,j])
#   }
# 
#   # UNDAMAGED FRUITS: POPULATION ------------------------------------------------------
# 
#   # for dlnorm(alpha,beta)
#   # mu = exp(alpha+((beta^2)/2))
#   # sigma2 = (exp(beta^2)-1))*exp(2*alpha+beta^2)
#   mu_p_und[i] ~ dlnorm(log(g_und[i]), tau0_und[i])
#   # beta_p_und[i] = 1/tau0_und[i]
#   # mean_p_und[i] = exp(log(g_und[i])+beta_p[i]/2)
# 
#   # DAMAGED FRUITS : POPULATION ------------------------------------------------------
#   mu_p_dam[i] ~ dlnorm(log(g_dam[i]), tau0_dam[i])
# 
#   # POPULATION AND YEAR------------------------------------------------------
# 
#   for(j in sparse_und[id_und[i]:(id_und[i+1]-1)]){
# 
#     # for dlnorm(alpha,beta)
#     # mu = exp(alpha+((beta^2)/2))
#     # sigma2 = (exp(beta^2)-1))*exp(2*alpha+beta^2)
#     mu_py_und[i,j] ~ dlnorm(log(mu_und[i,j]),tau_und[i,j])
#     # beta_py[i,j] = (1/tau_und[i,j])
#     # mean_py[i,j] = exp(log(mu_und[i,j])+beta_py[i,j]/2)
# 
#     #damaged and undamaged share index
#     mu_py_dam[i,j] ~ dlnorm(log(mu_dam[i,j]),tau_dam[i,j])
#   }
# 
#   # SEEDS PER UNDAMAGED FRUITS -------------------------------------------------
#   mu_p_seeds[i] ~ dlnorm(log(g_seeds[i]), tau0_seeds[i])
# 
#   # POPULATION AND YEAR------------------------------------------------------
#   for(j in sparse_seeds[id_seeds[i]:(id_seeds[i+1]-1)]){
#     mu_py_seeds[i,j] ~ dlnorm(log(mu_seeds[i,j]),tau_seeds[i,j])
#   }
#   # SEEDS PER DAMAGED FRUITS -------------------------------------------------
#   mu_p_dam_seeds[i] ~ dlnorm(log(g_dam_seeds[i]), tau0_dam_seeds[i])
# 
#   # POPULATION AND YEAR------------------------------------------------------
#   for(j in sparse_dam_seeds[id_dam_seeds[i]:(id_dam_seeds[i+1]-1)]){
#     mu_py_dam_seeds[i,j] ~ dlnorm(log(mu_dam_seeds[i,j]),tau_dam_seeds[i,j])
#  #}
# 
#   # calculate the fraction of seeds in a damaged versus undamaged fruit (only for 2013-2018)
#   # RATIO -------------------------------------------------------------------
#   #for(j in sparse_dam_seeds[id_dam_seeds[i]:(id_dam_seeds[i+1]-1)]){
#       ratio[i,j] = mu_py_dam_seeds[i,j]/mu_py_seeds[i,j]
#       # composite tfe
#       # only DLW doesn't have data on seeds from damaged fruits but has damaged fruits
#       mu_py_tfe_comp[i,j] = mu_py_und[i,j] + mu_py_dam[i,j]*ratio[i,j]
#     }
# 
#   }

}
