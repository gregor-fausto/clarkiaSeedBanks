
model {
  
  # PRIORS ------------------------------------------------
  
  for(i in 1:n_site){
    
    # TOTAL FRUIT EQUIVALENTS -----------------------------------------
    
    nu_tfe[i] ~ dgamma(1, 1)
    g_tfe[i] = exp(nu_tfe[i])
    
    # page 76 https://nwfsc-timeseries.github.io/atsa-2017/Labs/Week%203%20intro%20to%20jags/intro-to-jags.pdf
    sigma0_tfe[i] ~ dnorm(0, 1) T(0,)
    tau0_tfe[i] <- 1/(sigma0_tfe[i]*sigma0_tfe[i])
    
    for(j in 1:n_year){
      
      mu_tfe[i,j] ~ dlnorm(log(g_tfe[i]), tau0_tfe[i])
      sigma_tfe[i,j] ~ dnorm(0, 1) T(0,)
      tau_tfe[i,j] <- 1/(sigma_tfe[i,j]*sigma_tfe[i,j])

    }
  }
  
  
  for(i in 1:n_site2){
    
    # UNDAMAGED/DAMAGED FRUITS ------------------------------------------------
    
    ## UNDAMAGED
    nu_und[i] ~ dgamma(1, 1)
    g_und[i] = exp(nu_und[i])
    
    # page 76 https://nwfsc-timeseries.github.io/atsa-2017/Labs/Week%203%20intro%20to%20jags/intro-to-jags.pdf
    sigma0_und[i] ~ dnorm(0, 1) T(0,)
    tau0_und[i] <- 1/(sigma0_und[i]*sigma0_und[i])
    
    ## DAMAGED
    nu_dam[i] ~ dgamma(1, 1)
    g_dam[i] = exp(nu_dam[i])
    
    # page 76 https://nwfsc-timeseries.github.io/atsa-2017/Labs/Week%203%20intro%20to%20jags/intro-to-jags.pdf
    sigma0_dam[i] ~ dnorm(0, 1) T(0,)
    tau0_dam[i] <- 1/(sigma0_dam[i]*sigma0_dam[i])
    
    for(j in 1:n_year2){
      
      mu_und[i,j] ~ dlnorm(log(g_und[i]), tau0_und[i])
      sigma_und[i,j] ~ dnorm(0, 1) T(0,)
      tau_und[i,j] <- 1/(sigma_und[i,j]*sigma_und[i,j])
      
      mu_dam[i,j] ~ dlnorm(log(g_dam[i]), tau0_dam[i])
      sigma_dam[i,j] ~ dnorm(0, 1) T(0,)
      tau_dam[i,j] <- 1/(sigma_dam[i,j]*sigma_dam[i,j])
      
    }
  }
  
  for(i in 1:n_site){
    
    # SEEDS PER UNDAMAGED FRUITS  -----------------------------------------
    
    nu_seeds[i] ~ dgamma(1, 1)
    g_seeds[i] = exp(nu_seeds[i])
    
    # page 76 https://nwfsc-timeseries.github.io/atsa-2017/Labs/Week%203%20intro%20to%20jags/intro-to-jags.pdf
    sigma0_seeds[i] ~ dnorm(0, 1) T(0,)
    tau0_seeds[i] <- 1/(sigma0_seeds[i]*sigma0_seeds[i])
    
    for(j in 1:n_year3){
      
      mu_seeds[i,j] ~ dlnorm(log(g_seeds[i]), tau0_seeds[i])
      sigma_seeds[i,j] ~ dnorm(0, 1) T(0,)
      tau_seeds[i,j] <- 1/(sigma_seeds[i,j]*sigma_seeds[i,j])
      
    }
  }
  
  for(i in 1:n_site){
    
    # SEEDS PER DAMAGED FRUITS  -----------------------------------------
    
    nu_dam_seeds[i] ~ dgamma(1, 1)
    g_dam_seeds[i] = exp(nu_dam_seeds[i])
    
    # page 76 https://nwfsc-timeseries.github.io/atsa-2017/Labs/Week%203%20intro%20to%20jags/intro-to-jags.pdf
    sigma0_dam_seeds[i] ~ dnorm(0, 1) T(0,)
    tau0_dam_seeds[i] <- 1/(sigma0_dam_seeds[i]*sigma0_dam_seeds[i])
    
    for(j in 1:n_year4){
      
      mu_dam_seeds[i,j] ~ dlnorm(log(g_dam_seeds[i]), tau0_dam_seeds[i])
      sigma_dam_seeds[i,j] ~ dnorm(0, 1) T(0,)
      tau_dam_seeds[i,j] <- 1/(sigma_dam_seeds[i,j]*sigma_dam_seeds[i,j])      
    }
  }
  
  # LIKELIHOODS -------------------------------------------------------------
  
  for (i in 1:n){
    #  TOTAL FRUIT EQUIVALENTS -------------------------------------------------
    z_tfe[i] ~ dlnorm(log(mu_tfe[site[i],year[i]]),tau_tfe[site[i],year[i]])
    y_tfe[i] ~ dpois(z_tfe[i])
    
    y_tfe.sim[i] ~ dpois(z_tfe[i])
    chi2.tfe.obs[i] <- pow((y_tfe[i] - z_tfe[i]),2) / (z_tfe[i])         # obs.
    chi2.tfe.sim[i] <- pow((y_tfe.sim[i]- z_tfe[i]),2) / (z_tfe[i])         # obs.
    
  }
  
  for (i in 1:n2){
    #  UNDAMAGED FRUITS -------------------------------------------------
    z_und[i] ~ dlnorm(log(mu_und[site2[i],year2[i]]),tau_und[site2[i],year2[i]])
    y_und[i] ~ dpois(z_und[i])
    
    y_und.sim[i] ~ dpois(z_und[i])
    chi2.und.obs[i] <- pow((y_und[i] - z_und[i]),2) / (z_und[i])         # obs.
    chi2.und.sim[i] <- pow((y_und.sim[i]- z_und[i]),2) / (z_und[i])         # obs.
  }
  
  for (i in 1:n2){
    #  DAMAGED FRUITS -------------------------------------------------
    z_dam[i] ~ dlnorm(log(mu_dam[site2[i],year2[i]]),tau_dam[site2[i],year2[i]])
    y_dam[i] ~ dpois(z_dam[i])
    
    y_dam.sim[i] ~ dpois(z_dam[i])
    chi2.dam.obs[i] <- pow((y_dam[i] - z_dam[i]),2) / (z_dam[i])         # obs.
    chi2.dam.sim[i] <- pow((y_dam.sim[i]- z_dam[i]),2) / (z_dam[i])         # obs.
  }
  
  # SEEDS PER UNDAMAGED FRUITS -------------------------------------------------
  for (i in 1:n3){
    z_sd[i] ~ dlnorm(log(mu_seeds[site3[i],year3[i]]),tau_seeds[site3[i],year3[i]])
    sdno[i] ~ dpois(z_sd[i])
    
    y_sd.sim[i] ~ dpois(z_sd[i])
    chi2.sd.obs[i] <- pow((sdno[i] - z_sd[i]),2) / (z_sd[i])         # obs.
    chi2.sd.sim[i] <- pow((y_sd.sim[i]- z_sd[i]),2) / (z_sd[i])         # obs.
  }
  
  # SEEDS PER DAMAGED FRUITS -------------------------------------------------
  for (i in 1:n4){
    z_sd_dam[i] ~ dlnorm(log(mu_dam_seeds[site4[i],year4[i]]),tau_dam_seeds[site4[i],year4[i]])
    sdno_dam[i] ~ dpois(z_sd_dam[i])
    
    y_sd_dam.sim[i] ~ dpois(z_sd_dam[i])
    chi2.sd_dam.obs[i] <- pow((sdno_dam[i] - z_sd_dam[i]),2) / (z_sd_dam[i])         # obs.
    chi2.sd_dam.sim[i] <- pow((y_sd_dam.sim[i]- z_sd_dam[i]),2) / (z_sd_dam[i])         # obs.
  }
  
  
  # DERIVED QUANTITIES ------------------------------------------------------
  
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
  #     #}
  #     
  #     # calculate the fraction of seeds in a damaged versus undamaged fruit (only for 2013-2018)
  #     # RATIO -------------------------------------------------------------------
  #     #for(j in sparse_dam_seeds[id_dam_seeds[i]:(id_dam_seeds[i+1]-1)]){
  #     ratio[i,j] = mu_py_dam_seeds[i,j]/mu_py_seeds[i,j]
  #     # composite tfe
  #     # only DLW doesn't have data on seeds from damaged fruits but has damaged fruits
  #     mu_py_tfe_comp[i,j] = mu_py_und[i,j] + mu_py_dam[i,j]*ratio[i,j]
  #   }
  #   
  # }
  
}
