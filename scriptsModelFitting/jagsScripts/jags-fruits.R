
model {
  
  # PRIORS ------------------------------------------------
  
  for(i in 1:n_site){
    
    # TOTAL FRUIT EQUIVALENTS -----------------------------------------
    
    nu_tfe[i] ~ dgamma(1, 1)
    #g_tfe[i] = exp(nu_tfe[i])
    
    # page 76 https://nwfsc-timeseries.github.io/atsa-2017/Labs/Week%203%20intro%20to%20jags/intro-to-jags.pdf
    sigma0_tfe[i] ~ dnorm(0, 2) T(0,)
    tau0_tfe[i] <- 1/(sigma0_tfe[i]*sigma0_tfe[i])
    
    for(j in 1:n_year){
      
      mu_tfe[i,j] ~ dlnorm(nu_tfe[i], tau0_tfe[i])
      mu.log_tfe[i,j] <- log(mu_tfe[i,j])
      sigma_tfe[i,j] ~ dnorm(0, 2) T(0,)
      tau_tfe[i,j] <- 1/(sigma_tfe[i,j]*sigma_tfe[i,j])
      
    }
  }
  
  
  for(i in 1:n_site2){
    
    # UNDAMAGED/DAMAGED FRUITS ------------------------------------------------
    
    ## UNDAMAGED
    nu_und[i] ~ dgamma(1, 1)
    #g_und[i] = exp(nu_und[i])
    
    # page 76 https://nwfsc-timeseries.github.io/atsa-2017/Labs/Week%203%20intro%20to%20jags/intro-to-jags.pdf
    sigma0_und[i] ~ dnorm(0, 2) T(0,)
    tau0_und[i] <- 1/(sigma0_und[i]*sigma0_und[i])
    
    ## DAMAGED
    nu_dam[i] ~ dgamma(1, 1)
    #g_dam[i] = exp(nu_dam[i])
    
    # page 76 https://nwfsc-timeseries.github.io/atsa-2017/Labs/Week%203%20intro%20to%20jags/intro-to-jags.pdf
    sigma0_dam[i] ~ dnorm(0, 2) T(0,)
    tau0_dam[i] <- 1/(sigma0_dam[i]*sigma0_dam[i])
    
    for(j in 1:n_year2){
      
      mu_und[i,j] ~ dlnorm(nu_und[i], tau0_und[i])
      mu.log_und[i,j] <- log(mu_und[i,j])
      sigma_und[i,j] ~ dnorm(0, 2) T(0,)
      tau_und[i,j] <- 1/(sigma_und[i,j]*sigma_und[i,j])
      
      mu_dam[i,j] ~ dlnorm(nu_dam[i], tau0_dam[i])
      mu.log_dam[i,j] <- log(mu_dam[i,j])
      sigma_dam[i,j] ~ dnorm(0, 2) T(0,)
      tau_dam[i,j] <- 1/(sigma_dam[i,j]*sigma_dam[i,j])
      
    }
  }
  
  # LIKELIHOODS -------------------------------------------------------------
  
  for (i in 1:n){
    
    #  TOTAL FRUIT EQUIVALENTS -------------------------------------------------
    z_tfe[i] ~ dlnorm(mu.log_tfe[site[i],year[i]],tau_tfe[site[i],year[i]])
    y_tfe[i] ~ dpois(z_tfe[i])
    
    # Posterior predictive and chi-squared
    y_tfe.sim[i] ~ dpois(z_tfe[i])
    chi2.tfe.obs[i] <- pow((y_tfe[i] - z_tfe[i]),2) / (z_tfe[i])
    chi2.tfe.sim[i] <- pow((y_tfe.sim[i]- z_tfe[i]),2) / (z_tfe[i])
    
  }
  
  for (i in 1:n2){
    #  UNDAMAGED FRUITS -------------------------------------------------
    z_und[i] ~ dlnorm(mu.log_und[site2[i],year2[i]],tau_und[site2[i],year2[i]])
    y_und[i] ~ dpois(z_und[i])
    
    # Posterior predictive and chi-squared
    y_und.sim[i] ~ dpois(z_und[i])
    chi2.und.obs[i] <- pow((y_und[i] - z_und[i]),2) / (z_und[i])
    chi2.und.sim[i] <- pow((y_und.sim[i]- z_und[i]),2) / (z_und[i])
  }
  
  for (i in 1:n2){
    #  DAMAGED FRUITS -------------------------------------------------
    z_dam[i] ~ dlnorm(mu.log_dam[site2[i],year2[i]],tau_dam[site2[i],year2[i]])
    y_dam[i] ~ dpois(z_dam[i])
    
    # Posterior predictive and chi-squared
    y_dam.sim[i] ~ dpois(z_dam[i])
    chi2.dam.obs[i] <- pow((y_dam[i] - z_dam[i]),2) / (z_dam[i])
    chi2.dam.sim[i] <- pow((y_dam.sim[i]- z_dam[i]),2) / (z_dam[i])
  }
  
  
}
