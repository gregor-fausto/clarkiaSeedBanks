
model {
  
  # PRIORS ------------------------------------------------    
  
  for(i in 1:n_site3){
    
    # SEEDS PER UNDAMAGED FRUITS  -----------------------------------------
    
    for(j in 1:n_year3){
      
      mu.log_seeds[i,j] ~ dgamma(1, 1)
      sigma_seeds[i,j] ~ dnorm(0, 2) T(0,)
      tau_seeds[i,j] <- 1/(sigma_seeds[i,j]*sigma_seeds[i,j])
      
    }
  }
  
  for(i in 1:n_site3){
    
    # SEEDS PER DAMAGED FRUITS  -----------------------------------------
    
    for(j in 1:n_year4){
      
      mu.log_dam_seeds[i,j] ~ dgamma(1, 1)
      sigma_dam_seeds[i,j] ~ dnorm(0, 2) T(0,)
      tau_dam_seeds[i,j] <- 1/(sigma_dam_seeds[i,j]*sigma_dam_seeds[i,j])      
    }
  }
  
  # LIKELIHOODS -------------------------------------------------------------
  
  # SEEDS PER UNDAMAGED FRUITS -------------------------------------------------
  for (i in 1:n3){
    z_sd[i] ~ dlnorm(mu.log_seeds[site3[i],year3[i]],tau_seeds[site3[i],year3[i]])
    sdno[i] ~ dpois(z_sd[i])

    # Posterior predictive and chi-squared
    y_sd.sim[i] ~ dpois(z_sd[i])
    chi2.sd.obs[i] <- pow((sdno[i] - z_sd[i]),2) / (z_sd[i])
    chi2.sd.sim[i] <- pow((y_sd.sim[i]- z_sd[i]),2) / (z_sd[i])
  }
  
  # # SEEDS PER DAMAGED FRUITS -------------------------------------------------
  for (i in 1:n4){
    z_sd_dam[i] ~ dlnorm(mu.log_dam_seeds[site4[i],year4[i]],tau_dam_seeds[site4[i],year4[i]])
    sdno_dam[i] ~ dpois(z_sd_dam[i])
    
    # Posterior predictive and chi-squared
    y_sd_dam.sim[i] ~ dpois(z_sd_dam[i])
    chi2.sd_dam.obs[i] <- pow((sdno_dam[i] - z_sd_dam[i]),2) / (z_sd_dam[i])
    chi2.sd_dam.sim[i] <- pow((y_sd_dam.sim[i]- z_sd_dam[i]),2) / (z_sd_dam[i])
  }
  
  
}
