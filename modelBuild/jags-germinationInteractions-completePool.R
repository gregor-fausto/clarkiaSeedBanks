model {
  
  # Likelihood --------------------------------------------------------------
  
  # for (i in 1:length(temp_pred_std)) {
  #   for(j in 1:length(precip_pred_std)){
  #     logit(g_pred[i,j]) <- alpha + beta[1]*temp_pred_std[i] + beta[2]*precip_pred_std[j] + beta[3]*temp_pred_std[i]*precip_pred_std[j]
  #   }
  #  
  # }
  
  # seed germination --------------------------------------------------------------
  for(i in 1:n){
    
    # LIKELIHOOD
    seedlingJan[i] ~ dbinom(g[i], totalJan[i])
    
    logit(g[i]) <- alpha + beta[1]*clim.std[i] #+ beta[2]*p.std[i] + beta[3]*t.std[i]*p.std[i]
    #alpha_g[i] ~ dnorm(mu_g[siteGermination[i]],tau_g[siteGermination[i]])
    
  }
  
 
  # for(i in 1:n_yearGermination){
  #   mu0_y[i] ~ dnorm(0, 1)
  # }
  
 

  # for(i in 1:n_siteGermination){
  #   
  #   # i indexes site; j indexes year; k indexes germination variable
  #   mu_g[i] ~ dnorm(mu0_g, tau0_g)
  #   
  #   # Weakly informative prior (Rosenbaum et al. 2019)
  #   sigma_g[i] ~ dnorm(0,1) T(0,)
  #   tau_g[i] <- 1/(sigma_g[i]*sigma_g[i])
  # }
  
  # Priors --------------------------------------------------------------
  for (i in 1){
    beta[i] ~ dnorm(0, 1)  
  }
  
  alpha ~ dnorm(0,1)
  
  # Weakly informative prior (Rosenbaum et al. 2019)
  # slightly more informative prior because few groups https://statmodeling.stat.columbia.edu/2015/12/08/hierarchical-modeling-when-you-have-only-2-groups-i-still-think-its-a-good-idea-you-just-need-an-informative-prior-on-the-group-level-variation/
  #sigma0_g ~ dnorm(0,1) T(0,)
  #tau0_g <- 1/(sigma0_g*sigma0_g)
  
  # seed germination ---------------------------------------------------------------
  # population*year parameter, hierarchical logit
  #mu0_g ~ dnorm(0, 1)
  
 
  
  
  # for(i in 1:n_siteGermination){
  #   
  #   alpha_pred[i] ~ dnorm(mu_g[siteGermination[i]],tau_g[siteGermination[i]])
  #   
  #   for (j in 1:length(temp_pred_std)) {
  #     # gets the germination probability across temps at mean precipitation (0, centered) 
  #     logit(psiPredict[i,j]) <- alpha_pred[i] + beta[1] * temp_pred_std[j]
  #   }
  #   
  #   for (j in 1:length(precip_pred_std)) {
  #     # gets the germination probability across temps at mean precipitation (0, centered) 
  #     logit(phiPredict[i,j]) <- alpha_pred[i] + beta[2] * precip_pred_std[j]
  #   }
  #   
  # }
  
   
}



