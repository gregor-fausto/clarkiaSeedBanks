model {
 # Priors ------------------------------------------------------------------
  
  # * Global ------------------------------------------------------------------
      for (i in 1:n_site) {
    # * Survival parameters ----------------------------------------------------------------
    mu0_s[i] ~  dnorm(0, 0.001)
    sigma0_s[i] ~ dunif(0, 1.5)
    tau0_s[i] <- 1 / (sigma0_s[i] * sigma0_s[i])

    # * Fecundity parameters ----------------------------------------------------------------
    mu0_f[i] ~  dnorm(0, 0.001)
    sigma0_f[i] ~ dunif(0,1.5)
    tau0_f[i] <- 1/(sigma0_f[i]*sigma0_f[i])        

    # ** Group -------------------------------------------------------------------
    for (j in 1:n_year) {
      # ** Seed survival parameters ----------------------------------------------------------------
      mu_s[i, j] ~ dnorm(mu0_s[i], tau0_s[i])
      sigma_s[i, j] ~ dunif(0, 1.5)
      tau_s[i, j] <- 1 / (sigma_s[i, j] * sigma_s[i, j])
      
      # ** Fecundity parameters ----------------------------------------------------------------
      gamma_f[i,j] ~ dnorm(mu0_f[i], tau0_f[i])
      r_f[i,j] ~ dgamma(.001,.001)     
    }
  }

  # Likelihoods --------------------------------------------------------------
 
   # * Survival --------------------------------------------------------------
  
  for (i in 1:n) {
    # ** Inverse logit ----------------------------------------------------------------
    alpha_s[i] ~ dnorm(mu_s[site[i], year[i]], tau_s[site[i], year[i]])
    logit(theta_s[i]) <- alpha_s[i]
    
    # ** Binomial likelihood ----------------------------------------------------------------
    fruitplNumber[i] ~ dbinom(theta_s[i], seedlingNumber[i])
  }
  
# * Fecundity --------------------------------------------------------------

    for (i in 1:n){
    # ** Transform ----------------------------------------------------------------
    lambda_f[i] = exp(gamma_f[site[i],year[i]])
    # ** Negative binomial likelihood ----------------------------------------------------------------
    sdno[i] ~ dnegbin(r_f[site[i],year[i]]/(r_f[site[i],year[i]]+lambda_f[i]),r_f[site[i],year[i]])
  }  

  # Derived quantities --------------------------------------------------------------

  # * Suvival --------------------------------------------------------------
  
  for (i in 1:n_site) {
    p.i0_s[i] ~ dnorm(mu0_s[i], tau0_s[i])
    logit(p0_s[i]) <- p.i0_s[i]
    
    for (j in 1:n_year) {
      p.i_s[i, j] ~ dnorm(mu_s[i, j], tau_s[i, j])
      logit(p_s[i, j]) <- p.i_s[i, j]
      
    }
  }
  
  # * Fecundity --------------------------------------------------------------
  
  for(i in 1:n_site){
    
    p.i0_f[i] ~ dnorm(mu0_f[i],tau0_f[i])
    p0_f[i] <- exp(p.i0_f[i])
    
    for (j in 1:n_year) {
    p_f[i,j] <- exp(gamma_f[i,j])
      
    }
    
  }
}
