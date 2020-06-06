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
