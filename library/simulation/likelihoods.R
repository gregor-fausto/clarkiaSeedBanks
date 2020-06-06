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
