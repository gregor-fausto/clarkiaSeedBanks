model{
  
  # LIKELIHOOD
  for(i in 1:n){
    seedlingJan[i] ~ dbinom(g[i], totalJan[i])
    logit(g[i]) <- alpha[siteGermination[i]] + beta[siteGermination[i]]*t.std[i]
  }
  
  # model for group intercept and slope
  for(j in 1:n_siteGermination){
    # get intercept and slopes
    alpha[j] <- B[j, 1]  # group level intercept
    beta[j]  <- B[j, 2]  # group level slope
    # draw parameters from MVN
    B[j, 1:2] ~ dmnorm(B.hat[j, ], Tau.B[,])
    # means come from vectors
    B.hat[j, 1] <- mu.alpha  # required by JAGS syntax
    B.hat[j, 2] <- mu.beta   # required by JAGS syntax
  }
  
  # priors on means
  mu.alpha ~ dnorm(0, 1)
  mu.beta ~ dnorm(0, 1)
  
  # get correlations
  rho.a.b1 <- rho.B[1,2]
 #rho.a.b2 <- rho.B[1,3]
 # rho.b1.b2 <- rho.B[2,3]
  
  for (k in 1:2){
    for (k.prime in 1:2){
      rho.B[k,k.prime] <- Sigma.B[k,k.prime]/sqrt(Sigma.B[k,k]*Sigma.B[k.prime,k.prime])
    }
    # get standard deviations
    sigma.B[k] <- sqrt(Sigma.B[k,k])
  }
  
  #set K = 2
  ### Model variance-covariance with Wishert disttribution
  Sigma.B[1:2,1:2] <- inverse(Tau.B[,])
  Tau.B[1:2,1:2] ~ dscaled.wishart(s, df)
  df <- 2
  
}
