model{
  
  # LIKELIHOOD
  for(i in 1:n){
    seedlingJan[i] ~ dbinom(g[i], totalJan[i])
    logit(g[i]) <- alpha[siteGermination[i]] + beta1[siteGermination[i]]*t.std[i] + beta2[siteGermination[i]]*p.std[i]
  }
 
  # model for group intercept and slope
  for(j in 1:n_siteGermination){
    # get slopes for each site
    alpha[j] <- B[j, 1]  # group level intercept
    beta1[j]  <- B[j, 2]  # group level slope 1
    beta2[j]  <- B[j, 3]  # group level slope 2
   
    # generate correlated intercept and slopes
    B[j, 1:3] ~ dmnorm(B.hat[j, 1:3], Tau.B)  
    B.hat[j, 1] <- mu.alpha  # required by JAGS syntax
    B.hat[j, 2] <- mu.beta1   # required by JAGS syntax
    B.hat[j, 3] <- mu.beta2   # required by JAGS syntax
  }
  
  # priors
  mu.alpha ~ dnorm(0, .0001) 
  mu.beta1 ~ dnorm(0, .0001) 
  mu.beta2 ~ dnorm(0, .0001) 
  
  # inverse of covariance matrix 
  Tau.B[1:3, 1:3] <- inverse(Sigma.B[1:3, 1:3])
  
  # diagonal elements of covariance matrix
  Sigma.B[1, 1] <- sigma.alpha^2
  sigma.alpha ~ dunif(0,200)
  Sigma.B[2, 2] <- sigma.beta1^2
  sigma.beta1 ~ dunif(0,200)
  Sigma.B[3, 3] <- sigma.beta2^2
  sigma.beta2 ~ dunif(0,200)
  
  # covariance is correlation coef. x product of standard deviations
  Sigma.B[1, 2] <- rho.a.b1 * sigma.alpha * sigma.beta1 
  Sigma.B[2, 1] <- Sigma.B[1,2]
  Sigma.B[1, 3] <- rho.a.b2 * sigma.alpha * sigma.beta2 
  Sigma.B[3, 1] <- Sigma.B[1,3]
  Sigma.B[2, 3] <- rho.b1.b2 * sigma.beta1 * sigma.beta2 
  Sigma.B[3, 2] <- Sigma.B[2,3]
  rho.a.b1 ~ dunif(-1, 1)
  rho.a.b2 ~ dunif(-1, 1)
  rho.b1.b2 ~ dunif(-1, 1)
  
}