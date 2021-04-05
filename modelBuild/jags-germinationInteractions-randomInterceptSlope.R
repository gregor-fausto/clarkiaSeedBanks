model{
  
  # LIKELIHOOD
  for(i in 1:n){
    seedlingJan[i] ~ dbinom(g[i], totalJan[i])
    logit(g[i]) <- alpha[siteGermination[i]] + beta[siteGermination[i]]*t.std[i]
  }
 
  # model for group intercept and slope
  for(j in 1:n_siteGermination){
    # get slopes
    alpha[j] <- B[j, 1]  # group level intercept
    beta[j]  <- B[j, 2]  # group level slope
    # generate slopes from multivariatenormal
    B[j, 1:2] ~ dmnorm(B.hat[j,1:2], Tau.B)  
    # means of multivariate normal
    B.hat[j, 1] <- mu.alpha  # required by JAGS syntax
    B.hat[j, 2] <- mu.beta   # required by JAGS syntax
  }
  
  # priors
  mu.alpha ~ dnorm(0, 1) 
  mu.beta ~ dnorm(0, 1)
  
  # inverse of covariance matrix 
  Tau.B[1:2, 1:2] <- inverse(Sigma.B[1:2, 1:2])
  
  # diagonal elements of covariance matrix
  Sigma.B[1, 1] <- sigma.alpha^2
  sigma.alpha ~ dunif(0,200)
  Sigma.B[2, 2] <- sigma.beta^2
  sigma.beta ~ dunif(0,200)
  
  # covariance is correlation coef. x product of standard deviations
  Sigma.B[1, 2] <- rho * sigma.alpha * sigma.beta 
  Sigma.B[2, 1] <- Sigma.B[1,2]
  rho ~ dunif(-1, 1)
  
}