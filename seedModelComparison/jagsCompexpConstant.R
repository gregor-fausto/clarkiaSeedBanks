model {
  
  # seed bags ---------------------------------------------------------------
  
  a ~ dgamma(0.001,0.001)
  b ~ dgamma(0.001,0.001)
  
  beta ~ dgamma(0.001,0.001) 
  beta2 ~ dgamma(0.001,0.001) 
  beta3 ~ dgamma(0.001,0.001) 
  beta4 ~ dgamma(0.001,0.001) 
  
  g1 ~ dbeta(1,1)
  
  # # Inverse gamma is not a built in distribution
  # # so place a gamma distribution on the precision
  #tau ~ dgamma(0.001,0.001)
  #sigma = 1/tau
  # sigma ~ dunif(0.0001,.5)
  #sigma2 ~ dunif(0.0001,.1)
  # mu ~ dunif(0,1)
  
  # alpha ~ dbeta(1,1)
  # beta ~ dbeta(1,1)
  
  # likelihood (seed bags) --------------------------------------------------------------
  
  for(i in 1:n1){
    
    mu[i] <- 1/((1+b*months[i])^a)
    alpha[i] <- (mu[i]*beta)/(1-mu[i])
    theta[i] ~ dbeta(alpha[i],beta)
    
    mu2[i] <- (1-g1)*(1/((1+b*months_2[i])^a))
    alpha2[i] <- (mu2[i]*beta2)/(1-mu2[i])
    theta2[i] ~ dbeta(alpha2[i],beta2)
    
    # likelihood
    totalJan[i] ~ dbinom(theta[i], seedStart[i])
    seedlingJan[i] ~ dbinom(g1, totalJan[i])
    intactOct_2[i] ~dbinom(theta2[i], seedStart_2[i])
  }
  
  for(i in 1:n2){
    
    mu3[i] <- (1-g1)*(1/((1+b*months_3[i])^a))
    alpha3[i] <- (mu3[i]*beta3)/(1-mu3[i])
    theta3[i] ~ dbeta(alpha3[i],beta3)
    
    mu4[i] <- (1-g1)*(1-g1)*(1/((1+b*months_4[i])^a))
    alpha4[i] <- (mu4[i]*beta4)/(1-mu4[i])
    theta4[i] ~ dbeta(alpha4[i],beta4)
    
    # likelihood
    totalJan_3[i] ~ dbinom(theta3[i], seedStart_3[i])
    seedlingJan_3[i] ~ dbinom(g1, totalJan_3[i])
    intactOct_4[i] ~dbinom(theta4[i], seedStart_4[i])
  }
  
  
  # derived quantities --------------------------------------------------------------
  
  
  
}
