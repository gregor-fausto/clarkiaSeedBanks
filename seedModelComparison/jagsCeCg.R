model {
  
  # seed bags ---------------------------------------------------------------
  
  a ~ dgamma(0.001,0.001)
  b ~ dgamma(0.001,0.001)
  
  for(i in 1:6){
    beta[i] ~ dgamma(0.001,0.001)
  }
  
  g ~ dbeta(1,1)
  
  # likelihood (seed bags) --------------------------------------------------------------
  
  for(i in 1:n1){
    
    mu[i] <- 1/((1+b*months1[i])^a)
    alpha[i] <- (mu[i]*beta[1])/(1-mu[i])
    theta[i] ~ dbeta(alpha[i]+.01,beta[1]+.01)
    
    mu2[i] <- (1-g)*(1/((1+b*months2[i])^a))
    alpha2[i] <- (mu2[i]*beta[2])/(1-mu2[i])
    theta2[i] ~ dbeta(alpha2[i]+.01,beta[2]+.01)
    
    # likelihood
    totalJan1[i] ~ dbinom(theta[i], seedStart1[i])
    seedlingJan1[i] ~ dbinom(g, totalJan1[i])
    intactOct2[i] ~ dbinom(theta2[i], seedStart2[i])
    
    # log-likelihood
    logLik[i] <- logdensity.bin(totalJan1[i],theta[i],seedStart1[i])
    
  }
  
  for(i in 1:n2){
    
    mu3[i] <- (1-g)*(1/((1+b*months3[i])^a))
    alpha3[i] <- (mu3[i]*beta[3])/(1-mu3[i])
    theta3[i] ~ dbeta(alpha3[i]+.01,beta[3]+.01)
    
    mu4[i] <- (1-g)*(1-g)*(1/((1+b*months4[i])^a))
    alpha4[i] <- (mu4[i]*beta[4])/(1-mu4[i])
    theta4[i] ~ dbeta(alpha4[i]+.01,beta[4]+.01)
    
    # likelihood
    totalJan3[i] ~ dbinom(theta3[i], seedStart3[i])
    seedlingJan3[i] ~ dbinom(g, totalJan3[i])
    intactOct4[i] ~dbinom(theta4[i], seedStart4[i])
  }
  
  for(i in 1:n3){
    
    mu5[i] <- (1-g)*(1-g)*(1/((1+b*months5[i])^a))
    alpha5[i] <- (mu5[i]*beta[5])/(1-mu5[i])
    theta5[i] ~ dbeta(alpha5[i]+.01,beta[5]+.01)
    
    mu6[i] <- (1-g)*(1-g)*(1-g)*(1/((1+b*months6[i])^a))
    alpha6[i] <- (mu6[i]*beta[6])/(1-mu6[i])
    theta6[i] ~ dbeta(alpha6[i]+.01,beta[6]+.01)
    
    # likelihood
    totalJan5[i] ~ dbinom(theta5[i], seedStart5[i])
    seedlingJan5[i] ~ dbinom(g, totalJan5[i])
    intactOct6[i] ~dbinom(theta6[i], seedStart6[i])
  }
  
  
  # derived quantities --------------------------------------------------------------
  
  
  
}
