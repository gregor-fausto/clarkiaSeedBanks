

model{

  # prior
  
  lambda ~ dgamma(0.001,0.001)
 
  # data model
  y[i] ~ dpois(lambda)  
  }# end of data model

} #end of model


