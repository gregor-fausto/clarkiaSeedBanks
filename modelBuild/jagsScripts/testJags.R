model { 
## Priors
##hyperprior for intercept alpha

  for(i in 1:nsites){
    sigma.site[i]~dunif(0,100)
    tau.site[i]<-1/(sigma.site[i]*sigma.site[i])
   # mu.alpha[i]~dnorm(0,0.001)   
  }
  
  for(i in 1:nsiteyears){
    sigma.year[i]~dunif(0,100)
    tau.year[i]<-1/(sigma.year[i]*sigma.year[i])
  }

  # prior
  for(j in 1:nsites){
    #site intercepts
    alpha.i[j] ~ dnorm(0, tau.site[site[j]])
    # alpha.s[j] ~ dnorm(mu.s[site[j]], tau.s[site[j]])
  }
  for(j in 1:nsiteyears){
      beta.i[j] ~ dnorm(alpha.i[site_year[j]],tau.year[j])
      #beta.s[j] ~ dnorm(mu.b.s[siteyear[j]], tau.b.s[siteyear[j]])
  }

## Likelihood
for(i in 1:N_burial){
    pi[i] <- ilogit(beta.i[siteyear[i]])#+beta.i[siteyear[i]])
  y_total[i] ~ dbin(pi[i], n_buried[i])

}
}


