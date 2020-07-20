
model {

  # for(i in 1:n_site2){
   for(i in 1:n_site2){
  
    # UNDAMAGED FRUITS

    nu_und[i] ~ dgamma(1, 1)
    g_und[i] = exp(nu_und[i])
    
    # page 76 https://nwfsc-timeseries.github.io/atsa-2017/Labs/Week%203%20intro%20to%20jags/intro-to-jags.pdf
    tau0_und[i] ~ dgamma(0.001,0.001)  # Inverse gamma is not a built in distribution

    for(j in 1:n_year2){

      # theta 1
      mu_und[i,j] ~ dlnorm(log(g_und[i]), tau0_und[i])
      
      tau_und[i,j] ~ dgamma(0.001,0.001)  # Inverse gamma is not a built in distribution

    }
  }

  ## LIKELIHOOD
  for (i in 1:n2){
    # for lognormal(a,b); first parameter is mu;
    # median of lognormal is exp(mu)
    # log(median) = mu
    z_und[i] ~ dlnorm(log(mu_und[site2[i],year2[i]]),tau_und[site2[i],year2[i]])
    y_und[i] ~ dpois(z_und[i])
   # y_und_sim[i] ~ dpois(lambda_und[i])
    }

### UNDAMAGED PLANTS
for(i in 1:n_site2){
  # for dlnorm(alpha,beta)
  # mu = exp(alpha+((beta^2)/2))
  # sigma2 = (exp(beta^2)-1))*exp(2*alpha+beta^2)
  mu_p[i] ~ dlnorm(log(g_und[i]), tau0_und[i])

 # beta_p[i] = 1/tau0_und[i]
 # mean_p[i] = exp(log(g_und[i])+beta_p[i]/2)
  
  # for (k in 1:length(item)) {
  #   sparse[id[item[k]]:(id[item[k]+1] - 1)]
  # }
  
  for(j in sparse[id[i]:(id[i+1]-1)]){

    # for dlnorm(alpha,beta)
    # mu = exp(alpha+((beta^2)/2))
    # sigma2 = (exp(beta^2)-1))*exp(2*alpha+beta^2)
    mu_py[i,j] ~ dlnorm(log(mu_und[i,j]),tau_und[i,j])

    beta_py[i,j] = (1/tau_und[i,j])
    mean_py[i,j] = exp(log(mu_und[i,j])+beta_py[i,j]/2)

  }
}
  
  
  

 }
