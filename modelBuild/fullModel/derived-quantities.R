## write code to work with derived quantities for seed bags

samples = readRDS("/Users/Gregor/Dropbox/dataLibrary/workflow/samples/belowgroundSamples.RDS")
MCMCsummary(samples.rjags,"p_g")
hist(MCMCchains(samples.rjags,"p_g")[,2],breaks=50)

mu0_g <- MCMCchains(samples,params="mu0_g")
mu0_v <- MCMCchains(samples,params="mu0_v")
sigma0_g <- MCMCchains(samples,params="sigma0_g")
sigma0_v <- MCMCchains(samples,params="sigma0_v")

mu0_g1 = mu0_g[,grep("mu0_g[",colnames(mu0_g),fixed=TRUE)]
#mu0_g2 = mu0_g[,grep("mu0_g[2",colnames(mu0_g),fixed=TRUE)]
#mu0_g3 = mu0_g[,grep("mu0_g[3",colnames(mu0_g),fixed=TRUE)]

sigma0_g1 = sigma0_g[,grep("sigma0_g[",colnames(sigma0_g),fixed=TRUE)]
sigma0_g2 = sigma0_g[,grep("sigma0_g[2",colnames(sigma0_g),fixed=TRUE)]
sigma0_g3 = sigma0_g[,grep("sigma0_g[3",colnames(sigma0_g),fixed=TRUE)]

p_g1 = prob.inverse_g1 =matrix(NA, ncol=20,nrow=1000)

for(i in 1:2){
  
  # age 1 viabilities
  prob.inverse_g<- rnorm(length(mu0_g[,2]),mu0_g[,2],sigma0_g[,2])
  p_g <- boot::inv.logit(prob.inverse_g)
  hist(p_g,breaks=50,add=TRUE,col='red')
  
  prob.inverse_v1[i] ~ dnorm(mu0_v[1,i],tau0_v[1,i])
  logit(p_v1[i]) <- prob.inverse_v1[i]
  
  nu_1[i] = p_g1[i] + p_v1[i]*(1-p_g1[i])
  
  # age 2 viabilities
  prob.inverse_g2[i] ~ dnorm(mu0_g[2,i],tau0_g[2,i])
  logit(p_g2[i]) <- prob.inverse_g2[i]
  
  prob.inverse_v2[i] ~ dnorm(mu0_v[2,i],tau0_v[2,i])
  logit(p_v2[i]) <- prob.inverse_v2[i]
  
  nu_2[i] = p_g2[i] + p_v2[i]*(1-p_g2[i])
  
  # age 3 viabilities
  prob.inverse_g3[i] ~ dnorm(mu0_g[3,i],tau0_g[3,i])
  logit(p_g3[i]) <- prob.inverse_g3[i]
  
  prob.inverse_v3[i] ~ dnorm(mu0_v[3,i],tau0_v[3,i])
  logit(p_v3[i]) <- prob.inverse_v3[i]
  
  nu_3[i] = p_g3[i] + p_v3[i]*(1-p_g3[i])
  
  ## age 1 seed bag probabilities 
  # probability intact
  prob.inverse_1[i] ~ dnorm(mu0_1[1,i],tau0_1[1,i])
  logit(p_1[i]) <- prob.inverse_1[i]
  
  # probability emergence
  prob.inverse_2[i] ~ dnorm(mu0_2[1,i],tau0_2[1,i])
  logit(p_2[i]) <- prob.inverse_2[i]
  
  # probability conditional intact
  prob.inverse_3[i] ~ dnorm(mu0_3[1,i],tau0_3[1,i])
  logit(p_3[i]) <- prob.inverse_3[i]
  
  ## age 2 seed bag probabilities 
  # probability intact
  prob.inverse_4[i] ~ dnorm(mu0_1[2,i],tau0_1[2,i])
  logit(p_4[i]) <- prob.inverse_4[i]
  
  # probability emergence
  prob.inverse_5[i] ~ dnorm(mu0_2[2,i],tau0_2[2,i])
  logit(p_5[i]) <- prob.inverse_5[i]
  
  # probability conditional intact
  prob.inverse_6[i] ~ dnorm(mu0_3[2,i],tau0_3[2,i])
  logit(p_6[i]) <- prob.inverse_6[i]
  
  ## age 3 seed bag probabilities 
  # probability intact
  prob.inverse_7[i] ~ dnorm(mu0_1[3,i],tau0_1[3,i])
  logit(p_7[i]) <- prob.inverse_7[i]
  
  # probability emergence
  prob.inverse_8[i] ~ dnorm(mu0_2[3,i],tau0_2[3,i])
  logit(p_8[i]) <- prob.inverse_8[i]
  
  # probability conditional intact
  prob.inverse_9[i] ~ dnorm(mu0_3[3,i],tau0_3[3,i])
  logit(p_9[i]) <- prob.inverse_9[i]    
  
  ## obtain population-level estimates 
  s1[i] = p_1[i]*(p_2[i] + (1-p_2[i])*(nu_1[i])^(1/3))
  g1[i] = p_2[i]/(1-(1-(nu_1[i]^(1/3)))*(1-p_2[i]))
  s2[i] = p_3[i]*(nu_1[i]^(2/3))
  
  # conditional interpolation of viability
  nu_2c[i] = ifelse(nu_2[i] < nu_1[i], nu_1[i]^(2/3)*nu_2[i]^(1/3), nu_2[i]^(1/3)) 
  
  s3[i] = (p_4[i]*(p_5[i] + (1-p_5[i])*(nu_2c[i])^(1/3)))/(s1[i]*(1-g1[i]*s2[i]))    
  
  for(k in 1:n_yearBags){
    
    p.i0_g[i,k] ~ dnorm(mu_g[i,k],tau_g[i,k])
    logit(p0_g[i,k]) <- p.i0_g[i,k]
    
    p.i0_v[i,k] ~ dnorm(mu_v[i,k],tau_v[i,k])
    logit(p0_v[i,k]) <- p.i0_v[i,k]
    
    nu0_1[i,k] = p0_g[i,k] + p0_v[i,k]*(1-p0_g[i,k])
    
    p.i0_1[i,k] ~ dnorm(mu_1[i,k],tau_1[i,k])
    logit(p0_1[i,k]) <- p.i0_1[i,k]
    
    p.i0_2[i,k] ~ dnorm(mu_2[i,k],tau_2[i,k])
    logit(p0_2[i,k]) <- p.i0_2[i,k]
    
    p.i0_3[i,k] ~ dnorm(mu_3[i,k],tau_3[i,k])
    logit(p0_3[i,k]) <- p.i0_3[i,k]
    
    s1.0[i,k] = p0_1[i,k]*(p0_2[i,k] + (1-p0_2[i,k])*(nu0_1[i,k])^(1/3))
    g1.0[i,k] = p0_2[i,k]/(1-(1-(nu0_1[i,k]^(1/3)))*(1-p0_2[i,k]))
    s2.0[i,k] = p0_3[i,k]*(nu0_1[i,k]^(2/3))
    
  }
}

# interpolation
nu0_2 = ifelse(nu0_2 < nu0_1, nu0_1^(2/3)*nu0_2^(1/3), nu0_2^(1/3))

numerator = p0_4[i,k]*(p0_5[i,k] + (1-p0_5[i,k])*(nu0_2[i,k])^(1/3))
denominator = s1.0*(1-g1.0)*s2.0
s3.0[i,k] = numerator/denominator