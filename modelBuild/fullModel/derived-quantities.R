## write code to work with derived quantities for seed bags
library(MCMCvis)
samples = readRDS("/Users/Gregor/Dropbox/dataLibrary/workflow/samples/belowgroundSamples.RDS")

mu0_g <- MCMCchains(samples,params="mu0_g")
mu0_v <- MCMCchains(samples,params="mu0_v")
mu0_1 <- MCMCchains(samples,params="mu0_1")
mu0_2 <- MCMCchains(samples,params="mu0_2")
mu0_3 <- MCMCchains(samples,params="mu0_3")

sigma0_g <- MCMCchains(samples,params="sigma0_g")
sigma0_v <- MCMCchains(samples,params="sigma0_v")
sigma0_1 <- MCMCchains(samples,params="sigma0_1")
sigma0_2 <- MCMCchains(samples,params="sigma0_2")
sigma0_3 <- MCMCchains(samples,params="sigma0_3")

p_g1 = prob.inverse_g1 = matrix(NA, ncol=20,nrow=1000)
p_v1 = prob.inverse_v1 = matrix(NA, ncol=20,nrow=1000)
nu_1 = matrix(NA,ncol=20,nrow=1000)
p_1 = prob.inverse_1 = matrix(NA, ncol=20,nrow=1000)
p_2 = prob.inverse_2 = matrix(NA, ncol=20,nrow=1000)
p_3 = prob.inverse_3 = matrix(NA, ncol=20,nrow=1000)

s1 = g1 = s2 = matrix(NA, ncol=20,nrow=1000)


for(i in 1:20){
  
  # age 1 viabilities
  prob.inverse_g1[,i] <- rnorm(length(mu0_g[,i]),mu0_g[,i],sigma0_g[,i])
  p_g1[,i] <- boot::inv.logit(prob.inverse_g1[,i])

  prob.inverse_v1[,i] <- rnorm(length(mu0_v[,i]),mu0_v[,i],sigma0_v[,i])
  p_v1[,i] <- boot::inv.logit(prob.inverse_v1[,i])
  
  nu_1[,i] = p_g1[,i] + p_v1[,i]*(1-p_g1[,i])
  
  ## age 1 seed bag probabilities 
  # probability intact
  prob.inverse_1[,i] <- rnorm(length(mu0_1[,i]),mu0_1[,i],sigma0_1[,i])
  p_1[,i] <-  boot::inv.logit(prob.inverse_1[,i])
  
  # probability emergence
  prob.inverse_2[,i] <- rnorm(length(mu0_2[,i]),mu0_2[,i],sigma0_2[,i])
  p_2[,i] <-  boot::inv.logit(prob.inverse_2[,i])
  
  # probability conditional intact
  prob.inverse_3[,i] <- rnorm(length(mu0_3[,i]),mu0_3[,i],sigma0_3[,i])
  p_3[,i] <-  boot::inv.logit(prob.inverse_3[,i])

  ## obtain population-level estimates 
  s1[,i] = p_1[,i]*(p_2[,i] + (1-p_2[,i])*(nu_1[,i])^(1/3))
  g1[,i] = p_2[,i]/(1-(1-(nu_1[,i]^(1/3)))*(1-p_2[,i]))
  s2[,i] = p_3[,i]*(nu_1[,i]^(2/3))
  
  # # conditional interpolation of viability
  # nu_2c[i] = ifelse(nu_2[i] < nu_1[i], nu_1[i]^(2/3)*nu_2[i]^(1/3), nu_2[i]^(1/3)) 
  # 
  # s3[i] = (p_4[i]*(p_5[i] + (1-p_5[i])*(nu_2c[i])^(1/3)))/(s1[i]*(1-g1[i]*s2[i]))    
  # 
  # for(k in 1:n_yearBags){
  #   
  #   p.i0_g[i,k] ~ dnorm(mu_g[i,k],tau_g[i,k])
  #   logit(p0_g[i,k]) <- p.i0_g[i,k]
  #   
  #   p.i0_v[i,k] ~ dnorm(mu_v[i,k],tau_v[i,k])
  #   logit(p0_v[i,k]) <- p.i0_v[i,k]
  #   
  #   nu0_1[i,k] = p0_g[i,k] + p0_v[i,k]*(1-p0_g[i,k])
  #   
  #   p.i0_1[i,k] ~ dnorm(mu_1[i,k],tau_1[i,k])
  #   logit(p0_1[i,k]) <- p.i0_1[i,k]
  #   
  #   p.i0_2[i,k] ~ dnorm(mu_2[i,k],tau_2[i,k])
  #   logit(p0_2[i,k]) <- p.i0_2[i,k]
  #   
  #   p.i0_3[i,k] ~ dnorm(mu_3[i,k],tau_3[i,k])
  #   logit(p0_3[i,k]) <- p.i0_3[i,k]
  #   
  #   s1.0[i,k] = p0_1[i,k]*(p0_2[i,k] + (1-p0_2[i,k])*(nu0_1[i,k])^(1/3))
  #   g1.0[i,k] = p0_2[i,k]/(1-(1-(nu0_1[i,k]^(1/3)))*(1-p0_2[i,k]))
  #   s2.0[i,k] = p0_3[i,k]*(nu0_1[i,k]^(2/3))
  #   
  # }
}

hist(s1[,1],breaks=30)
hist(MCMCchains(samples,params='s1')[,1],breaks=30)

plot(density(s1[,1]),type='n',ylim=c(0,15))
for(i in 1:20){
  lines(density(nu_1[,i]),col='lightgray')
}

# interpolation
nu0_2 = ifelse(nu0_2 < nu0_1, nu0_1^(2/3)*nu0_2^(1/3), nu0_2^(1/3))

numerator = p0_4[i,k]*(p0_5[i,k] + (1-p0_5[i,k])*(nu0_2[i,k])^(1/3))
denominator = s1.0*(1-g1.0)*s2.0
s3.0[i,k] = numerator/denominator



mu_g <- MCMCchains(samples,params="mu_g")
mug <- list()
mug[[1]] = mu_g[,grep(",1]",colnames(mu_g),fixed=TRUE)]
mug[[2]] = mu_g[,grep(",2]",colnames(mu_g),fixed=TRUE)]
mug[[3]] = mu_g[,grep(",3]",colnames(mu_g),fixed=TRUE)]

sigma_g <- MCMCchains(samples,params="sigma_g")
sigmag <- list()
sigmag[[1]] = sigma_g[,grep(",1]",colnames(sigma_g),fixed=TRUE)]
sigmag[[2]] = sigma_g[,grep(",2]",colnames(sigma_g),fixed=TRUE)]
sigmag[[3]] = sigma_g[,grep(",3]",colnames(sigma_g),fixed=TRUE)]
