##############################################################
#
# Rat tumor problem
#
#
rm(list=ls(all=TRUE)) # clear R environment

fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/singleSite/")
dir.create(file.path(fileDirectory), showWarnings = FALSE)

# seed data

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/"
simFiles <- paste0(directory,list.files(directory))

simData <- readRDS(simFiles[[17]])

selectedSite<-sample(unique(simData$site),1)
simData=simData[simData$site==selectedSite,]
selectedYear<-sample(unique(simData$year),1)
simData=simData[simData$year==selectedYear,]

y = simData$fruitingPlantNumber
n = simData$seedlingNumber
testDataset=data.frame(cbind(y,n))
saveRDS(testDataset, file=paste0(fileDirectory,"testDataset.rds"))


library(rjags)
library(MCMCvis)
# pass data to list for JAGS
data = list(
  y = y,
  n = n,
  N = length(y)
)

# -------------------------------------------------------------------
# Code for JAGS model
# -------------------------------------------------------------------

modelString1 <- "model { 

  # priors

   vfreq ~ dbeta(1,1)
   vcount ~ dgamma(1,1/20)
   alpha <- vfreq*vcount
   beta <- (1-vfreq)*vcount
   
  # Likelihood
  for(i in 1:N){
    theta[i] ~ dbeta(alpha,beta)
    y[i] ~ dbinom(theta[i], n[i])
  }
    
    p ~ dbeta(alpha,beta)
  
}"

modelString2 <- "model { 
  
   mu ~ dnorm(0,.001)
   sigma ~ dt(0, .04, 1)I(0,)
   alpha <- (mu^2-mu^3-mu*sigma^2)/(sigma^2)
   beta <- (mu-2*mu^2+mu^3-sigma^2+mu*sigma^2)/(sigma^2)
  
  # Likelihood
  for(i in 1:N){
    theta[i] ~ dbeta(alpha+0.01,beta+0.01)
    y[i] ~ dbinom(theta[i], n[i])
  }
    
    p ~ dbeta(alpha,beta)

  
}"


modelString3 <- "model { 

  # priors
  
   mu ~ dnorm(0,.001)
   sigma ~ dunif(0,20)
   tau <- 1/(sigma*sigma)
   
  
  # Likelihood
  for(i in 1:N){
    alpha[i] ~ dnorm(mu,tau)
    logit(theta[i]) <- alpha[i]
    y[i] ~ dbinom(theta[i], n[i])
    }
  
    p.i ~ dnorm(mu,tau)
    logit(p) <- p.i

}"


modelString4 <- "model { 

  # priors
  
   mu ~ dnorm(0,.001)
   sigma ~ dunif(0,20)
  
  # Likelihood
  for(i in 1:N){
    alpha.std[i] ~ dnorm(0,1)
    alpha[i]<- mu + sigma*alpha.std[i]
    logit(theta[i]) <- alpha[i]
    y[i] ~ dbinom(theta[i], n[i])
  }
    
    tau = 1/(sigma*sigma)
    p.i ~ dnorm(mu,tau)
    logit(p) <- p.i
  
}"

# -------------------------------------------------------------------
# Initial values
# -------------------------------------------------------------------
# set inits for rjags

mu.prior<-boot::logit(sum(y)/sum(n))
inits = list(list(vfreq = .5,vcount = 20),list(vfreq = .5,vcount = 20),list(vfreq = .5,vcount = 20))
inits2 = list(list(mu = .05,sigma = .05),list(mu = .05,sigma = .05),list(mu = .05,sigma = .05))
inits3 = list(list(mu = 0, sigma = 5),list(mu = 0, sigma = 5),list(mu = 0, sigma = 5))
inits4 = list(list(mu = mu.prior, sigma = 5),list(mu = mu.prior, sigma = 5),list(mu = mu.prior, sigma = 5))
# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# -------------------------------------------------------------------
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 1000
n.update = 5000
n.iterations = 9000
n.thin = 1
parsToMonitor = c("theta","vfreq","vcount","p")
parsToMonitor2 = c("theta","mu","sigma","p")
parsToMonitor3 = c("theta","alpha","mu","sigma","p")
parsToMonitor4 = c("theta","alpha","alpha.std","mu","sigma","p")

# -------------------------------------------------------------------
# Call to JAGS via rjags (model 1)
# -------------------------------------------------------------------
#set.seed(2)
# tuning (n.adapt)
jm1 = jags.model(textConnection(modelString1), inits=inits,
                 n.chains = length(inits), data = data, n.adapt = n.adapt)
jm2 = jags.model(textConnection(modelString2), inits=inits2,
                 n.chains = length(inits2), data = data, n.adapt = n.adapt)
jm3 = jags.model(textConnection(modelString3), inits=inits3,
                 n.chains = length(inits3), data = data, n.adapt = n.adapt)
jm4 = jags.model(textConnection(modelString4), inits=inits4,
                 n.chains = length(inits4), data = data, n.adapt = n.adapt)
# burn-in (n.update)
update(jm1, n.iterations = n.update)
update(jm2, n.iterations = n.update)
update(jm3, n.iterations = n.update)
update(jm4, n.iterations = n.update)

# chain (n.iter)
samples.rjags1 = coda.samples(jm1, variable.names = c(parsToMonitor), n.iter = n.iterations, thin = n.thin)
samples.rjags2 = coda.samples(jm2, variable.names = c(parsToMonitor2), n.iter = n.iterations, thin = n.thin)
samples.rjags3 = coda.samples(jm3, variable.names = c(parsToMonitor3), n.iter = n.iterations, thin = n.thin)
samples.rjags4 = coda.samples(jm4, variable.names = c(parsToMonitor4), n.iter = n.iterations, thin = n.thin)

saveRDS(samples.rjags1,file=paste0(fileDirectory,"samples.rjags1.rds"))
saveRDS(samples.rjags2,file=paste0(fileDirectory,"samples.rjags2.rds"))
saveRDS(samples.rjags3,file=paste0(fileDirectory,"samples.rjags3.rds"))
saveRDS(samples.rjags4,file=paste0(fileDirectory,"samples.rjags4.rds"))

