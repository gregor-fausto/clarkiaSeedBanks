# with fewer than n groups, mu.alpha becomes poorly constrained regardless
# of whether it's in the model mean or not

# model 2 > model 1 for few groups

# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(rjags) 
library(R2jags)
library(MCMCvis)

# -------------------------------------------------------------------
# Simulate data
# -------------------------------------------------------------------
set.seed(10)

N.plots = 4
trials=400

n = rep(100,trials)
N=length(n)
plotReps=N/N.plots

# Block 1
p<-c(.2,.4,.6,.8)
b0<-rep(c(.2,.4,.6,.8),each=plotReps)

y = rbinom(trials, 100, p = c(b0))

# vectors indexing plots and blocks
plot = rep(1:N.plots,each=plotReps)

# pass data to list for JAGS
data = list(
  y = y,
  n = n,
  N = length(n),
  plot = plot,
  N.plots = N.plots
)

plot<-as.factor(plot)

# create your model matrix
my_matrix <- model.matrix(y ~ 0 + plot )
my_matrix %>% View

dat <- data.frame(y,n,plot=as.factor(plot))

# -------------------------------------------------------------------
# Code for JAGS model
# -------------------------------------------------------------------

modelString1 <- "model { 
  ## Priors
  
  # hyperpriors
  mu.alpha ~ dnorm(0, 0.0001)
  
  sigma.plot ~ dunif(0,100) 
  tau.plot <- 1 / sigma.plot^2
  
  # priors 
  for(i in 1:N.plots){     
    eps.plot[i]~dnorm(0,tau.plot)
  }
  
  # Likelihood
  for(i in 1:N){
    logit(p[i]) <- mu.alpha + eps.plot[plot[i]]
    y[i] ~ dbin(p[i], n[i])
    
  }
}"

modelString2 <- "model { 
  ## Priors
  
  # hyperpriors
  mu.alpha ~ dnorm(0, 0.0001)
  
  sigma.plot ~ dunif(0,100) 
  tau.plot <- 1 / sigma.plot^2
  
  # priors 
  for(i in 1:N.plots){     
    eps.plot[i]~dnorm(mu.alpha,tau.plot)
  }
  
  # Likelihood
  for(i in 1:N){
    logit(p[i]) <- eps.plot[plot[i]]
    y[i] ~ dbin(p[i], n[i])
    
  }
}"

# https://github.com/stan-dev/example-models/blob/master/ARM/Ch.12/radon_intercept.stan
modelString3 <- "model { 
  ## Priors
  
  # hyperpriors
  mu.alpha ~ dnorm(0, 0.0001)
  #sigma.alpha ~ dunif(0,100)
  #tau.alpha <- 1/sigma.alpha^2
  
  # priors 
  for(i in 2:N.plots){     
    b0[i]~dnorm(0,0.0001)
  }
  
  # sum to zero constraint
  b0[1] <- -sum( b0[2:N.plots])  
  
  # Likelihood
  for(i in 1:N){
    logit(p[i]) <- mu.alpha + b0[plot[i]]
    y[i] ~ dbin(p[i], n[i])
    
  }
}"



# -------------------------------------------------------------------
# Initial values
# -------------------------------------------------------------------
# set inits for rjags
inits = list(list(mu.alpha = 0,sigma.plot=2),
             list(mu.alpha = 0,sigma.plot=2),
             list(mu.alpha = 0,sigma.plot=2)) 

# set inits function for R2jags
initsFun<-function(){list(
  mu.alpha=0,
  sigma.plot=2
)}

# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# -------------------------------------------------------------------
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 500
n.update = 5000
n.iterations = 10000
n.thin = 1
parsToMonitor = c("mu.alpha","sigma.plot","eps.plot")

# -------------------------------------------------------------------
# Call to JAGS via rjags (model 1)
# -------------------------------------------------------------------
set.seed(2)
# tuning (n.adapt)
jm1 = jags.model(textConnection(modelString1), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm1, n.iterations = n.update)

# chain (n.iter)
samples.rjags1 = coda.samples(jm1, variable.names = c(parsToMonitor), n.iter = n.iterations, thin = n.thin)

# -------------------------------------------------------------------
# Call to JAGS via rjags (model 2)
# -------------------------------------------------------------------
set.seed(2)
# tuning (n.adapt)
jm2 = jags.model(textConnection(modelString2), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm2, n.iterations = n.update)

# chain (n.iter)
samples.rjags2 = coda.samples(jm2, variable.names = c(parsToMonitor), n.iter = n.iterations, thin = n.thin)

# -------------------------------------------------------------------
# Call to JAGS via rjags (model 3)
# -------------------------------------------------------------------
# set inits for rjags
inits = list(list(mu.alpha = 0,b0 = rep(0,data$N.plots)),
             list(mu.alpha = 0,b0 = rep(0,data$N.plots)),
             list(mu.alpha = 0,b0 = rep(0,data$N.plots)) )

parsToMonitor = c("mu.alpha","sigma.alpha","b0")

set.seed(2)
# tuning (n.adapt)
jm3 = jags.model(textConnection(modelString3), data = data, inits = inits,
                 n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm3, n.iterations = n.update)

# chain (n.iter)
samples.rjags3 = coda.samples(jm3, variable.names = c(parsToMonitor), n.iter = n.iterations, thin = n.thin)

# -------------------------------------------------------------------
# Create data frame of summaries and rhat values
# -------------------------------------------------------------------
sum.rjags1 <- MCMCvis::MCMCsummary(samples.rjags1,params=c("mu.alpha","eps.plot","sigma.plot"))
sum.rjags1

sum.rjags2 <- MCMCvis::MCMCsummary(samples.rjags2,params=c("mu.alpha","eps.plot","sigma.plot"))
sum.rjags2

sum.rjags3 <- MCMCvis::MCMCsummary(samples.rjags3,params=c("mu.alpha","sigma.alpha","b0"))
sum.rjags3

# -------------------------------------------------------------------
# Check convergence
# -------------------------------------------------------------------
MCMCvis::MCMCtrace(samples.rjags1,params=c("mu.alpha"),pdf=FALSE)
MCMCvis::MCMCtrace(samples.rjags2,params=c("mu.alpha"),pdf=FALSE)
MCMCvis::MCMCtrace(samples.rjags3,params=c("mu.alpha"),pdf=FALSE)

par(mfrow=c(1,1))
# model 1
mat <- matrix(NA,nrow=dim(MCMCvis::MCMCchains(samples.rjags1,params=c("mu.alpha")))[1],ncol=N.plots,)
for(i in 1:N.plots){
  mat[,i]<-MCMCvis::MCMCchains(samples.rjags1,params=c("mu.alpha"))+MCMCvis::MCMCchains(samples.rjags1,params=c("eps.plot"))[,i]
}

plot(apply(apply(mat,2,boot::inv.logit),2,mean),rep(c(.2,.4,.6,.8),N.plots/length(p)))
abline(a=0,b=1)
# model 2
mat <- matrix(NA,nrow=dim(MCMCvis::MCMCchains(samples.rjags2,params=c("mu.alpha")))[1],ncol=N.plots,)
for(i in 1:N.plots){
  mat[,i]<-MCMCvis::MCMCchains(samples.rjags2,params=c("mu.alpha"))+MCMCvis::MCMCchains(samples.rjags2,params=c("eps.plot"))[,i]
}

plot(apply(apply(mat,2,boot::inv.logit),2,mean),rep(c(.2,.4,.6,.8),N.plots/length(p)))
abline(a=0,b=1)

# model 3
mat <- matrix(NA,nrow=dim(MCMCvis::MCMCchains(samples.rjags3,params=c("b0")))[1],ncol=N.plots,)
for(i in 1:N.plots){
  mat[,i]<-MCMCvis::MCMCchains(samples.rjags3,params=c("mu.alpha"))+MCMCvis::MCMCchains(samples.rjags3,params=c("b0"))[,i]
}

plot(apply(apply(mat,2,boot::inv.logit),2,mean),rep(c(.2,.4,.6,.8),N.plots/length(p)))
abline(a=0,b=1)

