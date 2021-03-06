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

N.plots = 8
N.blocks = 2
trials=400

n = rep(100,trials)
N=length(n)
plotReps=N/N.plots
blockReps=N/N.blocks

# Block 1
b1<-rep(c(.2,.4,.6,.8),each=plotReps)
# Block 2
b2<-rep(c(.2,.4,.6,.8),each=plotReps)

y = rbinom(trials, 100, p = c(b1,b2))

# vectors indexing plots and blocks
plot = rep(1:8,each=plotReps)
#block = rep(1:2,each=blockReps)

# pass data to list for JAGS
data = list(
  y = y,
  n = n,
  N = length(n),
  plot = plot,
  block= block,
  N.plots = N.plots,
  N.blocks = N.blocks
)

# -------------------------------------------------------------------
# Code for JAGS model
# -------------------------------------------------------------------

modelString <- "model { 
  ## Priors
  
  # hyperpriors
  mu.alpha ~ dnorm(0, 0.0001)
  
  sigma.plot ~ dunif(0,100) 
  tau.plot <- 1 / sigma.plot^2
  
 # sigma.block ~ dunif(0,100) 
 # tau.block <- 1 / sigma.block^2
  
  # priors 
  for(i in 1:N.plots){     
    eps.plot[i]~dnorm(0,tau.plot)
  }
  
 # for(i in 1:N.blocks){
#    eps.block[i]~dnorm(0,tau.block)
 # }
  
  # Likelihood
  for(i in 1:N){
    logit(p[i]) <- mu.alpha + eps.plot[plot[i]]#+ eps.block[block[i]]
    y[i] ~ dbin(p[i], n[i])
    
  }
}"

# -------------------------------------------------------------------
# Initial values
# -------------------------------------------------------------------
# set inits for rjags
inits = list(list(mu.alpha = 0,sigma.plot=2,sigma.block=2),
             list(mu.alpha = 0,sigma.plot=2,sigma.block=2),
             list(mu.alpha = 0,sigma.plot=2,sigma.block=2)) 

# set inits function for R2jags
initsFun<-function(){list(
  mu.alpha=0,
  sigma.plot=2,
  sigma.block=2
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
parsToMonitor = c("mu.alpha","sigma.plot","sigma.block","eps.plot","eps.block")

# -------------------------------------------------------------------
# Call to JAGS via rjags
# -------------------------------------------------------------------
set.seed(2)
# tuning (n.adapt)
jm = jags.model(textConnection(modelString), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iterations = n.update)

# chain (n.iter)
samples.rjags = coda.samples(jm, variable.names = c(parsToMonitor), n.iter = n.iterations, thin = n.thin)

# -------------------------------------------------------------------
# Call to JAGS via R2jags
# -------------------------------------------------------------------
set.seed(2)
samples.R2jags <-jags(data=data,inits=initsFun,parameters.to.save=parsToMonitor,model.file=textConnection(modelString),
                   n.thin=n.thin,n.chains=length(inits),n.burnin=n.adapt,n.iter=n.iterations,DIC=T)

# -------------------------------------------------------------------
# Create data frame of summaries and rhat values
# -------------------------------------------------------------------
sum.rjags <- MCMCvis::MCMCsummary(samples.rjags,params=c("mu.alpha","eps.plot","sigma.plot","sigma.block","eps.block"))
sum.rjags

sum.R2jags2 <- MCMCvis::MCMCsummary(samples.R2jags,params=c("mu.alpha","eps.plot","sigma.plot","sigma.block","eps.block"))
sum.R2jags2


rhat.rjags <- data.frame(variable = rownames(sum.rjags),Rhat.rjags = sum.rjags$Rhat)

sum.R2jags <- data.frame(samples.R2jags$BUGSoutput$summary)[-1,]
rhat.R2jags <- data.frame(variable = rownames(sum.R2jags),Rhat.R2jags = sum.R2jags$Rhat)
sum.R2jags



library(dplyr)
compareDf<-rhat.rjags %>% dplyr::left_join(rhat.R2jags,by="variable")
compareDf

# -------------------------------------------------------------------
# Compare traceplots
# -------------------------------------------------------------------
library(MCMCvis)
par(mfrow=c(2,2))
MCMCvis::MCMCtrace(samples.rjags,params="mu.alpha",pdf=FALSE)
MCMCvis::MCMCtrace(samples.R2jags,params="mu.alpha",pdf=FALSE)

apply(apply(MCMCchains(samples.rjags,params="mu.alpha"),2,boot::inv.logit),2,mean)
apply(MCMCchains(samples.rjags,params="eps.plot"),2,mean)
apply(MCMCchains(samples.rjags,params="eps.block"),2,mean)

model<-lme4::glmer(cbind(y,n) ~ 1 + (1|plot) + (1|block),family="binomial")
lme4::ranef(model)

library(ggplot2)
dat<-data.frame(y,n,plot,block)
ggplot(dat,aes(x=as.factor(plot),y=y/n,col=as.factor(block))) + 
  geom_jitter() + geom_boxplot(alpha=0.2)
