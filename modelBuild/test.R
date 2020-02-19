# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(rjags) # jags interface
library(MCMCvis)

# data for testing

N.plots = 8

trials=400

set.seed(10)
n = rep(100,trials)
N=length(n)
reps=N/N.plots

y = rbinom(trials, 100, p = c(c(rep(.25,reps),rep(.75,reps),rep(.9,reps),rep(.1,reps))+.05,c(rep(.25,reps),rep(.75,reps),rep(.9,reps),rep(.1,reps))-.05))
plot = c(rep(1,reps),rep(2,reps),rep(3,reps),rep(4,reps),rep(5,reps),rep(6,reps),rep(7,reps),rep(8,reps))
block = c(rep(1,trials/2),rep(2,trials/2))
# pass data to list for JAGS
data = list(
  y = y,
  n = n,
  N = length(n),
  plot = plot,
  block= block,
  N.plots = N.plots,
  N.blocks = 2
)


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 3000
n.update = 5000
n.iter = 1000
n.thin = 1

set.seed(10)
dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/jagsScripts/")


# set inits for JAGS
inits = list(list(mu.alpha = rnorm(1),sigma.plot=rlnorm(1)),
             list(mu.alpha = rnorm(1),sigma.plot=rlnorm(1)),
             list(mu.alpha = rnorm(1),sigma.plot=rlnorm(1))
)

# Call to JAGS

# tuning (n.adapt)
jm = jags.model(paste0(dir,"test2JAGS.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

parsToMonitor = c("mu.alpha","eps.plot","sigma.plot",)
# chain (n.iter)
zc = coda.samples(jm, variable.names = c(parsToMonitor), n.iter = n.iter, thin = n.thin)

#save(zc,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBagsEvansModelFit.rds")

MCMCsummary(zc,params=c("mu.alpha","eps.plot","sigma.plot","sigma.block","eps.block"))
MCMCtrace(zc,params=c("mu.alpha"))

hist(boot::inv.logit(MCMCchains(zc,params=c("mu.alpha"))+MCMCchains(zc,params=c("eps.plot"))[,2]) ,breaks=100); abline(v=.75,col='red')

apply(boot::inv.logit(MCMCchains(zc,params=c("eps.plot"))),2,mean)
apply(boot::inv.logit(MCMCchains(zc,params=c("mu.alpha"))),2,mean)
      
apply(boot::inv.logit(MCMCchains(zc,params=c("eps.plot"))),2,mean)

