# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(rjags) # jags interface
library(MCMCvis)

# data for testing

set.seed(10)
n = rep(100,100)
reps=length(n)/2

y = rbinom(100, 100, p = c(rep(.25,50),rep(.75,50)))
N.plots = 2
plot = c(rep(1,reps),rep(2,reps))
# pass data to list for JAGS
data = list(
  y = y,
  n = n,
  N = length(n),
  plot = plot,
  N.plots = 2
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
inits = list(list(mu.alpha = 0,sigma.plot=10,eps.plot = rep(0,data$N.plots)),
             list(mu.alpha = 0,sigma.plot=20,eps.plot = rep(1,data$N.plots)),
             list(mu.alpha = 0,sigma.plot=12,eps.plot = rep(-1,data$N.plots))
)

# Call to JAGS

# tuning (n.adapt)
jm = jags.model(paste0(dir,"test2JAGS.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

parsToMonitor = c("mu.alpha","eps.plot","sigma.plot")
# chain (n.iter)
zc = coda.samples(jm, variable.names = c(parsToMonitor), n.iter = n.iter, thin = n.thin)

#save(zc,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBagsEvansModelFit.rds")

MCMCsummary(zc,params=c("mu.alpha","eps.plot","sigma.plot"))
MCMCtrace(zc,params=c("mu.alpha"))

hist(boot::inv.logit(MCMCchains(zc,params=c("mu.alpha"))+MCMCchains(zc,params=c("eps.plot"))[,2]) ,breaks=100); abline(v=.75,col='red')

apply(boot::inv.logit(MCMCchains(zc,params=c("eps.plot"))),2,mean)
apply(boot::inv.logit(MCMCchains(zc,params=c("mu.alpha"))),2,mean)
      
apply(boot::inv.logit(MCMCchains(zc,params=c("eps.plot"))),2,mean)

