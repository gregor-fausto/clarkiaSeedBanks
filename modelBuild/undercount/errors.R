library(tidybayes)
library(rjags)
library(magrittr) # for %<>% with recover_types
library(MCMCvis)
library(HDInterval)
library(ggplot2)

censusSeedlingsFruitingPlants <- readRDS("~/Dropbox/dataLibrary/postProcessingData/censusSeedlingsFruitingPlants.RDS")

data <- censusSeedlingsFruitingPlants %>% 
  dplyr::filter(site=="BG"&year==2007) 

data %>% dplyr::group_by(transect) %>%
  dplyr::summarise(sum(fruitplNumber)/sum(seedlingNumber))

# seedlingNumber = rpois(1000,20)
# y = rbinom(n,size=n,prob=.25)
# hist(y/n,breaks=20);abline(v=0.25,col='red')

dataJags <- tidybayes::compose_data(data)

inits = list(
  list( pi = rep(.25, 4) , z = rep( 100, 30), alpha = .1, upsilon = .1),
  list( pi = rep(.5, 4) , z = rep( 100, 30), alpha = 2, upsilon = .5),
  list( pi = rep(.75, 4) ,z = rep( 100, 30), alpha = 10, upsilon = .9) 
)

setwd("~/Dropbox/clarkiaSeedBanks/modelBuild/undercount")
sink("model-survival.R")
  cat("

model{
  
    # hyperpriors
    alpha ~ dgamma(0.001, 0.001)
    upsilon ~ dgamma(0.001, 0.001)
    
    # priors
    # place a (0,.5) prior on pi
    # assuming that > half of seedlings are counted
    # as parameter decreases, true number of seedlings = observed number
    # as parameter increases, true number of seedlings > observed number (max p = 0)
    # as parameter = .5 half of seeldings are observed
    for(i in 1:n_transect){ pi[i] ~ dbeta(1,1) }
    phi ~ dbeta(1, 1) 
    for(i in 1:n_transect){ lambda[i] ~ dgamma(alpha, upsilon) } 
    
    # likelihood
    for (i in 1:n){ 
    # process model for true number of seedlings drawn from lambda
    z[i] ~ dpois(lambda[transect[i]])
    
    # data model
    # counts of seedling number related to true number of seedlings (z[i])
    # via detection probability (1-pi)
    seedlingNumber[i] ~ dbinom(1-pi[transect[i]], z[i])
    
    # data model 
    # counts of fruiting plants relating true number of seedlings (z[i])
    # via survival probability (phi)
    fruitplNumber[i] ~ dbinom(phi, z[i])
    } 
    
  # end of data model



} #end of model

",fill=TRUE)
  sink()

  # scalars that specify the 
  # number of iterations in the chain for adaptation
  # number of iterations for burn-in
  # number of samples in the final chain
  n.adapt = 1000
  n.update = 10000
  n.iter = 10000
  
  # Call to JAGS
  

  
  # tuning (n.adapt)
  jm = jags.model("~/Dropbox/clarkiaSeedBanks/modelBuild/undercount/model-survival.R", 
                  data = dataJags,  n.adapt = n.adapt, inits = inits,
                  n.chains = length(inits))
  
  # burn-in (n.update)
  update(jm, n.iter = n.update)
  
  # chain (n.iter)
  zc = coda.samples(jm, variable.names = c("z","lambda","alpha","pi","phi","upsilon"), n.iter = n.iter)
  
  # summary table
  MCMCsummary(zc, n.eff = TRUE, params = c("z","lambda","alpha","pi","phi","upsilon"))
  
  cbind(data,p=(apply(MCMCchains(zc,params='pi'),2,median)))
  
  z.sim<-MCMCchains(zc,params="z.sim")
  hist(z.sim[,1],breaks=100);abline(v=dataJags$seedlingNumber[1],col='red')
  
  z.50<-apply(MCMCchains(zc, params="z"),2,median)
  lambda.50<-apply(MCMCchains(zc, params="lambda"),2,median)
  
plot(z.50,dataJags$seedlingNumber,ylim=c(0,100),xlim=c(0,100)) 
points(lambda.50,dataJags$seedlingNumber,pch=16) 
abline(a=0,b=1)

sum<-data %>%
  dplyr::mutate(surv = fruitplNumber/seedlingNumber)


est <- sum %>%
  dplyr::summarise(n = sum(seedlingNumber,na.rm=TRUE),
                   y = sum(fruitplNumber,na.rm=TRUE)) %>%
  dplyr::mutate(y/n)

hist(sum$surv,breaks=30);abline(v=est,col='red',lwd=2)


plot(sum$seedlingNumber,sum$surv)

cbind(lambda.50,dataJags$seedlingNumber)


HPDI <- MCMCpstr(zc, params = c("z"), func = function(x) hdi(x, .95))
med <- MCMCpstr(zc, params = c("z"), func = median)
lambda <- MCMCpstr(zc, params = c("lambda"), func = median)
lambda.hpdi <- MCMCpstr(zc, params = c("lambda"), func = function(x)
  hdi(x, .95))

df<-data.frame(index = seq(1,length(dataJags$seedlingNumber),1),
               n = dataJags$seedlingNumber, 
               z = med$z,
               lo = HPDI$z[,1], hi = HPDI$z[,2])
               # lambda = lambda$lambda,
               # lambda.lo = lambda.hpdi$lambda[,1],
               # lambda.hi = lambda.hpdi$lambda[,2]) 

ggplot(df) +
 # geom_pointrange(aes(x = index, y = lambda, ymin = lambda.lo, ymax = lambda.hi), col = "orange",alpha = 0.3) +
  geom_pointrange(aes(x = index, y = z, ymin = lo, ymax = hi), col = "blue",alpha = 0.3) +
  geom_point(aes(x = index, y = n)) + 
  theme_minimal() 

ggplot(df) +
  geom_abline(intercept = 0, slope =1 ) +
  geom_pointrange(aes(x = n, y = lambda, ymin = lambda.lo, ymax = lambda.hi), col = "orange",alpha = 0.3) +
  geom_pointrange(aes(x = n, y = z, ymin = lo, ymax = hi), col = "blue",alpha = 0.3) +
  theme_minimal() 


### UNDERCOUNT SEEDLINGS MODEL 2

# use ALL observations if modeling seedling counts as a latent state
dat <- fecDF %>% 
  dplyr::filter(year == 2009) %>%
  dplyr::select(c(year,site,transect,position,noSeedlings,noFruitingPlants)) %>%
  unique %>% 
  as.data.frame

dat$site.index <- as.integer(dat$site)

data = list(
  n = as.double(dat$noSeedlings),
  y = as.double(dat$noFruitingPlants),
  site = as.double(dat$site.index),
  N = nrow(dat),
  n.sites = nlevels(dat$site)
)

inits = list(
  list( pi = rep(.1, 20) , z = rep( 528, nrow(dat)), alpha = .1, upsilon = .1),
  list( pi = rep(.15, 20) , z = rep( 528, nrow(dat)), alpha = 2, upsilon = .5),
  list( pi = rep(.25, 20) ,z = rep( 528, nrow(dat)), alpha = 10, upsilon = .9) 
)


setwd("~/Dropbox/projects/clarkiaScripts/bayesian-attempts")
sink("fecUndercountSeedlings2JAGS.R")
cat("
    model { 
    
    # hyperpriors
    alpha ~ dgamma(0.001, 0.001)
    upsilon ~ dgamma(0.001, 0.001)
    
    # priors
    for(j in 1:n.sites){ pi[j] ~ dunif(0,.5) }
    for(j in 1:n.sites){ phi[j] ~ dbeta(1, 1) }
    for(i in 1:N){ lambda[i] ~ dgamma(alpha, upsilon) } 
    
    # likelihood
    for (i in 1:N){ 
    z[i] ~ dpois(lambda[i])
    n[i] ~ dbinom(1-pi[site[i]], z[i])
    
    # data nidek
    y[i] ~ dbinom(phi[site[i]], z[i])
    } 
    
    }
    ", fill = TRUE)
sink()


# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 1000
n.update = 1000
n.iter = 10000

# Call to JAGS

# tuning (n.adapt)
jm.latent = jags.model("fecUndercountSeedlings2JAGS.R", data = data, inits = inits,
                       n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm.latent, n.iter = n.update)

# chain (n.iter)
zc.latent = coda.samples(jm.latent, variable.names = c("pi", "alpha", "upsilon", "phi", "lambda","z"), n.iter = n.iter)

# summary table
MCMCsummary(zc.latent, n.eff = TRUE, params = c("pi", "phi", "alpha", "upsilon"))

# trace plot
MCMCtrace(zc.latent, n.eff = TRUE, params = c("alpha", "upsilon"), pdf = FALSE)

HPDI <- MCMCpstr(zc.latent, params = c("z"), func = function(x)
  hdi(x, .95))
med <- MCMCpstr(zc.latent, params = c("z"), func = median)
lambda <- MCMCpstr(zc.latent, params = c("lambda"), func = median)
lambda.hpdi <- MCMCpstr(zc.latent, params = c("lambda"), func = function(x)
  hdi(x, .95))

df<-data.frame(index = seq(1,length(dat$noSeedlings),1),
               n = dat$noSeedlings, 
               y = dat$noFruitingPlants,
               z = med$z,
               lo = HPDI$z[,1], hi = HPDI$z[,2], 
               lambda = lambda$lambda,
               lambda.lo = lambda.hpdi$lambda[,1],
               lambda.hi = lambda.hpdi$lambda[,2],
               site = dat$site)

ggplot(df) +
  geom_pointrange(aes(x = index, y = lambda, ymin = lambda.lo, ymax = lambda.hi), col = "red",alpha = 0.3) +
  geom_pointrange(aes(x = index, y = z, ymin = lo, ymax = hi), col = "blue",alpha = 0.3) +
  geom_point(aes(x = index, y = n)) + 
  theme_minimal() +
  facet_wrap(~site, scales = "free")

ggplot(df) +
  geom_pointrange(aes(x = n, y = lambda, ymin = lambda.lo, ymax = lambda.hi), col = "red",alpha = 0.3) +
  geom_pointrange(aes(x = n, y = z, ymin = lo, ymax = hi), col = "blue",alpha = 0.3) +
  geom_abline(intercept = 0, slope =1 ) +
  theme_minimal() +
  facet_wrap(~site, scales = "free")
