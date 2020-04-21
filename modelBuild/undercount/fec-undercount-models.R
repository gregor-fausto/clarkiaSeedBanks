# -------------------------------------------------------------------
# Fruits per plant (fec)
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)
# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(rjags) # jags interface
library(MCMCvis)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(HDInterval)

set.seed(10)
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
setwd("~/Dropbox/projects/clarkiaScripts/data/cleanData")
load(file="fecDF.RData")

### UNDERCOUNT SEEDLINGS MODEL

# use ALL observations if modeling seedling counts as a latent state
dat <- fecDF %>% 
  dplyr::filter(year == 2009) %>%
  dplyr::select(c(year,site,transect,position,noSeedlings,noFruitingPlants)) %>%
  unique

dat$site.index <- as.integer(dat$site)

data = list(
  n = as.double(dat$noSeedlings),
  #y = as.double(dat$noFruitingPlants),
  N = nrow(dat)
)

inits = list(
  list( pi = runif(1, 0, 1) , z = rep( 600, nrow(dat)), lambda = rep( 1, nrow(dat)), alpha = .1, upsilon = 5),
  list( pi = runif(1, 0, 1) , z = rep( 600, nrow(dat)), lambda = rep( 5, nrow(dat)), alpha = 2, upsilon = .5),
  list( pi = runif(1, 0, 1) ,z = rep( 600, nrow(dat)), lambda = rep( 3, nrow(dat)), alpha = 10, upsilon = 5) 
)


setwd("~/Dropbox/projects/clarkiaScripts/bayesian-attempts")
sink("fecUndercountSeedlingsJAGS.R")
cat("
    model { 
    
    # hyperpriors
    alpha ~ dgamma(0.001, 0.001)
    upsilon ~ dgamma(0.001, 0.001)
    
    # priors
    pi ~ dbeta(1, 1) 
    for(i in 1:N){ lambda[i] ~ dgamma(alpha, upsilon) } 
    
    # likelihood
    for (i in 1:N){ 
    z[i] ~ dpois(lambda[i])
    n[i] ~ dbinom(1-pi, z[i])
    } 
    
    }
    ", fill = TRUE)
sink()


# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 1000
n.update = 10000
n.iter = 100000

# Call to JAGS

# tuning (n.adapt)
jm.latent = jags.model("fecUndercountSeedlingsJAGS.R", data = data, inits = inits,
                       n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm.latent, n.iter = n.update)

# chain (n.iter)
zc.latent = coda.samples(jm.latent, variable.names = c("pi", "alpha", "upsilon", "lambda","z"), n.iter = n.iter)

# summary table
MCMCsummary(zc.latent, n.eff = TRUE, params = c("pi", "alpha", "upsilon"))

# trace plot
MCMCtrace(zc.latent, n.eff = TRUE, params = c("pi", "alpha", "upsilon"), pdf = FALSE)

HPDI <- MCMCpstr(zc.latent, params = c("z"), func = function(x)
  hdi(x, .95))
med <- MCMCpstr(zc.latent, params = c("z"), func = median)
lambda <- MCMCpstr(zc.latent, params = c("lambda"), func = median)
lambda.hpdi <- MCMCpstr(zc.latent, params = c("lambda"), func = function(x)
  hdi(x, .95))

df<-data.frame(index = seq(1,length(dat$noSeedlings),1),
               n = dat$noSeedlings, 
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
