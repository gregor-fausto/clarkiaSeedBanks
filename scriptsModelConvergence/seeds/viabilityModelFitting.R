# list of JAGS scripts used in this file
# seedBagsCompletePoolingViabilityPartialPoolingBurialJAGS.R
# seedBagsPartialPoolingViabilityPartialPoolingBurialJAGS.R

# -------------------------------------------------------------------
# Models for joint estimates of year 1 below ground rates
# Seed survival, germination, and viability
# Models use log-odds parameterization
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
library(magrittr)
library(tidybayes)
library(parallel)
library(stringr)

set.seed(10)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Read data
# -------------------------------------------------------------------
# -------------------------------------------------------------------

#data = readRDS("~/Dropbox/dataLibrary/workflow/data/belowgroundDataAgeOneTwo.RDS")

data = readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/viabilityData.RDS")

#data$intactOctDup = data$intactOct

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.chain = 3
n.adapt = 3000
n.update = 5000
n.iterations = 10000
n.thin = 1

dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Complete pooling of germination and viability trials
# Partial pooling of seed burial experiment (site level)
# -------------------------------------------------------------------
# -------------------------------------------------------------------

initsMu0 <- function(samps = data$n_siteBags){
  rnorm(n = samps, mean = 0, sd = 1)
}

initsSigma0 <- function(samps = data$n_siteBags){
  extraDistr::rhnorm(n = samps, sigma = 1/2)
}

initsSigma <- function(rows = data$n_siteBags, cols = data$n_yearBags){
  matrix(extraDistr::rhnorm(n = rows*cols, sigma = 1/2), rows, cols)
}

# set inits for JAGS
inits <- list()
for(i in 1:3){
  inits[[i]] <- list(initsMu0(), initsMu0(),
                     initsSigma0(), initsSigma0(), 
                     initsSigma(cols=data$n_yearBags), initsSigma(cols=data$n_yearBags) )
  
  names(inits[[i]]) = c(paste(rep("mu0",2),c("g","v"),sep="_"),
                        paste(rep("sigma0",2),c("g","v"),sep="_"),
                        paste(rep("sigma",2),c("g","v"),sep="_"))
  
}


# # tuning (n.adapt)
jm = jags.model(paste0(dir,"jagsModelViability.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

parsToMonitor = c(paste(rep("mu0",2),c("g","v"),sep="_"),
                  paste(rep("mu",2),c("g","v"),sep="_"),
                  paste(rep("sigma0",2),c("g","v"),sep="_"),
                  paste(rep("sigma",2), c("g","v"),sep="_"),
                  "p_g","p_v",
                  "nu_1",
                  "p0_g","p0_v",
                  "nu0_1",
                  "theta_g","theta_v",
                  "nu_bag"
                  )

# chain (n.iter)
samples.rjags = coda.samples(jm,
                             variable.names = c(parsToMonitor),
                             n.iter = n.iterations, thin = n.thin)

fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/posteriors/modelChecking/")
 dir.create(file.path(fileDirectory), showWarnings = FALSE)
 
 saveRDS(samples.rjags,file=paste0(fileDirectory,"viabilitySamples.RDS"))

