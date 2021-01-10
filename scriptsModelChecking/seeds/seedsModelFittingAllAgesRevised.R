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

initsMu0 <- function(samps = data$n_siteBags, rows = data$n_siteBags, cols = data$n_ageBags){
  matrix(rnorm(n = samps, mean = 0, sd = 1), rows, cols)
}

initsSigma0 <- function(samps = data$n_siteBags, rows = data$n_siteBags, cols = data$n_ageBags){
  matrix(extraDistr::rhnorm(n = samps, sigma = 1/2), rows, cols)
}

initsSigma <- function(rows = data$n_siteBags, cols=6){
  matrix(extraDistr::rhnorm(n = rows*cols, sigma = 1/2), rows, cols)
}

# set inits for JAGS
inits <- list()
for(i in 1:3){
  inits[[i]] <- list(initsMu0(), initsMu0(), initsMu0(), initsMu0(),
                     initsSigma0(), initsSigma0(), initsSigma0(),  initsSigma0(), 
                     initsSigma(cols=6), initsSigma(cols=6), initsSigma(cols=6), initsSigma(cols=6))
  
  names(inits[[i]]) = c(paste(rep("mu0",4),rep(c(1:4),1),sep="_"),
                        paste(rep("sigma0",4),rep(c(1:4),1),sep="_"),
                        paste(rep("sigma",4),rep(c(1:4),1),sep="_"))
  
}


# # tuning (n.adapt)
jm = jags.model(paste0(dir,"jagsModelSeedsAllAgesRevised.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

parsToMonitor = c(paste(rep("mu0",4),c(1:4),sep="_"),
                  paste(rep("mu",4),c(1:4),sep="_"),
                  paste(rep("sigma0",4),c(1:4),sep="_"),
                  paste(rep("sigma",4),c(1:4),sep="_"),
                  "p_1","p_2","p_3","p_4",
                  "p0_1","p0_2","p0_3","p0_4",
                  "theta_1","theta_2","theta_3","theta_4"
)

# chain (n.iter)
samples.rjags = coda.samples(jm,
                             variable.names = c(parsToMonitor),
                             n.iter = n.iterations, thin = n.thin)


fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/posteriors/modelChecking/")
dir.create(file.path(fileDirectory), showWarnings = FALSE)

saveRDS(samples.rjags,file=paste0(fileDirectory,"seedSamplesAllAges.RDS"))

