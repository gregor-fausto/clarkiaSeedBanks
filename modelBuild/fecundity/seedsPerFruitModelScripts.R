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

set.seed(10)
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
countSeedPerFruit <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countSeedPerFruit.RDS")

countSeedPerFruit <- countSeedPerFruit %>%
  dplyr::filter(damaged==0) %>%
  dplyr::filter(demography==1) %>%
  dplyr::select(site,year,sdno)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------
# -------------------------------------------------------------------
countSeedPerFruit$year <- as.character(countSeedPerFruit$year)
data <- tidybayes::compose_data(countSeedPerFruit)


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
n.iterations = 10000
n.thin = 10

set.seed(10)
dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/jagsScriptsFecundity/")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Complete pooling of germination and viability trials
# Partial pooling of seed burial experiment (site level)
# -------------------------------------------------------------------
# -------------------------------------------------------------------


inits = list(
  list( 
    mu0 = rep(-1, data$n_site), 
    sigma0 = rep(.5, data$n_site), 
    r = matrix(10, nrow=data$n_site,ncol=data$n_year)),
  list( 
    mu0 = rep(0, data$n_site), 
    sigma0 = rep(1, data$n_site), 
    r = matrix(10, nrow=data$n_site,ncol=data$n_year)),
  list( 
    mu0 = rep(1, data$n_site), 
    sigma0 = rep(1.25, data$n_site), 
    r = matrix(10, nrow=data$n_site,ncol=data$n_year))
)


# # Call to JAGS
# 
# # tuning (n.adapt)
jm = jags.model(paste0(dir,"seedsJags.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

parsToMonitor = c("mu0","sigma0","gamma","r","lambda")

# chain (n.iter)
samples.rjags = coda.samples(jm, 
                             variable.names = c(parsToMonitor), 
                             n.iter = n.iterations, thin = n.thin)

fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedsPerFruit/")
dir.create(file.path(fileDirectory), showWarnings = FALSE)
# 
saveRDS(samples.rjags,file=paste0(fileDirectory,"seedsPerFruitSamples.rds"))
saveRDS(data,file=paste0(fileDirectory,"data.rds"))
saveRDS(countSeedPerFruit,file=paste0(fileDirectory,"countSeedPerFruit.rds"))


