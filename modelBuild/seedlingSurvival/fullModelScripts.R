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
# setwd and read data files
censusSeedlingsFruitingPlants <- readRDS("~/Dropbox/dataLibrary/postProcessingData/censusSeedlingsFruitingPlants.RDS")

# -------------------------------------------------------------------
# Clean and organize data
# -------------------------------------------------------------------
censusSeedlingsFruitingPlants <- censusSeedlingsFruitingPlants %>% 
  dplyr::filter(!(fruitplNumber>seedlingNumber))

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------
# -------------------------------------------------------------------
censusSeedlingsFruitingPlants$year <- as.character(censusSeedlingsFruitingPlants$year)
data <- tidybayes::compose_data(censusSeedlingsFruitingPlants)

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

set.seed(10)
dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/jagsScriptsSeedlingSurvival/")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Complete pooling of germination and viability trials
# Partial pooling of seed burial experiment (site level)
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# set inits for JAGS
inits = list(list(mu0_1 = rep(0,data$n_site), sigma0_1 = rep(.5,data$n_site),
                  sigma_1 = matrix(rep(.5,data$n_site*data$n_year),nrow=data$n_site,ncol=data$n_year)),
             list(mu0_1 = rep(-1,data$n_site), sigma0_1 = rep(1,data$n_site),
                  sigma_1 = matrix(rep(1,data$n_site*data$n_year),nrow=data$n_site,ncol=data$n_year)),
             list(mu0_1 = rep(1,data$n_site), sigma0_1 = rep(1.25,data$n_site),
                  sigma_1 = matrix(rep(1.25,data$n_site*data$n_year),nrow=data$n_site,ncol=data$n_year)))

# # Call to JAGS
# 
# # tuning (n.adapt)
jm = jags.model(paste0(dir,"hierarchicalLogitCentered.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

parsToMonitor_1 = c("theta_1","mu0_1","sigma0_1","mu_1","sigma_1")
parsToMonitor_deriv = c("p0_1","p_1")

# chain (n.iter)
samples.rjags = coda.samples(jm, 
                             variable.names = c(parsToMonitor_1, parsToMonitor_deriv), 
                             n.iter = n.iterations, thin = n.thin)

fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedlingSurvival/")
dir.create(file.path(fileDirectory), showWarnings = FALSE)
# 
saveRDS(samples.rjags,file=paste0(fileDirectory,"seedSurvivalSamples.rds"))
saveRDS(samples.rjags,file=paste0("/Users/Gregor/Dropbox/dataLibrary/posteriors/seedSurvivalSamples.rds"))
saveRDS(data,file=paste0(fileDirectory,"data.rds"))
saveRDS(censusSeedlingsFruitingPlants,file=paste0(fileDirectory,"censusSeedlingsFruitingPlants.rds"))
