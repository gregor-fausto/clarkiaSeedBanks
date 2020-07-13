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

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
countFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countFruitsPerPlantAllPlots.RDS")

countFruitsPerPlantAllPlots <- countFruitsPerPlantAllPlots %>%
  dplyr::rename(countFruitsPerPlant = countFruitNumberPerPlant) %>%
  dplyr::select(site,year,countFruitsPerPlant)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------
# -------------------------------------------------------------------
countFruitsPerPlantAllPlots$year <- as.character(countFruitsPerPlantAllPlots$year)
data <- tidybayes::compose_data(countFruitsPerPlantAllPlots)

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
jm = jags.model(paste0(dir,"fecJags.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

parsToMonitor = c("mu0","sigma0","gamma","r","lambda","p","p0")

# chain (n.iter)
samples.rjags = coda.samples(jm, 
                             variable.names = c(parsToMonitor), 
                             n.iter = n.iterations, thin = n.thin)

fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/fruitsPerPlantAllPlots/")
dir.create(file.path(fileDirectory), showWarnings = FALSE)
# 
saveRDS(samples.rjags,file=paste0(fileDirectory,"fruitsPerPlantSamples.rds"))
saveRDS(samples.rjags,file=paste0("/Users/Gregor/Dropbox/dataLibrary/posteriors/fruitsPerPlantSamples.rds"))
saveRDS(data,file=paste0(fileDirectory,"data.rds"))
saveRDS(countFruitsPerPlantAllPlots,file=paste0(fileDirectory,"countFruitsPerPlantAllPlots.rds"))


