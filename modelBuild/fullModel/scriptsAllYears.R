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

data = readRDS("~/Dropbox/dataLibrary/workflow/data/belowgroundDataAgeOneTwo.RDS")

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
n.thin = 10

dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/jagsScriptsSeedBags/")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Complete pooling of germination and viability trials
# Partial pooling of seed burial experiment (site level)
# -------------------------------------------------------------------
# -------------------------------------------------------------------

initsMu0 <- function(samps = data$n_siteBags){
  rnorm(n = samps, mean = 0, sd =1)
}

initsSigma0 <- function(samps = data$n_siteBags){
  extraDistr::rhnorm(n = samps, sigma = 1)
}

initsSigma <- function(rows = data$n_siteBags, cols = data$n_yearBags){
  matrix(extraDistr::rhnorm(n = rows*cols, sigma = 1), rows, cols)
}

# set inits for JAGS
inits <- list()
for(i in 1:1){
  inits[[i]] <- list(initsMu0(), initsMu0(), initsMu0(), initsMu0(), initsMu0(),
                     initsSigma0(), initsSigma0(), initsSigma0(), initsSigma0(), initsSigma0(),
                     initsSigma(cols=data$n_yearBags2), initsSigma(cols=data$n_yearBags2), initsSigma(cols=data$n_yearBags2), initsSigma(cols=data$n_yearBags2), initsSigma(cols=data$n_yearBags2))

  names(inits[[i]]) = c(paste(rep("mu0",5),c(4:6,"g2","v2"),sep="_"),
                        paste(rep("sigma0",5),c(4:6,"g2","v2"),sep="_"),
                        paste(rep("sigma",5), c(4:6,"g2","v2"),sep="_"))
}



# # tuning (n.adapt)
jm = jags.model(paste0(dir,"hierarchicalLogitCenteredAgeOneTwo.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

parsToMonitor = c(paste(rep("mu0",5),c(1:3,"g","v"),sep="_"),
  paste(rep("mu0",5),c(4:6,"g2","v2"),sep="_"),
  paste(rep("mu",5),c(1:3,"g","v"),sep="_"),
  paste(rep("mu",5),c(4:6,"g2","v2"),sep="_"),
  paste(rep("sigma0",5),c(1:3,"g","v"),sep="_"),
  paste(rep("sigma0",5),c(4:6,"g2","v2"),sep="_"),
  paste(rep("sigma",5), c(1:3,"g","v"),sep="_"),
  paste(rep("sigma",5), c(4:6,"g2","v2"),sep="_"),
  "s1","g1","s2")

# parsToMonitor_1 = c("mu0_1","sigma0_1","mu_1","sigma_1")
# parsToMonitor_2 = c("mu0_2","sigma0_2","mu_2","sigma_2")
# parsToMonitor_3 = c("mu0_3","sigma0_3","mu_3","sigma_3")
# parsToMonitor_g = c("mu0_g","sigma0_g","mu_g","sigma_g")
# parsToMonitor_v = c("mu0_v","sigma0_v","mu_v","sigma_v") 



# chain (n.iter)
samples.rjags = coda.samples(jm,
                             variable.names = c(parsToMonitor),
                             n.iter = n.iterations, thin = n.thin)

fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/workflow/samples/")
dir.create(file.path(fileDirectory), showWarnings = FALSE)

saveRDS(samples.rjags,file=paste0(fileDirectory,"belowgroundSamplesAllYears.RDS"))

#MCMCsummary(samples.rjags,params="sigma0_v")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Parameters to monitor
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# 
# parsToMonitor_1 = c("theta_1","mu0_1","sigma0_1","mu_1","sigma_1","p_1")
# parsToMonitor_2 = c("theta_2","mu0_2","sigma0_2","mu_2","sigma_2","p_2")
# parsToMonitor_3 = c("theta_3","mu0_3","sigma0_3","mu_3","sigma_3","p_3")
# parsToMonitor_g = c("theta_g","mu0_g","sigma0_g","mu_g","sigma_g","p_g")
# parsToMonitor_v = c("theta_v","mu0_v","sigma0_v","mu_v","sigma_v","p_v")
# parsToMonitor_deriv = c("nu_1","s1","g1","s2")
# parsToMonitor_deriv2 = c("nu0_1","s1.0","g1.0","s2.0")
# 
# pars = c(parsToMonitor_1,parsToMonitor_2,parsToMonitor_3,parsToMonitor_g,parsToMonitor_v,parsToMonitor_deriv,parsToMonitor_deriv2)
# 
# 
# # -------------------------------------------------------------------
# # -------------------------------------------------------------------
# # Set up JAGS for parallel run
# # -------------------------------------------------------------------
# # -------------------------------------------------------------------
# 
# cl <- makeCluster(3)
# 
# myWorkers <- NA
# for(i in 1:3) myWorkers[i] <- stringr::word(capture.output(cl[[i]]), -1)
# 
# initsP = function() {temp = list(initsMu0(), initsMu0(), initsMu0(), initsMu0(), initsMu0(),
#                           initsSigma0(), initsSigma0(), initsSigma0(), initsSigma0(), initsSigma0(),
#                           initsSigma(), initsSigma(), initsSigma(), initsSigma(), initsSigma(),
#                           "base::Mersenne-Twister", runif(1, 1, 2000))
#   
#   names(temp) = c(paste(rep("mu0",5),c(1:3,"g","v"),sep="_"), 
#                         paste(rep("sigma0",5),c(1:3,"g","v"),sep="_"),
#                         paste(rep("sigma",5), c(1:3,"g","v"),sep="_"),
#                         ".RNG.name",".RNG.seed")
#   return(temp)
# }
# initsP = list(initsP(),initsP(),initsP())
# 
# path=paste0(dir,"hierarchicalLogitCentered.R")
# 
# parallel::clusterExport(cl, c("myWorkers","data", "initsP", "n.adapt", "n.update",
#                               "n.iterations","path"))
# 
# 
# out <- clusterEvalQ(cl, {
#   library(rjags)
#   jm = jags.model(file=path, data = data, n.chains = 1,
#                   n.adapt = n.adapt, inits = initsP[[which(myWorkers==Sys.getpid())]])
#   update(jm, n.iter = n.update)
#   zm = coda.samples(jm, variable.names = c("nu0_1","s1.0","g1.0","s2.0"),
#                     n.iter = n.iterations, thin = 1)
#   return(as.mcmc(zm))
# })
# stopCluster(cl)
# zmP = mcmc.list(out)
# 
# saveRDS(samples.rjags,file=paste0(fileDirectory,"seedBurialSamples.rds"))
