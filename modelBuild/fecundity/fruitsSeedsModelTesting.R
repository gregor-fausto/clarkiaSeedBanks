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
  dplyr::rename(y_tfe = countFruitNumberPerPlant) %>%
  dplyr::select(site,year,y_tfe)

countFruitsPerPlantAllPlots$year <- as.character(countFruitsPerPlantAllPlots$year)


countUndamagedDamagedFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countUndamagedDamagedFruitsPerPlantAllPlots.RDS")

countUndamagedDamagedFruitsPerPlantAllPlots <- countUndamagedDamagedFruitsPerPlantAllPlots %>%
  dplyr::rename(y_und = countUndamagedFruitNumberPerPlant) %>%
  dplyr::rename(y_dam = countDamagedFruitNumberPerPlant) %>%
  dplyr::rename(site2 = site) %>%
  dplyr::rename(year2 = year) %>%
  dplyr::select(site2,year2,y_und,y_dam)

countUndamagedDamagedFruitsPerPlantAllPlots$year2 <- as.character(countUndamagedDamagedFruitsPerPlantAllPlots$year2)



countSeedPerFruit <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countSeedPerFruit.RDS")

countSeedPerUndamagedFruit <- countSeedPerFruit %>%
  dplyr::filter(demography==1) %>%
  dplyr::filter(damaged==0) %>%
  dplyr::rename(site3 = site) %>%
  dplyr::rename(year3 = year) %>%
  dplyr::select(site3,year3,sdno)

countSeedPerUndamagedFruit$year3 <- as.character(countSeedPerUndamagedFruit$year3)

countSeedPerDamagedFruit <- countSeedPerFruit %>%
  dplyr::filter(demography==1) %>%
  dplyr::filter(damaged==1) %>%
  dplyr::rename(site4 = site) %>%
  dplyr::rename(year4 = year) %>%
  dplyr::rename(sdno_dam = sdno) %>%
  dplyr::select(site4,year4,sdno_dam)

countSeedPerDamagedFruit$year4 <- as.character(countSeedPerDamagedFruit$year4)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------
# -------------------------------------------------------------------
data <- tidybayes::compose_data(countFruitsPerPlantAllPlots,
                                countUndamagedDamagedFruitsPerPlantAllPlots,
                                countSeedPerUndamagedFruit,
                                countSeedPerDamagedFruit)

data$n = dim(countFruitsPerPlantAllPlots)[1]
data$n2 = dim(countUndamagedDamagedFruitsPerPlantAllPlots)[1]
data$n3 = dim(countSeedPerUndamagedFruit)[1]
data$n4 = dim(countSeedPerDamagedFruit)[1]

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
n.iterations = 5000
n.thin = 10

dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/jagsScriptsFecundity/")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Complete pooling of germination and viability trials
# Partial pooling of seed burial experiment (site level)
# -------------------------------------------------------------------
# -------------------------------------------------------------------

initsMu0 <- function(samps = data$n_site){
  rnorm(n = samps, mean = 0, sd = 1)
}

initsSigma0 <- function(samps = data$n_site){
  extraDistr::rhnorm(n = samps, sigma = 1)
}

initsR <- function(rows = data$n_site, cols = data$n_year){
  samps = rows*cols
  matrix(extraDistr::rhnorm(n = samps, sigma = 1), rows, cols)
}

# set inits for JAGS
inits <- list()
for(i in 1:3){
  inits[[i]] <- list(initsMu0(), initsMu0(), initsMu0(), initsMu0(), initsMu0(),
                     initsSigma0(), initsSigma0(), initsSigma0(), initsSigma0(), initsSigma0(),
                     initsR(rows = data$n_site,cols=data$n_year), 
                     initsR(rows = data$n_site,cols=data$n_year2), 
                     initsR(rows = data$n_site,cols=data$n_year2), 
                     initsR(rows = data$n_site,cols=data$n_year3), 
                     initsR(rows = data$n_site,cols=data$n_year4) )
  
  names(inits[[i]]) = c(paste(rep("mu0",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                        paste(rep("sigma0",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                        paste(rep("r",5),c("tfe.inv","und.inv","dam.inv","seeds.inv","dam_seeds.inv"),sep="_"))
  
}

# # Call to JAGS
# 
# # tuning (n.adapt)
jm = jags.model(paste0(dir,"seedsAllJagsTesting2.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)


parsToMonitor = c(paste(rep("mu0",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                  paste(rep("sigma0",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                  paste(rep("r",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                  paste(rep("r",5),c("tfe.inv","und.inv","dam.inv","seeds.inv","dam_seeds.inv"),sep="_"),
                  paste(rep("lambda",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                  "ratio",
                  paste(rep("p0",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                  "tfe","p_und","p_dam","tfe_comp",
                  "y_und_sim")
                                 
# chain (n.iter)
samples.rjags = coda.samples(jm, 
                             variable.names = c(parsToMonitor), 
                             n.iter = n.iterations, thin = n.thin)

fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/")
dir.create(file.path(fileDirectory), showWarnings = FALSE)
# 
saveRDS(samples.rjags,file=paste0(fileDirectory,"fitnessSamples.rds"))
saveRDS(samples.rjags,file=paste0("/Users/Gregor/Dropbox/dataLibrary/posteriors/fitnessSamplesTesting.rds"))
saveRDS(data,file=paste0(fileDirectory,"data.rds"))

saveRDS(countFruitsPerPlantAllPlots,file=paste0(fileDirectory,"countFruitsPerPlantAllPlots.rds"))
saveRDS(countUndamagedDamagedFruitsPerPlantAllPlots,file=paste0(fileDirectory,"countUndamagedDamagedFruitsPerPlantAllPlots.rds"))
saveRDS(countSeedPerUndamagedFruit,file=paste0(fileDirectory,"countSeedPerUndamagedFruit.rds"))
saveRDS(countSeedPerDamagedFruit,file=paste0(fileDirectory,"countSeedPerDamagedFruit.rds"))


