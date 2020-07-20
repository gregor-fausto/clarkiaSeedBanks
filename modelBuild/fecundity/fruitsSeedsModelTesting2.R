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
  dplyr::select(site2,year2,y_und,y_dam) %>%
  dplyr::filter(!is.na(y_und))

countUndamagedDamagedFruitsPerPlantAllPlots$year2 <- as.character(countUndamagedDamagedFruitsPerPlantAllPlots$year2)

# countUndamagedDamagedFruitsPerPlantAllPlots %>%
#   dplyr::filter(is.na(y_und)) %>% 
#   dplyr::select(site2,year2)

#countUndamagedDamagedFruitsPerPlantAllPlots$und_missing=ifelse(is.na(countUndamagedDamagedFruitsPerPlantAllPlots$y_und),0,1)

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
# For testing
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# keep<-data$site2==12
# data$site2<-data$site2[keep]-11
# data$year2<-data$year2[keep]
# data$y_und<-data$y_und[keep]
# data$n2 = length(data$y_und)
# data$n_site2 = 1

# index = 1:length(data$y_und)
# #und_missing = countUndamagedDamagedFruitsPerPlantAllPlots$und_missing[keep]
# present = index*und_missing
# data$present2 = present[present>0]

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.chain = 1
n.adapt = 3000
n.update = 5000
n.iterations = 10000
n.thin = 10

dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/jagsScriptsFecundity/")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Complete pooling of germination and viability trials
# Partial pooling of seed burial experiment (site level)
# -------------------------------------------------------------------
# -------------------------------------------------------------------
initsMu0 <- function(samps = data$n_site2){
  #rnorm(n = samps, mean = 0, sd = 1)
  #truncdist::rtrunc(n = samps,"gamma",a=.001,shape=0.01,rate=0.01)
  rgamma(n = samps, shape = 1, rate = 1)
}

initsSigma0 <- function(samps = data$n_site2){
  #extraDistr::rhnorm(n = samps, sigma = 1)
   truncdist::rtrunc(n = samps,"gamma",a=.001,shape=0.001,rate=0.001)
  # rgamma(n = samps, shape = .001, rate = .001)
  
}

initsR <- function(rows = data$n_site2, cols = data$n_year2){
  samps = rows*cols
  #matrix(truncdist::rtrunc(n = samps,"gamma",a=.001,shape=0.001,rate=0.001), rows, cols)
  #matrix(runif(n = samps,0,1), rows, cols)
  matrix(truncdist::rtrunc(n = samps,"gamma",a=.001,shape=0.001,rate=0.001), rows, cols)
}

# set inits for JAGS
inits <- list()
for(i in 1:3){
  inits[[i]] <- list(initsMu0(), initsMu0(), initsMu0(), initsMu0(), initsMu0(),
                     initsSigma0(), initsSigma0(), initsSigma0(), initsSigma0(), initsSigma0(),
                     initsR(rows = data$n_site2,cols=data$n_year), 
                     initsR(rows = data$n_site2,cols=data$n_year2), 
                     initsR(rows = data$n_site2,cols=data$n_year2), 
                     initsR(rows = data$n_site2,cols=data$n_year3), 
                     initsR(rows = data$n_site2,cols=data$n_year4) )
  
  names(inits[[i]]) = c(paste(rep("nu",3),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                        paste(rep("tau0",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                        paste(rep("tau",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"))

}

# # Call to JAGS
# 
# # tuning (n.adapt)
jm = jags.model(paste0(dir,"seedsAllJagsTesting2.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)


parsToMonitor = c(paste(rep("nu",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                  paste(rep("tau0",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                  paste(rep("sigma0",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                  paste(rep("tau",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                  paste(rep("sigma",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                  paste(rep("gamma",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                  paste(rep("lambda",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                  "ratio","nu_und","g_und","mu_und",
                  paste(rep("p0",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                  paste(rep("p",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                  "tfe","p_und","p_dam","tfe_comp",
                  "y_und_sim","mu_p","mu_py","beta_p","mean_p","beta_py","mean_py")
                                 
# chain (n.iter)
samples.rjags = coda.samples(jm, 
                             variable.names = c(parsToMonitor), 
                             n.iter = n.iterations, thin = n.thin)

fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/")
dir.create(file.path(fileDirectory), showWarnings = FALSE)
# 
#saveRDS(samples.rjags,file=paste0(fileDirectory,"fitnessSamples.rds"))
saveRDS(samples.rjags,file=paste0("/Users/Gregor/Dropbox/dataLibrary/posteriors/fitnessSamplesTesting2.rds"))
saveRDS(data,file=paste0(fileDirectory,"data.rds"))

# saveRDS(countFruitsPerPlantAllPlots,file=paste0(fileDirectory,"countFruitsPerPlantAllPlots.rds"))
# saveRDS(countUndamagedDamagedFruitsPerPlantAllPlots,file=paste0(fileDirectory,"countUndamagedDamagedFruitsPerPlantAllPlots.rds"))
# saveRDS(countSeedPerUndamagedFruit,file=paste0(fileDirectory,"countSeedPerUndamagedFruit.rds"))
# saveRDS(countSeedPerDamagedFruit,file=paste0(fileDirectory,"countSeedPerDamagedFruit.rds"))


