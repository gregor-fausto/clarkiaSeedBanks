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
# Import and organize seed bag data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# setwd and read data files
seedBagsData <- readRDS("~/Dropbox/dataLibrary/postProcessingData/seedBagsData.rds")

# -------------------------------------------------------------------
# Clean and organize seed bag data
# -------------------------------------------------------------------
seedBagsData$seedStart<-100

seedBagsData$seedStart <- as.integer(seedBagsData$seedStart)

## FOR NOW REMOVE MISSING DATA AND PROBLEMS
# `MC III 5 13 1 is a problem (91 intact, 17 seedlings)`
# the main issue will be bags that were not recovered in January
# but then recovered in October; there is no way of getting an estimate on how many 
# seeds 'started' those trials

seedBagsData <- seedBagsData %>%
  dplyr::filter(totalJan<=100)

seedBagsData<-subset(seedBagsData,!is.na(seedBagsData$totalJan))
seedBagsData<-subset(seedBagsData,!is.na(seedBagsData$intactOct))
seedBagsData<-subset(seedBagsData,!(seedBagsData$intactJan<seedBagsData$intactOct))

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize viability trial data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# setwd and read data files
viabilityRawData <- readRDS("~/Dropbox/dataLibrary/postProcessingData/viabilityRawData.rds")

viabilityRawData <- viabilityRawData %>% 
  dplyr::select(-c(germPerc,germNot,viabPerc,viabPerc2,condTest))

viabilityRawData$bag <-as.integer(as.numeric(viabilityRawData$bagNo))

# another check
seedBurial <- seedBagsData %>% dplyr::select(site,bagNo,round,age) %>%
  dplyr::mutate(seedBurial = 1)
viabilityTrial <- viabilityRawData %>% dplyr::select(site,bag,round,age) %>%
  dplyr::mutate(viabilityTrial = 1)

joinedDF <- seedBurial %>% 
  dplyr::full_join(viabilityTrial) 

# one row is coded differently so that 
# viabStart=NA and viabStain=NA
# all others have viabStart=NA and viabStain=NA
# recode
viabilityRawData[is.na(viabilityRawData$viabStart),]$viabStart = 0

## data check
viabilityRawData %>% dplyr::filter(germStart - germCount - viabStart<0) 

# filter out rows with problems
# these need to be corrected
viabilityRawData<-viabilityRawData %>% 
  dplyr::filter(germStart - germCount - viabStart >= 0)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------
# -------------------------------------------------------------------

## filter the dataset for testing purposes
filterData<-function(x) {
  x %>%
    dplyr::filter(age==1) 
}

seedBagsData<-filterData(seedBagsData)
viabilityRawData<-filterData(viabilityRawData)

# assign variable that combines site and bag; unique id for each bag
# for each dataset, create that unique identifier again and then
# use that to link it to the reference identifier created above
seedBagsData<-seedBagsData %>%
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="", remove=FALSE) %>%
  tidyr::unite(col='siteBag', c(site,bagNo), sep="", remove=FALSE) %>%
  dplyr::mutate(siteBag = as.factor(siteBag)) 

viabilityRawData<-viabilityRawData %>%
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="", remove=FALSE) %>%
  tidyr::unite(col='siteBag', c(site,bagNo), sep="", remove=FALSE) %>%
  dplyr::mutate(siteBag = as.factor(siteBag)) 

# once each identifier has been created and linked to the reference table
# and the dataset filtered, the dataset needs to be re-indexed
# this may be redundant?

# this line creates a unique id for the subsetted data that is then 
# used to index each of the 2 datasets
# and provides the reference set of bags that were included in the experiment

#https://community.rstudio.com/t/how-to-add-a-counter-to-each-group-in-dplyr/12986/2

referenceTable<-data.frame(id=union(seedBagsData$id, viabilityRawData$id)) %>%
  dplyr::mutate(idNo = 1:length(id)) 

seedBagsData<-seedBagsData %>%
  dplyr::left_join(referenceTable,by="id")

viabilityRawData<-viabilityRawData %>%
  dplyr::left_join(referenceTable,by="id")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------
# -------------------------------------------------------------------

seedBagsData = seedBagsData %>%
  dplyr::mutate(year = as.factor(yearStart)) %>%
  dplyr::select(site,year,totalJan,seedStart,seedlingJan,intactOct) %>%
  dplyr::rename(siteBags = site,
                yearBags = year)

viabilityRawData = viabilityRawData %>%
  dplyr::mutate(year = as.factor(round)) %>%
  dplyr::select(site, year, germStart, germCount, viabStart, viabStain, idNo) %>%
  dplyr::rename(siteViab = site,
                yearViab = year,
                bag = idNo) %>%
  dplyr::mutate(bag = as.factor(bag)) %>%
  dplyr::group_by(siteViab,yearViab,bag) %>%
  # sum observations in each bag; this is ignoring some variation
  dplyr::summarise(germStart = sum(germStart),
                   germCount = sum(germCount),
                   viabStart = sum(viabStart),
                   viabStain = sum(viabStain))

data <- tidybayes::compose_data(seedBagsData,viabilityRawData)

data$n = dim(seedBagsData)[1]
# data$nViab = dim(viabilityExperiment)[1]
# data$n <- NULL

# pass data to list for JAGS
# data = list(
# 
#   # Seed burial experiment, year one
#   y = as.double(seedBagExperiment$totalJan),
#   n = as.double(seedBagExperiment$seedStart),
#   
#   N = nrow(seedBagExperiment)
#   
# )

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
dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/jagsScriptsSeedBags/")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Complete pooling of germination and viability trials
# Partial pooling of seed burial experiment (site level)
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# set inits for JAGS
inits = list(list(mu0_1 = rep(0,data$n_siteBags), sigma0_1 = rep(.5,data$n_siteBags),
                  sigma_1 = matrix(rep(.5,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
                  mu0_2 = rep(0,data$n_siteBags), sigma0_2 = rep(.5,data$n_siteBags),
                  sigma_2 = matrix(rep(.5,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
                  mu0_3 = rep(0,data$n_siteBags), sigma0_3 = rep(.5,data$n_siteBags),
                  sigma_3 = matrix(rep(.5,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
                  mu0_g = rep(0,data$n_siteBags), sigma0_g = rep(.5,data$n_siteBags),
                  sigma_g = matrix(rep(.5,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
                  mu0_v = rep(0,data$n_siteBags), sigma0_v = rep(.5,data$n_siteBags),
                  sigma_v = matrix(rep(.5,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags)),
             list(mu0_1 = rep(-1,data$n_siteBags), sigma0_1 = rep(1,data$n_siteBags),
                  sigma_1 = matrix(rep(1,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
                  mu0_2 = rep(-1,data$n_siteBags), sigma0_2 = rep(1,data$n_siteBags),
                  sigma_2 = matrix(rep(1,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
                  mu0_3 = rep(-1,data$n_siteBags), sigma0_3 = rep(1,data$n_siteBags),
                  sigma_3 = matrix(rep(1,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
                  mu0_g = rep(-1,data$n_siteBags), sigma0_g = rep(1,data$n_siteBags),
                  sigma_g = matrix(rep(1,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
                  mu0_v = rep(-1,data$n_siteBags), sigma0_v = rep(1,data$n_siteBags),
                  sigma_v = matrix(rep(1,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags)),
             list(mu0_1 = rep(1,data$n_siteBags), sigma0_1 = rep(1.25,data$n_siteBags),
                  sigma_1 = matrix(rep(1.25,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
                  mu0_2 = rep(1,data$n_siteBags), sigma0_2 = rep(1.25,data$n_siteBags),
                  sigma_2 = matrix(rep(1.25,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
                  mu0_3 = rep(1,data$n_siteBags), sigma0_3 = rep(1.25,data$n_siteBags),
                  sigma_3 = matrix(rep(1.25,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
                  mu0_g = rep(1,data$n_siteBags), sigma0_g = rep(1.25,data$n_siteBags),
                  sigma_g = matrix(rep(1.25,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
                  mu0_v = rep(1,data$n_siteBags), sigma0_v = rep(1.25,data$n_siteBags),
                  sigma_v = matrix(rep(1.25,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags)))

# # Call to JAGS
# 
# # tuning (n.adapt)
jm = jags.model(paste0(dir,"hierarchicalLogitCentered.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

parsToMonitor_1 = c("theta_1","mu0_1","sigma0_1","mu_1","sigma_1","p_1")
parsToMonitor_2 = c("theta_2","mu0_2","sigma0_2","mu_2","sigma_2","p_2")
parsToMonitor_3 = c("theta_3","mu0_3","sigma0_3","mu_3","sigma_3","p_3")
parsToMonitor_g = c("theta_g","mu0_g","sigma0_g","mu_g","sigma_g","p_g")
parsToMonitor_v = c("theta_v","mu0_v","sigma0_v","mu_v","sigma_v","p_v")
parsToMonitor_deriv = c("nu_1","s1","g1","s2")
parsToMonitor_deriv2 = c("nu0_1","s1.0","g1.0","s2.0")


# chain (n.iter)
samples.rjags = coda.samples(jm, 
                             variable.names = c(parsToMonitor_1,parsToMonitor_2,
                                                parsToMonitor_3,
                                                parsToMonitor_g,parsToMonitor_v,
                                                parsToMonitor_deriv,parsToMonitor_deriv2), 
                             n.iter = n.iterations, thin = n.thin)

fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBurial/")
dir.create(file.path(fileDirectory), showWarnings = FALSE)
# 
saveRDS(samples.rjags,file=paste0(fileDirectory,"seedBurialSamples.rds"))
saveRDS(data,file=paste0(fileDirectory,"data.rds"))
saveRDS(seedBagsData,file=paste0(fileDirectory,"seedBagExperiment.rds"))
saveRDS(viabilityRawData,file=paste0(fileDirectory,"viabilityExperiment.rds"))


MCMCsummary(samples.rjags, params = c("s1","g1"))

