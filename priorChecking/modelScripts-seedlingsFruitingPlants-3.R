# list of JAGS scripts used in this file

# -------------------------------------------------------------------
# Models for ...
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
# REMOVE DATA WHERE NUMBER OF FRUITING PLANTS IS GREATER THAN NUMBER OF SEEDLINGS
# -------------------------------------------------------------------
censusSeedlingsFruitingPlants <- censusSeedlingsFruitingPlants %>% 
  dplyr::filter(!(fruitplNumber>seedlingNumber)) %>%
  # FOR PRIOR PREDICTIVE, SELECT ONE SITE
  dplyr::filter(site=="BG")

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
dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/priorChecking/jagsScripts/")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Describe pooling structure
# -------------------------------------------------------------------
# -------------------------------------------------------------------

initsMu0 <- function(samps = data$n_site){
  rnorm(n = samps, mean = 0, sd = 1)
}

initsSigma0 <- function(samps = data$n_site){
  extraDistr::rhnorm(n = samps, sigma = 1)
  # runif(n = samps, min = 0, max = 1)
}

initsSigma <- function(rows = data$n_site, cols = data$n_year){
  matrix(extraDistr::rhnorm(n = rows*cols, sigma = 1), rows, cols)
 # matrix(runif(n=rows*cols, min = 0, max = 1), rows, cols)
}

# set inits for JAGS
inits <- list()
for(i in 1:3){
  inits[[i]] <- list(initsMu0(), initsSigma0(), initsSigma() )
  
  names(inits[[i]]) = c("mu0","sigma0","sigma")
  
}

# set inits for JAGS
# inits = list(list(mu0_1 = rep(0,data$n_site), sigma0_1 = rep(.5,data$n_site),
#                   sigma_1 = matrix(rep(.5,data$n_site*data$n_year),nrow=data$n_site,ncol=data$n_year)),
#              list(mu0_1 = rep(-1,data$n_site), sigma0_1 = rep(1,data$n_site),
#                   sigma_1 = matrix(rep(1,data$n_site*data$n_year),nrow=data$n_site,ncol=data$n_year)),
#              list(mu0_1 = rep(1,data$n_site), sigma0_1 = rep(1.25,data$n_site),
#                   sigma_1 = matrix(rep(1.25,data$n_site*data$n_year),nrow=data$n_site,ncol=data$n_year)))

# # Call to JAGS
# 
# # tuning (n.adapt)
jm = jags.model(paste0(dir,"binomialLikelihood-logitLink-normalHierarchical-3.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

parsToMonitor = c("theta","mu0","sigma0","tau0","mu","sigma","tau")
parsToMonitor_predprior = c("y_prior")
#parsToMonitor_deriv = c("p0","p")
#parsToMonitor_fit = c("fruitplNumber.sim","p.value.mean","p.value.sd")

# chain (n.iter)
samples.rjags = coda.samples(jm, 
                             variable.names = c(parsToMonitor, parsToMonitor_predprior), 
                             n.iter = n.iterations, thin = n.thin)

# fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedlingSurvival/")
# dir.create(file.path(fileDirectory), showWarnings = FALSE)
# # 
# saveRDS(samples.rjags,file=paste0(fileDirectory,"seedSurvivalSamples.rds"))
# saveRDS(samples.rjags,file=paste0("/Users/Gregor/Dropbox/dataLibrary/posteriors/seedSurvivalSamples.rds"))
# saveRDS(data,file=paste0(fileDirectory,"data.rds"))
# saveRDS(censusSeedlingsFruitingPlants,file=paste0(fileDirectory,"censusSeedlingsFruitingPlants.rds"))

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Evaluate priors
# -------------------------------------------------------------------
# -------------------------------------------------------------------

hist.chain=function(x){
  name<-colnames(x)
  tmp = hist(x,breaks=100,plot=FALSE)
  tmp$xname = name
  plot(tmp,freq=FALSE)
}

## Priors on hyperparameters
par(mfrow=c(1,3))
hist.chain(MCMCchains(samples.rjags,params="mu0"))
hist.chain(MCMCchains(samples.rjags,params="sigma0"))
hist.chain(MCMCchains(samples.rjags,params="sigma")[,1,drop=FALSE])

## Priors on transformed parameters
n_params=dim(MCMCchains(samples.rjags,params=c("theta")))[2]
par(mfrow=c(1,3))
hist.chain(MCMCchains(samples.rjags,params=c("theta"))[,sample(1:n_params,1),drop=FALSE])
hist.chain(MCMCchains(samples.rjags,params=c("theta"))[,sample(1:n_params,1),drop=FALSE])
hist.chain(MCMCchains(samples.rjags,params=c("theta"))[,sample(1:n_params,1),drop=FALSE])

## Prior predictive
y_prior=MCMCchains(samples.rjags,params="y_prior")

par(mfrow=c(1,2))
plot(data$fruitplNumber,y_prior[sample(dim(y_prior)[1],1),],
     type='n',ylim=c(0,200),
     xlab="Observed data",ylab="Simulated data")
for(i in 1:100){
  points(data$fruitplNumber,y_prior[sample(dim(y_prior)[1],1),],
         pch=16,cex=.5)
}
abline(a=0,b=1)

# easier to see the full range from the prior in a proportion plot
# certain values are only observed in certain combinations because of denominator (number of trials)
plot(data$fruitplNumber/data$seedlingNumber,
     y_prior[sample(dim(y_prior)[1],1),]/data$seedlingNumber,
     type='n',xlim=c(0,1),ylim=c(0,1),
     xlab="Observed data",ylab="Simulated data")
for(i in 1:100){
  points(data$fruitplNumber/data$seedlingNumber,
         y_prior[sample(dim(y_prior)[1],1),]/data$seedlingNumber,
         pch=16,cex=.25)
}

# only 3 plots with >10 seedlings have survival>.8
censusSeedlingsFruitingPlants %>%
  dplyr::mutate(p = fruitplNumber/seedlingNumber) %>%
  dplyr::filter(p>0.8)

## SUMMARY
MCMCsummary(samples.rjags,params=c("mu0","sigma0","tau0"))
MCMCsummary(samples.rjags,params=c("mu","sigma","tau"))



