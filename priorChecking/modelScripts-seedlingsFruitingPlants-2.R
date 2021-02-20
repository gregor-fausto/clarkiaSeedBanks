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
  rnorm(n = samps, mean = 0, sd = 1000)
}

initsSigma0 <- function(samps = data$n_site){
  #extraDistr::rhnorm(n = samps, sigma = 1)
  runif(n = samps, min = 0, max = 1)
}

initsSigma <- function(rows = data$n_site, cols = data$n_year){
 # matrix(extraDistr::rhnorm(n = rows*cols, sigma = 1), rows, cols)
  matrix(runif(n=rows*cols, min = 0, max = 1), rows, cols)
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
jm = jags.model(paste0(dir,"binomialLikelihood-logitLink-normalHierarchical-2.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

parsToMonitor = c("theta","mu0","sigma0","mu","sigma")
parsToMonitor_predprior = c("y_pred")
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






MCMCsummary(samples.rjagsWrAg,params=c("mu0_g","mu_g"))

y_pred=MCMCchains(samples.rjagsWrAg,params="y_pred")
theta_g=MCMCchains(samples.rjagsWrAg,params="g")
theta_pred=MCMCchains(samples.rjagsWrAg,params="g_pred")

dim(y_pred)[1]

par(mfrow=c(3,3))
hist(y_pred[sample(dim(y_pred)[1],100),],breaks=100,freq=FALSE)

hist(theta_pred[,1:49],breaks=100,freq=FALSE)
#hist(theta_g[,1:49],breaks=100,freq=FALSE,add=TRUE)

#hist(data$seedlingJan,breaks=25,add=TRUE,col='red',freq=FALSE)

plot(data$seedlingJan,y_pred[sample(dim(y_pred)[1],1),],
     type='n',xlim=c(0,100),ylim=c(0,100),
     xlab="Observed data",ylab="Simulated data")
for(i in 1:100){
  points(data$seedlingJan,y_pred[sample(dim(y_pred)[1],1),],pch=16,cex=.5,
         col=data$yearGermination)
}
abline(a=0,b=1)

#sigma_g=MCMCchains(samples.rjagsWrAg,params="sigma_g")[,1]
#hist(sigma_g,breaks=100,freq=FALSE,ylim=c(0,2),xlim=c(0,5))



#par(mfrow=c(1,3))
plot(density(MCMCchains(samples.rjagsWrAg,params="mu0_g")[,1]),xlim=c(-4,4))
lines(density(MCMCchains(samples.rjagsWrAg,params="mu0_pred")[,1]),col='red')

lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,1]),col='black',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,2]),col='black',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,3]),col='black',lty='dotted')

lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,1]),col='red',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,2]),col='red',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,3]),col='red',lty='dotted')

plot(density(MCMCchains(samples.rjagsWrAg,params="mu0_g")[,2]),xlim=c(-4,4))
lines(density(MCMCchains(samples.rjagsWrAg,params="mu0_pred")[,2]),col='red')

lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,4]),col='black',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,5]),col='black',lty='dotted')

lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,4]),col='red',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,5]),col='red',lty='dotted')

plot(density(MCMCchains(samples.rjagsWrAg,params="mu0_g")[,3]),xlim=c(-4,4))
lines(density(MCMCchains(samples.rjagsWrAg,params="mu0_pred")[,3]),col='red')

lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,6]),col='black',lty='dotted')

lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,6]),col='red',lty='dotted')


## SIGMA

plot(density(MCMCchains(samples.rjagsWrAg,params="sigma0_g")[,1]))
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma0_pred")[,1]),col='red')

lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,1]),col='black',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,2]),col='black',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,3]),col='black',lty='dotted')

lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,1]),col='red',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,2]),col='red',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,3]),col='red',lty='dotted')

plot(density(MCMCchains(samples.rjagsWrAg,params="sigma0_g")[,2]))
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma0_pred")[,2]),col='red')

lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,4]),col='black',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,5]),col='black',lty='dotted')

lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,4]),col='red',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,5]),col='red',lty='dotted')

plot(density(MCMCchains(samples.rjagsWrAg,params="sigma0_g")[,3]))
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma0_pred")[,3]),col='red')

lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,6]),col='black',lty='dotted')

lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,6]),col='red',lty='dotted')


# 
# 
# hist(MCMCchains(samples.rjagsWrAg,params="mu0_g")[,2],breaks=100,freq=FALSE)
# lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,4]),col='red')
# lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,5]),col='red')
# 
# hist(MCMCchains(samples.rjagsWrAg,params="mu0_g")[,3],breaks=100,freq=FALSE)
# lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,6]),col='red')


#sigma0_g=MCMCchains(samples.rjagsWrAg,params="sigma0_g")[,1]

#par(mfrow=c(1,3))
# hist(MCMCchains(samples.rjagsWrAg,params="sigma0_g")[,1],breaks=100,freq=FALSE)
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,1]),col='red')
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,2]),col='red')
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,3]),col='red')
# 
# hist(MCMCchains(samples.rjagsWrAg,params="sigma0_g")[,2],breaks=100,freq=FALSE)
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,4]),col='red')
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,5]),col='red')
# 
# hist(MCMCchains(samples.rjagsWrAg,params="sigma0_g")[,3],breaks=100,freq=FALSE)
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,6]),col='red')
# 

## SUMMARY
MCMCsummary(samples.rjagsWrAg,params=c("tau0_g","tau_g"))
MCMCsummary(samples.rjagsWrAg,params=c("sigma0_g","sigma_g"))
MCMCsummary(samples.rjagsWrAg,params=c("mu0_g","mu_g"))


y_pred=MCMCchains(samples.rjagsWrAg,params="y_pred")
y_sim=MCMCchains(samples.rjagsWrAg,params="y_sim")


par(mfrow=c(1,2))
plot(data$seedlingJan,y_pred[sample(dim(y_pred)[1],1),],
     type='n',xlim=c(0,100),ylim=c(0,100),
     xlab="Observed data",ylab="Simulated data")
for(i in 1:100){
  points(data$seedlingJan,y_pred[sample(dim(y_pred)[1],1),],pch=16,cex=.5,
         col=data$yearGermination)
}
abline(a=0,b=1)

plot(data$seedlingJan,y_sim[sample(dim(y_sim)[1],1),],
     type='n',xlim=c(0,100),ylim=c(0,100),
     xlab="Observed data",ylab="Simulated data")
for(i in 1:100){
  points(data$seedlingJan,y_sim[sample(dim(y_sim)[1],1),],pch=16,cex=.5,
         col=data$yearGermination)
}
abline(a=0,b=1)
