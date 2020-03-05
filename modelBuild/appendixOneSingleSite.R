################################################################################
# Load packages
################################################################################
library(janitor)
library(tidyverse)
library(tidybayes)
library(rjags)
library(magrittr) # for %<>% with recover_types
library(dplyr)
library(MCMCvis)

#################################################################################
# get paths for JAGS script and figure directories
#################################################################################
dirJagsScripts = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/appendixOneScripts/")

# site = "TEST"
# year = as.character(rep(1:30,each=1))
# n = rep(100, length(year))
# # p = rep(seq(.275,.725,by=.05),each=100)
# p = rep(0.5, length(year))
# y = rbinom(length(year), size = n, prob = p)
# 
# simData =  data.frame( site, year, seedlingNumber = n, fruitingPlantNumber = y )
# 
# simDataBayesJags <- tidybayes::compose_data(simData)

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/"
simFiles <- paste0(directory,list.files(directory))

simData <- readRDS(simFiles[[16]])
# 
# ggplot(simData %>% dplyr::filter(site=="EC")) +
#   geom_histogram(aes(fruitingPlantNumber/seedlingNumber)) +
#   facet_wrap(~year) + theme_bw()

selectedSite<-sample(unique(simData$site),1)
simData=simData %>%
  dplyr::filter(site==selectedSite) %>%
  tidyr::unite("plot",c(transect:plot))

simData$year <- factor(simData$year)

simDataBayesJags <- tidybayes::compose_data(simData)

saveRDS(simData, file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/simDataBayes.rds")

################################################################################
################################################################################
# JAGS fit for model with binomial likelihood, beta prior, complete pooling
################################################################################
################################################################################

################################################################################
# Initial values
################################################################################

nSites = simDataBayesJags$n_site
nYears = simDataBayesJags$n_year
nPlot = simDataBayesJags$n_plot

inits = list(list(p = rep(.1,nYears)),
             list(p = rep(.5,nYears)),
             list(p = rep(.9,nYears))) 

################################################################################
# Set parameters to monitor
################################################################################

parsToMonitor = c("p","fruitingPlantNumberSim")

################################################################################
# Set JAGS parameters
################################################################################

# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 5000
n.update = 5000
n.iterations = 100000
n.thin = 10

################################################################################
# Fit model 
################################################################################

# set random seed
set.seed(2)

# tuning (n.adapt)
jm = jags.model(paste0(dirJagsScripts,"survivorshipModel-binLik-betaPrior-completePooling-singleSite.R"), 
                data = simDataBayesJags, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iterations = n.update)

# chain (n.iter)
samplesBinLinkBetaPriorComplete = coda.samples(jm, variable.names = c(parsToMonitor), n.iter = n.iterations, thin = n.thin)

saveRDS(samplesBinLinkBetaPriorComplete,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/singleBinLinkBetaPriorComplete.rds")

# ################################################################################
# ################################################################################
# # JAGS fit for model with binomial likelihood, beta prior, partial pooling
# # parameterized via mode
# ################################################################################
# ################################################################################
# 
# ################################################################################
# # Initial values
# ################################################################################
# 
# nSites = simDataBayesJags$n_site
# nYears = simDataBayesJags$n_year
# nPlot = simDataBayesJags$n_plot
# 
# kappaInit = 100;
# 
# inits = list( list( #theta = matrix(rep(.1,nSites*nYears),nrow=nSites,ncol=nYears),
#   omega=rep(.5,nYears) ,
#   kappaMinusTwo=rep(98,nYears) ),
#   list( #theta = matrix(rep(.5,nSites*nYears),nrow=nSites,ncol=nYears),
#     omega=rep(.5,nYears) ,
#     kappaMinusTwo=rep(98,nYears) ),
#   list( #theta = matrix(rep(.9,nSites*nYears),nrow=nSites,ncol=nYears),
#     omega=rep(.5,nYears) ,
#     kappaMinusTwo=rep(98,nYears) )
# )
# 
# ################################################################################
# # Set parameters to monitor
# ################################################################################
# 
# parsToMonitor = c("theta","omega","kappa","p")
# 
# ################################################################################
# # Set JAGS parameters
# ################################################################################
# 
# # scalars that specify the 
# # number of iterations in the chain for adaptation
# # number of iterations for burn-in
# # number of samples in the final chain
# n.adapt = 500
# n.update = 5000
# n.iterations = 100000
# n.thin = 1
# 
# ################################################################################
# # Fit model 
# ################################################################################
# 
# # set random seed
# set.seed(2)
# 
# # tuning (n.adapt)
# jm = jags.model(paste0(dirJagsScripts,"survivorshipModel-binLik-betaPriorMode-partialPooling-singleSite.R"), 
#                 data = simDataBayesJags, inits = inits,
#                 n.chains = length(inits), n.adapt = n.adapt)
# 
# # burn-in (n.update)
# update(jm, n.iterations = n.update)
# 
# # chain (n.iter)
# samplesBinLinkBetaPriorPartialMode = coda.samples(jm, variable.names = c(parsToMonitor), 
#                                                   n.iter = n.iterations, thin = n.thin)
# 
# saveRDS(samplesBinLinkBetaPriorPartialMode,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/singleSiteBinLinkBetaPriorPartialMode.rds")

################################################################################
################################################################################
# JAGS fit for model with binomial likelihood, beta prior, partial pooling
# parameterized via mean
################################################################################
################################################################################

################################################################################
# Initial values
################################################################################

nSites = simDataBayesJags$n_site
nYears = simDataBayesJags$n_year
nPlot = simDataBayesJags$n_plot

inits = list( list( phi=rep(.1,nYears) ,
                    kappa=rep(1.1,nYears) ),
              list( phi=rep(.5,nYears) ,
                    kappa=rep(1.5,nYears) ),
              list( phi=rep(.9,nYears) ,
                    kappa=rep(2,nYears) )
)

################################################################################
# Set parameters to monitor
################################################################################

parsToMonitor = c("theta","phi","kappa","p","fruitingPlantNumberSim")

################################################################################
# Set JAGS parameters
################################################################################

# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 5000
n.update = 5000
n.iterations = 100000
n.thin = 10

################################################################################
# Fit model 
################################################################################

# set random seed
set.seed(2)

# tuning (n.adapt)
jm = jags.model(paste0(dirJagsScripts,"survivorshipModel-binLik-betaPriorMean-partialPooling-singleSite.R"), 
                data = simDataBayesJags, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iterations = n.update)

# chain (n.iter)
samplesBinLinkBetaPriorPartialMean = coda.samples(jm, variable.names = c(parsToMonitor), 
                                                  n.iter = n.iterations, thin = n.thin)


saveRDS(samplesBinLinkBetaPriorPartialMean,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/singleSiteBinLinkBetaPriorPartialMean.rds")


################################################################################
################################################################################
# JAGS fit for model with binomial likelihood, beta prior, partial pooling
# parameterized via mean, hierarchical to site
# 
################################################################################
################################################################################

################################################################################
# Initial values
################################################################################

library(VGAM)
plot(density(dpareto(seq(0,100,by=0.1),scale=1,shape=1.5)))

plot(density(rpareto(1000,scale=1,shape=1.5)))
plot(density(rpareto(1000,scale=1,shape=1.5)))

nSites = simDataBayesJags$n_site
nYears = simDataBayesJags$n_year
nPlot = simDataBayesJags$n_plot

inits = list( list( phi0=.25,
                    kappa0=10,
                    kappa = matrix(rep(1.1,nYears*nSites),nrow=nSites,ncol=nYears)),
              list( phi0=.5,
                    kappa0=10 ,
                    kappa = matrix(rep(1.5,nYears*nSites),nrow=nSites,ncol=nYears)),
              list( phi0=.75 ,
                    kappa0=10,
                    kappa = matrix(rep(2,nYears*nSites),nrow=nSites,ncol=nYears))
)

################################################################################
# Set parameters to monitor
################################################################################

parsToMonitor = c("phi0","kappa0","phi","kappa","theta","p_site","p","fruitingPlantNumberSim")

################################################################################
# Set JAGS parameters
################################################################################

# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 5000
n.update = 5000
n.iterations = 100000
n.thin = 10

################################################################################
# Fit model 
################################################################################

# set random seed
set.seed(2)

# tuning (n.adapt)
jm = jags.model(paste0(dirJagsScripts,"survivorshipModel-binLik-betaPriorMean-partialPoolingHier-singleSite.R"), 
                data = simDataBayesJags, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iterations = n.update)

# chain (n.iter)
samplesBinLinkBetaPriorPartialMeanHier = coda.samples(jm, variable.names = c(parsToMonitor), 
                                                  n.iter = n.iterations, thin = n.thin)

# MCMCsummary(samplesBinLinkBetaPriorPartialMean,params=c("phi0","kappa0","kappa"))
# MCMCsummary(samplesBinLinkBetaPriorPartialMean,params="phi")
# MCMCsummary(samplesBinLinkBetaPriorPartialMean,params="theta")
# MCMCsummary(samplesBinLinkBetaPriorPartialMean,params="p_site")
# MCMCsummary(samplesBinLinkBetaPriorPartialMean,params="p")

saveRDS(samplesBinLinkBetaPriorPartialMeanHier,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/singleSiteBinLinkBetaPriorPartialHierMean.rds")


################################################################################
################################################################################
# JAGS fit for model with binomial likelihood, beta prior, partial pooling
# parameterized via mean
# with logit link
################################################################################
################################################################################

################################################################################
# Initial values
################################################################################

nSites = simDataBayesJags$n_site
nYears = simDataBayesJags$n_year
nPlot = simDataBayesJags$n_plot

inits = list( list( mu=rep(0,nYears) ,
                    sigma=rep(.5,nYears) ),
              list( mu=rep(-1,nYears) ,
                    sigma=rep(1,nYears) ),
              list( mu=rep(.9,nYears) ,
                    sigma=rep(1.5,nYears) )
)

################################################################################
# Set parameters to monitor
################################################################################

parsToMonitor = c("theta","mu","sigma","alpha","fruitingPlantNumberSim")

################################################################################
# Set JAGS parameters
################################################################################

# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 5000
n.update = 5000
n.iterations = 100000
n.thin = 10

################################################################################
# Fit model 
################################################################################

# set random seed
set.seed(2)

# tuning (n.adapt)
jm = jags.model(paste0(dirJagsScripts,"survivorshipModel-binLik-betaPriorMeanLogit-partialPooling-singleSite.R"), 
                data = simDataBayesJags, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iterations = n.update)

# chain (n.iter)
samplesBinLinkBetaPriorPartialMeanLogit = coda.samples(jm, variable.names = c(parsToMonitor), 
                                                  n.iter = n.iterations, thin = n.thin)
# 
# MCMCsummary(samplesBinLinkBetaPriorPartialMeanLogit,params=c("mu","sigma"))
# 
# MCMCtrace(samplesBinLinkBetaPriorPartialMeanLogit,params=c("mu"))
# MCMCtrace(samplesBinLinkBetaPriorPartialMeanLogit,params=c("sigma"))

 saveRDS(samplesBinLinkBetaPriorPartialMeanLogit,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/singleSiteBinLinkBetaPriorPartialMeanLogit.rds")

 
 ################################################################################
 ################################################################################
 # JAGS fit for model with binomial likelihood, beta prior, partial pooling
 # parameterized via mean
 # with logit link
 # hierarchical
 ################################################################################
 ################################################################################
 
 ################################################################################
 # Initial values
 ################################################################################
 
 
 nSites = simDataBayesJags$n_site
 nYears = simDataBayesJags$n_year
 nPlot = simDataBayesJags$n_plot
 
 inits = list( list( mu0=0 ,
                     sigma0=.5,
                     sigma=matrix(rep(.5,nYears*nSites),nrow=nSites,ncol=nYears)),
               list( mu0=-1 ,
                     sigma0=1,
                     sigma=matrix(rep(1,nYears*nSites),nrow=nSites,ncol=nYears)),
               list(  mu0=1 ,
                      sigma0=1.1,
                      sigma=matrix(rep(1.1,nYears*nSites),nrow=nSites,ncol=nYears)) 
 )
 
 
 
 ################################################################################
 # Set parameters to monitor
 ################################################################################
 
 parsToMonitor = c("theta","mu","sigma","mu0","sigma0","alpha","alpha.std","fruitingPlantNumberSim")
 
 ################################################################################
 # Set JAGS parameters
 ################################################################################
 
 # scalars that specify the 
 # number of iterations in the chain for adaptation
 # number of iterations for burn-in
 # number of samples in the final chain
 n.adapt = 5000
 n.update = 5000
 n.iterations = 100000
 n.thin = 10
 
 ################################################################################
 # Fit model 
 ################################################################################
 
 # set random seed
 set.seed(2)
 
 # tuning (n.adapt)
 jm = jags.model(paste0(dirJagsScripts,"survivorshipModel-binLik-betaPriorMeanHierLogit-partialPooling-singleSite.R"), 
                 data = simDataBayesJags, inits = inits,
                 n.chains = length(inits), n.adapt = n.adapt)
 
 # burn-in (n.update)
 update(jm, n.iterations = n.update)
 
 # chain (n.iter)
 samplesBinLinkBetaPriorPartialMeanHierLogit = coda.samples(jm, variable.names = c(parsToMonitor), 
                                                        n.iter = n.iterations, thin = n.thin)
 # 
 # MCMCsummary(samplesBinLinkBetaPriorPartialMeanLogit,params=c("mu","sigma"))
 # 
 # MCMCtrace(samplesBinLinkBetaPriorPartialMeanLogit,params=c("mu"))
 # MCMCtrace(samplesBinLinkBetaPriorPartialMeanLogit,params=c("sigma"))
 
 saveRDS(samplesBinLinkBetaPriorPartialMeanHierLogit,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/singleSiteBinLinkBetaPriorPartialMeanHierLogit.rds")
 

################################################################################
################################################################################
# JAGS fit for model with binomial likelihood, beta prior, partial pooling
# parameterized via mean
# with logit link, non-centered parameterization
################################################################################
################################################################################

################################################################################
# Initial values
################################################################################

nSites = simDataBayesJags$n_site
nYears = simDataBayesJags$n_year
nPlot = simDataBayesJags$n_plot

inits = list( list( mu=rep(0,nYears) ,
                    sigma=rep(.5,nYears) ),
              list( mu=rep(-1,nYears) ,
                    sigma=rep(1,nYears) ),
              list( mu=rep(.9,nYears) ,
                    sigma=rep(1.5,nYears) )
)

################################################################################
# Set parameters to monitor
################################################################################

parsToMonitor = c("theta","mu","sigma","alpha","alpha.std","fruitingPlantNumberSim")

################################################################################
# Set JAGS parameters
################################################################################

# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 5000
n.update = 5000
n.iterations = 100000
n.thin = 10

################################################################################
# Fit model 
################################################################################

# set random seed
set.seed(2)

# tuning (n.adapt)
jm = jags.model(paste0(dirJagsScripts,"survivorshipModel-binLik-betaPriorMeanLogitNC-partialPooling-singleSite.R"), 
                data = simDataBayesJags, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iterations = n.update)

# chain (n.iter)
samplesBinLinkBetaPriorPartialMeanLogitNC = coda.samples(jm, variable.names = c(parsToMonitor), 
                                                       n.iter = n.iterations, thin = n.thin)

saveRDS(samplesBinLinkBetaPriorPartialMeanLogitNC,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/singleSiteBinLinkBetaPriorPartialMeanLogitNC.rds")


################################################################################
################################################################################
# JAGS fit for model with binomial likelihood, beta prior, partial pooling
# parameterized via mean
# with logit link, non-centered parameterization
# hierarchical
################################################################################
################################################################################

################################################################################
# Initial values
################################################################################

nSites = simDataBayesJags$n_site
nYears = simDataBayesJags$n_year
nPlot = simDataBayesJags$n_plot

inits = list( list( mu0=0 ,
                    sigma0=.5,
                    sigma=matrix(rep(.5,nYears*nSites),nrow=nSites,ncol=nYears)),
              list( mu0=-1 ,
                    sigma0=1,
                    sigma=matrix(rep(1,nYears*nSites),nrow=nSites,ncol=nYears)),
              list(  mu0=1 ,
                      sigma0=1.1,
                      sigma=matrix(rep(1.1,nYears*nSites),nrow=nSites,ncol=nYears)) 
              )



################################################################################
# Set parameters to monitor
################################################################################

parsToMonitor = c("theta","mu","sigma","mu0","sigma0","alpha","alpha.std","fruitingPlantNumberSim")

################################################################################
# Set JAGS parameters
################################################################################

# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 5000
n.update = 5000
n.iterations = 100000
n.thin = 10

################################################################################
# Fit model 
################################################################################

# set random seed
set.seed(2)

# tuning (n.adapt)
jm = jags.model(paste0(dirJagsScripts,"survivorshipModel-binLik-betaPriorMeanHierLogitNC-partialPooling-singleSite.R"), 
                data = simDataBayesJags, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iterations = n.update)

# chain (n.iter)
samplesBinLinkBetaPriorPartialMeanHierLogitNC = coda.samples(jm, variable.names = c(parsToMonitor), 
                                                         n.iter = n.iterations, thin = n.thin)

saveRDS(samplesBinLinkBetaPriorPartialMeanHierLogitNC,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/singleSiteBinLinkBetaPriorPartialMeanHierLogitNC.rds")

# 
# MCMCsummary(samplesBinLinkBetaPriorPartialMeanHierLogitNC,params=c("mu","sigma"))
# 
# 
# # compare centered vs. noncentered
# 
# MCMCsummary(samplesBinLinkBetaPriorPartialMeanLogit,params=c("mu","sigma"))
# MCMCsummary(samplesBinLinkBetaPriorPartialMeanLogitNC,params=c("mu","sigma","mu0","sigma0"))
# 
# MCMCtrace(samplesBinLinkBetaPriorPartialMeanLogitNC,params=c("mu"))
# MCMCtrace(samplesBinLinkBetaPriorPartialMeanLogitNC,params=c("sigma"))
# 
# 
# ## Comparison 1
# par(mfrow=c(1,2))
# 
# # centered
# alpha=MCMCchains(samplesBinLinkBetaPriorPartialMeanLogit,params="alpha")
# sigma=MCMCchains(samplesBinLinkBetaPriorPartialMeanLogit,params="sigma")
# 
# plot(alpha[,1],log(sigma[,1]))
# 
# # noncentered
# sigma=MCMCchains(samplesBinLinkBetaPriorPartialMeanLogitNC,params="sigma")
# alpha.std=MCMCchains(samplesBinLinkBetaPriorPartialMeanLogitNC,params="alpha.std")
# 
# plot(alpha.std[,1],log(sigma[,1]))
# 
# ## Comparison 2
# par(mfrow=c(1,2))
# 
# # centered
# mu=MCMCchains(samplesBinLinkBetaPriorPartialMeanLogit,params="mu")
# sigma=MCMCchains(samplesBinLinkBetaPriorPartialMeanLogit,params="sigma")
# 
# plot(mu[,1],log(sigma[,1]))
# 
# # noncentered
# sigma=MCMCchains(samplesBinLinkBetaPriorPartialMeanLogitNC,params="sigma")
# mu=MCMCchains(samplesBinLinkBetaPriorPartialMeanLogitNC,params="mu")
# 
# plot(mu[,1],log(sigma[,1]))
# 
# #saveRDS(samplesBinLinkBetaPriorPartialMean,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/singleSiteBinLinkBetaPriorPartialMean.rds")
# 
# sigma=MCMCchains(samplesBinLinkBetaPriorPartialMeanLogitNC,params="sigma")
# 
# par(mfrow=c(1,1))
# plot(alpha[,1],log(sigma[,1]))
# plot(alpha[,3],log(sigma[,1]))
# 
# alpha=MCMCchains(samplesBinLinkBetaPriorPartialMeanLogit,params="alpha")
# sigma=MCMCchains(samplesBinLinkBetaPriorPartialMeanLogit,params="sigma")
# 
# par(mfrow=c(1,1))
# plot(alpha[,1],log(sigma[,1]))
# plot(alpha[,3],log(sigma[,1]))
# 
# 
# MCMCsummary(samplesBinLinkBetaPriorPartialMean,params=c("phi","kappa"))
# MCMCsummary(samplesBinLinkBetaPriorPartialMeanHier,params=c("phi","kappa"))
# MCMCsummary(samplesBinLinkBetaPriorPartialMeanLogit,params=c("mu","sigma"))
# MCMCsummary(samplesBinLinkBetaPriorPartialMeanLogitNC,params=c("mu","sigma"))
