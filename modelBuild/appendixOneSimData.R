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

site = "TEST"
year = as.character(rep(1:30,each=1))
n = rep(100, length(year))
# p = rep(seq(.275,.725,by=.05),each=100)
p = rep(0.5, length(year))
y = rbinom(length(year), size = n, prob = p)

simData =  data.frame( site, year, seedlingNumber = n, fruitingPlantNumber = y )

simDataBayesJags <- tidybayes::compose_data(simData)

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/"
simFiles <- paste0(directory,list.files(directory))

simData <- readRDS(simFiles[[9]])
# 
# ggplot(simData %>% dplyr::filter(site=="EC")) +
#   geom_histogram(aes(fruitingPlantNumber/seedlingNumber)) +
#   facet_wrap(~year) + theme_bw()
simData=simData %>%
  dplyr::filter(site=="EC"&year==2009) %>%
  dplyr::select(-year) %>%
  tidyr::unite("plot",c(transect:plot)) %>%
  dplyr::rename(year=plot)

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

inits = list(list(p = matrix(rep(.1,nSites*nYears),nrow=nSites,ncol=nYears)),
             list(p = matrix(rep(.5,nSites*nYears),nrow=nSites,ncol=nYears)),
             list(p = matrix(rep(.9,nSites*nYears),nrow=nSites,ncol=nYears))) 

################################################################################
# Set parameters to monitor
################################################################################

parsToMonitor = c("p")

################################################################################
# Set JAGS parameters
################################################################################

# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 500
n.update = 5000
n.iterations = 1000
n.thin = 1

################################################################################
# Fit model 
################################################################################

# set random seed
set.seed(2)

# tuning (n.adapt)
jm = jags.model(paste0(dirJagsScripts,"survivorshipModel-binLik-betaPrior-completePooling.R"), 
                data = simDataBayesJags, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iterations = n.update)

# chain (n.iter)
samplesBinLinkBetaPriorComplete = coda.samples(jm, variable.names = c(parsToMonitor), n.iter = n.iterations, thin = n.thin)

saveRDS(samplesBinLinkBetaPriorComplete,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/simBinLinkBetaPriorComplete.rds")

################################################################################
################################################################################
# JAGS fit for model with binomial likelihood, beta prior, partial pooling
# parameterized via mode
################################################################################
################################################################################

################################################################################
# Initial values
################################################################################

nSites = simDataBayesJags$n_site
nYears = simDataBayesJags$n_year

kappaInit = 100;

inits = list( list( #theta = matrix(rep(.1,nSites*nYears),nrow=nSites,ncol=nYears),
  omega=rep(.5,nSites) ,
  kappaMinusTwo=rep(98,nSites) ),
  list( #theta = matrix(rep(.5,nSites*nYears),nrow=nSites,ncol=nYears),
    omega=rep(.5,nSites) ,
    kappaMinusTwo=rep(98,nSites) ),
  list( #theta = matrix(rep(.9,nSites*nYears),nrow=nSites,ncol=nYears),
    omega=rep(.5,nSites) ,
    kappaMinusTwo=rep(98,nSites) )
)

################################################################################
# Set parameters to monitor
################################################################################

parsToMonitor = c("theta","omega","kappa")

################################################################################
# Set JAGS parameters
################################################################################

# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 500
n.update = 5000
n.iterations = 10000
n.thin = 1

################################################################################
# Fit model 
################################################################################

# set random seed
set.seed(2)

# tuning (n.adapt)
jm = jags.model(paste0(dirJagsScripts,"survivorshipModel-binLik-betaPriorMode-partialPooling.R"), 
                data = simDataBayesJags, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iterations = n.update)

# chain (n.iter)
samplesBinLinkBetaPriorPartialMode = coda.samples(jm, variable.names = c(parsToMonitor), 
                                                  n.iter = n.iterations, thin = n.thin)

saveRDS(samplesBinLinkBetaPriorPartialMode,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/simBinLinkBetaPriorPartialMode.rds")

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

inits = list( list( phi=rep(.1,nSites) ,
                    kappa=rep(1.1,nSites) ),
              list( phi=rep(.5,nSites) ,
                    kappa=rep(1.5,nSites) ),
              list( phi=rep(.9,nSites) ,
                    kappa=rep(2,nSites) )
)

################################################################################
# Set parameters to monitor
################################################################################

parsToMonitor = c("theta","phi","kappa")

################################################################################
# Set JAGS parameters
################################################################################

# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 500
n.update = 5000
n.iterations = 10000
n.thin = 1

################################################################################
# Fit model 
################################################################################

# set random seed
set.seed(2)

# tuning (n.adapt)
jm = jags.model(paste0(dirJagsScripts,"survivorshipModel-binLik-betaPriorMean-partialPooling.R"), 
                data = simDataBayesJags, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iterations = n.update)

# chain (n.iter)
samplesBinLinkBetaPriorPartialMean = coda.samples(jm, variable.names = c(parsToMonitor), 
                                                  n.iter = n.iterations, thin = n.thin)


saveRDS(samplesBinLinkBetaPriorPartialMean,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/simBinLinkBetaPriorPartialMean.rds")


#################################################################################
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/"
simFiles <- paste0(directory,list.files(directory))

samplesBinLinkBetaPriorComplete <- readRDS(simFiles[[1]])
samplesBinLinkBetaPriorPartialMean <- readRDS(simFiles[[2]])
samplesBinLinkBetaPriorPartialMode <- readRDS(simFiles[[3]])
simData <- readRDS(simFiles[[4]])

nSites = length(unique(simData$site))
nYears = length(unique(simData$year))

#MCMCsummary(samplesBinLinkBetaPriorPartialMean)
#MCMCsummary(samplesBinLinkBetaPriorPartialMode)

stat.summary <- function(x){
  
  x %<>% recover_types(simData)
  
  tmp <- x %>% 
    tidybayes::spread_draws(theta[site,year]) %>%
    dplyr::summarise(ci.lo = quantile(theta, prob = c(.025)),
                     med = quantile(theta, prob = c(.5)),
                     ci.hi = quantile(theta, prob = c(.975)))
  return(tmp)
}

pBinLinkBetaPriorComplete <- MCMCchains(samplesBinLinkBetaPriorComplete,params="p")

par(mfrow=c(3,1))
plot(density(pBinLinkBetaPriorComplete[,1]), type='n',
     xlim=c(-.2,1.2),ylim=c(0,15))
for(i in 1:nYears){
  lines(density(pBinLinkBetaPriorComplete[,i]))
}

thetaBinLinkBetaPriorPartialMode <- MCMCchains(samplesBinLinkBetaPriorPartialMode,params="theta")

plot(density(thetaBinLinkBetaPriorPartialMode[,1]), type='n', 
     xlim=c(-.2,1.2),ylim=c(0,15))
for(i in 1:nYears){
  lines(density(thetaBinLinkBetaPriorPartialMode[,i]))
}
 
omega<-MCMCchains(samplesBinLinkBetaPriorPartialMode,params="omega")
 kappa<-MCMCchains(samplesBinLinkBetaPriorPartialMode,params="kappa")

alpha = omega*(kappa-2)+1
beta = (1-omega)*(kappa-2)+1
#plot(seq(0,1,by=0.001),dbeta(seq(0,1,by=0.001),median(alpha[1,]),median(beta[1,])),type='n')
for(i in 1:20){
lines(seq(0,1,by=0.001),
      dbeta(seq(0,1,by=0.001),median(alpha[i,]),median(beta[i,])),
      type='l',lwd=.5,col='red')
}


thetaBinLinkBetaPriorPartialMean <- MCMCchains(samplesBinLinkBetaPriorPartialMean,params="theta")

plot(density(thetaBinLinkBetaPriorPartialMean[,1]), type='n', 
     xlim=c(-.2,1.2),ylim=c(0,15))
for(i in 1:nYears){
  lines(density(thetaBinLinkBetaPriorPartialMean[,i]))
}

phi<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params="phi")
kappa<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params="kappa")
alpha = kappa*phi
beta = kappa*(1-phi)
# plot(seq(0,1,by=0.001),dbeta(seq(0,1,by=0.001),median(alpha[1,]),median(beta[1,])),type='n')
for(i in 1:20){
  lines(seq(0,1,by=0.001),
        dbeta(seq(0,1,by=0.001),median(alpha[i,]),median(beta[i,])),
        type='l',lwd=.5,col='red')
}

MCMCsummary(samplesBinLinkBetaPriorComplete)
MCMCsummary(samplesBinLinkBetaPriorPartialMean)


thetaBinLinkBetaPriorPartialMode <- stat.summary(samplesBinLinkBetaPriorPartialMode)
thetaBinLinkBetaPriorPartialMean <- stat.summary(samplesBinLinkBetaPriorPartialMean)

plot(thetaBinLinkBetaPriorPartialMode$med,thetaBinLinkBetaPriorPartialMean$med, 
     pch=16,cex = 0.5,xlim=c(0,1),ylim=c(0,1))
# segments(x0=thetaBinLinkBetaPriorPartialMode$med,y0=thetaBinLinkBetaPriorPartialMean$ci.lo,
#          x1=thetaBinLinkBetaPriorPartialMode$med,y1=thetaBinLinkBetaPriorPartialMean$ci.hi)
# segments(x0=thetaBinLinkBetaPriorPartialMode$ci.lo,y0=thetaBinLinkBetaPriorPartialMean$med,
#          x1=thetaBinLinkBetaPriorPartialMode$ci.hi,y1=thetaBinLinkBetaPriorPartialMean$med)
abline(a=0,b=1)

mle<-simData %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(y = sum(fruitingPlantNumber), n = sum(seedlingNumber)) %>%
  dplyr::mutate(pHat = y/n) %>%
  dplyr::left_join(thetaBinLinkBetaPriorPartialMode)

points(x=mle$med,y=mle$pHat,pch=2)

## summary
MCMCsummary(samplesBinLinkBetaPriorPartialMean,params=c("phi","kappa"))
MCMCsummary(samplesBinLinkBetaPriorPartialMode,params=c("omega","kappa"))

# calculate difference of posterior of mode - posterior of ean
post.chain.pars.mode <-MCMCvis::MCMCchains(samplesBinLinkBetaPriorPartialMode,params="theta")
post.chain.pars.mean <-MCMCvis::MCMCchains(samplesBinLinkBetaPriorPartialMean,params="theta")

# for plotting

samplesizeData <- simData %>%
  dplyr::rename(n=seedlingNumber)

samplesizeData$site <- as.factor(samplesizeData$site)
samplesizeData$year <- as.numeric(as.character(samplesizeData$year))

estimateComparisonDataArranged <- samplesizeData %>%
  dplyr::arrange(year)

tmp<-matrix(NA,ncol=(nSites*nYears),nrow=length(samplesBinLinkBetaPriorPartialMode)*dim(samplesBinLinkBetaPriorPartialMode[[1]]))
for(i in 1:(nSites*nYears)){
  tmp[,i]<-post.chain.pars.mode[,i]-post.chain.pars.mean[i]
}
delta.conf <- apply(tmp,2,quantile,probs=c(.025,.5,.975))

plot(estimateComparisonDataArranged$n,t(delta.conf)[,2],
     pch=16,cex=0.5,
     ylim=c(-1,1),
     xlab="Sample size (n)",
     ylab="Comparison of the beta-binomial parameterized via mode\n and the beta-binomial parameterized via mean\n (median and 95% CI)")
segments(x0=estimateComparisonDataArranged$n,y0=t(delta.conf)[,1],
         x1=estimateComparisonDataArranged$n,y1=t(delta.conf)[,3])
abline(h=0)

estimateDistributionComparison<-data.frame(cbind(estimateComparisonDataArranged,t(delta.conf)))
estimateDistributionComparison <- estimateDistributionComparison %>%
  dplyr::select(-pHat)
names(estimateDistributionComparison) = c("site","year","n","ci.lo","med","ci.hi")

library(ggplot2)

# replicate graph above
ggplot(data=estimateDistributionComparison) +
  geom_point(aes(x=n,y=med)) +
  geom_errorbar(aes(x=n, ymin=ci.lo, ymax=ci.hi), width=.1) +
  theme_classic()

# facet by site
ggplot(data=estimateDistributionComparison) +
  geom_point(aes(x=n,y=med)) +
  geom_errorbar(aes(x=n, ymin=ci.lo, ymax=ci.hi), width=.1) +
  facet_wrap(~site) +
  geom_hline(yintercept=0) +
  theme_classic()

# facet by year
ggplot(data=estimateDistributionComparison) +
  geom_point(aes(x=n,y=med)) +
  geom_errorbar(aes(x=n, ymin=ci.lo, ymax=ci.hi), width=.1) +
  facet_wrap(~year) +
  geom_hline(yintercept=0) +
  theme_classic()