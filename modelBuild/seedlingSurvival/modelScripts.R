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
load("/Users/Gregor/Dropbox/clarkiaSeedBanks/library/dataFromWorkflowFile/seedBagsData.rda")
df <- seedBags

# -------------------------------------------------------------------
# Clean and organize seed bag data
# -------------------------------------------------------------------
df$seedStart<-100

df$seedStart <- as.integer(df$seedStart)

df$totalJan <- ifelse(df$totalJan>100,100,df$totalJan)

# `MC III 5 13 1 is a problem (91 intact, 17 seedlings)`

# df <- df %>% dplyr::rename(bagNo=bag)

## FOR NOW REMOVE MISSING DATA AND PROBLEMS
# the main issue will be bags that were not recovered in January
# but then recovered in October; there is no way of getting an estimate on how many 
# seeds 'started' those trials

df<-subset(df,!is.na(df$totalJan))
df<-subset(df,!is.na(df$intactOct))
df<-subset(df,!(intactJan<intactOct))

seedBagExperiment <- df

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------
# -------------------------------------------------------------------

## filter the dataset for testing purposes
filterData<-function(x) {
  x %>%
    dplyr::filter(age==1) %>%
    dplyr::filter(site %in% unique(df$site)[1:2])
}


seedBagExperiment<-filterData(seedBagExperiment)

seedBagExperiment = seedBagExperiment %>%
  dplyr::mutate(year = as.factor(yearStart)) %>%
  dplyr::select(site,year,totalJan,seedStart)

data <- tidybayes::compose_data(seedBagExperiment)

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
n.iterations = 40000
n.thin = 10

set.seed(10)
dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/jagsScriptsBinomialModelMultilevel/")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Complete pooling of germination and viability trials
# Partial pooling of seed burial experiment (site level)
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# set inits for JAGS
inits = list(list(mu0 = rep(0,data$n_site), sigma0 = rep(.5,data$n_site),
                  sigma = matrix(rep(.5,data$n_site*data$n_year),nrow=data$n_site,ncol=data$n_year)),
             list(mu0 = rep(-1,data$n_site), sigma0 = rep(1,data$n_site),
                  sigma = matrix(rep(1,data$n_site*data$n_year),nrow=data$n_site,ncol=data$n_year)),
             list(mu0 = rep(1,data$n_site), sigma0 = rep(1.25,data$n_site),
                  sigma = matrix(rep(1.25,data$n_site*data$n_year),nrow=data$n_site,ncol=data$n_year)))

# # Call to JAGS
# 
# # tuning (n.adapt)
jm = jags.model(paste0(dir,"hierarchicalLogitCentered.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

parsToMonitor = c("theta","alpha","mu0","sigma0","mu","sigma","p")

# chain (n.iter)
samples.rjags = coda.samples(jm, variable.names = c(parsToMonitor), n.iter = n.iterations, thin = n.thin)


MCMCsummary(samples.rjags, params = c("mu0","sigma0"))
MCMCsummary(samples.rjags, params = c("mu","sigma","theta"))
# 
# par(mfrow=c(3,1))
# hist(MCMCchains(samples.rjags,params="p")[,1],breaks=100)
# hist(MCMCchains(samples.rjags,params="p")[,2],breaks=100)
# hist(MCMCchains(samples.rjags,params="p")[,3],breaks=100)

# fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBurial/")
# dir.create(file.path(fileDirectory), showWarnings = FALSE)
# 
# saveRDS(samples.rjags,file=paste0(fileDirectory,"samples.rjags.rds"))

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Complete pooling of germination and viability trials
# Partial pooling of seed burial experiment (site level)
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# set inits for JAGS
inits = list(list(mu0 = rep(0,data$n_site), sigma0 = rep(.5,data$n_site),
                  sigma = matrix(rep(.5,data$n_site*data$n_year),nrow=data$n_site,ncol=data$n_year)),
             list(mu0 = rep(-1,data$n_site), sigma0 = rep(1,data$n_site),
                  sigma = matrix(rep(1,data$n_site*data$n_year),nrow=data$n_site,ncol=data$n_year)),
             list(mu0 = rep(1,data$n_site), sigma0 = rep(1.25,data$n_site),
                  sigma = matrix(rep(1.25,data$n_site*data$n_year),nrow=data$n_site,ncol=data$n_year)))

# Call to JAGS

# tuning (n.adapt)
jm = jags.model(paste0(dir,"hierarchicalLogitNoncentered.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

parsToMonitor = c("theta","alpha","mu0","sigma0","mu","sigma","p_pop")

# chain (n.iter)
samples.rjags2 = coda.samples(jm, variable.names = c(parsToMonitor), n.iter = n.iterations, thin = n.thin)

MCMCsummary(samples.rjags2, params = c("mu0","sigma0"))
MCMCsummary(samples.rjags, params = c("mu0","sigma0"))

MCMCsummary(samples.rjags2, params = c("mu","sigma"))

# par(mfrow=c(1,1))
# hist(MCMCchains(samples.rjags2,params="p_pop")[,1],breaks=100)
# hist(MCMCchains(samples.rjags2,params="p_pop")[,2],breaks=100)
# hist(MCMCchains(samples.rjags2,params="p_pop")[,3],breaks=100)
# hist(MCMCchains(samples.rjags2,params="p_pop")[,4],breaks=100)
# hist(MCMCchains(samples.rjags2,params="p_pop")[,5],breaks=100)

par(mfrow=c(1,1))
plot(density((boot::inv.logit(MCMCchains(samples.rjags2,params="mu0")[,1]))),ylim=c(0,20))
lines(density(MCMCchains(samples.rjags2,params="p_pop")[,1]),col='red')
lines(density(boot::inv.logit(MCMCchains(samples.rjags2,params="mu")[,c(1)])),lty='dashed')
lines(density(boot::inv.logit(MCMCchains(samples.rjags2,params="mu")[,c(3)])),lty='dashed')
lines(density(boot::inv.logit(MCMCchains(samples.rjags2,params="mu")[,c(5)])),lty='dashed')

par(mfrow=c(1,1))
plot(density((boot::inv.logit(MCMCchains(samples.rjags2,params="mu0")[,2]))),ylim=c(0,20))
lines(density(MCMCchains(samples.rjags2,params="p_pop")[,2]),col='red')
lines(density(boot::inv.logit(MCMCchains(samples.rjags2,params="mu")[,c(2)])),lty='dashed')
lines(density(boot::inv.logit(MCMCchains(samples.rjags2,params="mu")[,c(4)])),lty='dashed')
lines(density(boot::inv.logit(MCMCchains(samples.rjags2,params="mu")[,c(6)])),lty='dashed')

plot(density((boot::inv.logit(MCMCchains(samples.rjags2,params="p_pop")[,1]))),
     ylim=c(0,10),type='n',xlim=c(0,1))
for(i in 1:2){
  lines(density(MCMCchains(samples.rjags2,params="p_pop")[,i]),col=i)

}


# fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBurial/")
# dir.create(file.path(fileDirectory), showWarnings = FALSE)
# 
# saveRDS(samples.rjags,file=paste0(fileDirectory,"samples.rjags.rds"))

data$freq <- data$totalJan/data$seedStart

summary.LogitNoncentered<-apply(MCMCchains(samples.rjags2,params="theta"),2,quantile, probs = c(.05,.5,.95))

summaryList<-list(summary.LogitNoncentered)

#pdf(paste0(dirFigures,"shrinkage.pdf"), width=8, height=6)
est=seedBagExperiment %>%
  dplyr::group_by(site) %>% 
  dplyr::summarise(totalJan=sum(totalJan),seedStart=sum(seedStart)) %>%
  dplyr::mutate(prop = totalJan/seedStart)
par(mfrow=c(1,1))
for(i in 1:1){
  plot(x = data$freq,
       y = t(summaryList[[i]])[,2],
       xlab= "Proportion of seeds remaining in bag",
       ylab="MCMC estimate and 95% CI for probability",
       pch=16,cex=.5,
       xlim=c(0,1),ylim=c(0,1),
       col=data$site)
  segments(x0=data$freq,
           y0= t(summaryList[[i]])[,1],
           y1= t(summaryList[[i]])[,3],
           col=data$site)
  abline(a=0,b=1)
  abline(h=est$prop,col=as.factor(est$site))
}
#dev.off()

apply(boot::inv.logit(MCMCchains(samples.rjags2,params="mu0")),2,quantile, probs = c(.05,.5,.95))
seedBagExperiment %>%
  dplyr::group_by(site) %>% 
  dplyr::summarise(totalJan=sum(totalJan),seedStart=sum(seedStart)) %>%
  dplyr::mutate(prop = totalJan/seedStart)

apply(boot::inv.logit(MCMCchains(samples.rjags2,params="mu")),2,quantile, probs = c(.05,.5,.95))
seedBagExperiment %>%
  dplyr::group_by(site,year) %>% 
  dplyr::summarise(totalJan=sum(totalJan),seedStart=sum(seedStart)) %>%
  dplyr::mutate(prop = totalJan/seedStart)

plot((MCMCchains(samples.rjags2,params="mu0"))[,2],pch=16,cex=.1)
hist(boot::inv.logit(MCMCchains(samples.rjags2,params="mu0"))[,1],breaks=1000)
hist(boot::inv.logit(MCMCchains(samples.rjags2,params="mu0"))[,2],breaks=1000)
