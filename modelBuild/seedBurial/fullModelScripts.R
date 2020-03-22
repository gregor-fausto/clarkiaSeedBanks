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
# Import and organize viability trial data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# setwd and read data files
load("/Users/Gregor/Dropbox/clarkiaSeedBanks/library/dataFromWorkflowFile/viabilityRawData.rds")
df <- viabilityRawData

df <- df %>% 
  dplyr::select(-c(germPerc,germNot,viabPerc,viabPerc2,condTest))

#df <- df %>% dplyr::rename(bagNo=bag)
df$bag <-as.integer(as.numeric(df$bagNo))

viabilityExperiment <- df

# another check
seedBurial <- seedBagExperiment %>% dplyr::select(site,bagNo,round,age) %>%
  dplyr::mutate(seedBurial = 1)
viabilityTrial <- viabilityExperiment %>% dplyr::select(site,bag,round,age) %>%
  dplyr::mutate(viabilityTrial = 1)

joinedDF<-seedBurial %>% 
  dplyr::full_join(viabilityTrial) 

# dplyr::group_by(site,round,age) %>%
#   dplyr::count(seedBurial,viabilityTrial)

# one row is coded differently so that 
# viabStart=NA and viabStain=NA
# all others have viabStart=NA and viabStain=NA
# recode
viabilityExperiment[is.na(viabilityExperiment$viabStart),]$viabStart = 0

## data check
viabilityExperiment %>% dplyr::filter(germStart - germCount - viabStart<0) 

# filter out rows with problems
# these need to be corrected
viabilityExperiment<-viabilityExperiment %>% 
  dplyr::filter(germStart - germCount - viabStart >= 0)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------
# -------------------------------------------------------------------

## filter the dataset for testing purposes
filterData<-function(x) {
  x %>%
    dplyr::filter(age==1) #%>%
   # dplyr::filter(site %in% unique(df$site)[1:2])
}

seedBagExperiment<-filterData(seedBagExperiment)
viabilityExperiment<-filterData(viabilityExperiment)

# assign variable that combines site and bag; unique id for each bag
# for each dataset, create that unique identifier again and then
# use that to link it to the reference identifier created above
seedBagExperiment<-seedBagExperiment %>%
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="", remove=FALSE) %>%
  tidyr::unite(col='siteBag', c(site,bagNo), sep="", remove=FALSE) %>%
  dplyr::mutate(siteBag = as.factor(siteBag)) #%>%
# dplyr::left_join(referenceTable,by="id")

viabilityExperiment<-viabilityExperiment %>%
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="", remove=FALSE) %>%
  tidyr::unite(col='siteBag', c(site,bagNo), sep="", remove=FALSE) %>%
  dplyr::mutate(siteBag = as.factor(siteBag)) #%>%
# dplyr::left_join(referenceTable,by="id")

# once each identifier has been created and linked to the reference table
# and the dataset filtered, the dataset needs to be re-indexed
# this may be redundant?

# this line creates a unique id for the subsetted data that is then 
# used to index each of the 2 datasets
# and provides the reference set of bags that were included in the experiment

#https://community.rstudio.com/t/how-to-add-a-counter-to-each-group-in-dplyr/12986/2

referenceTable<-data.frame(id=union(seedBagExperiment$id, viabilityExperiment$id)) %>%
  dplyr::mutate(idNo = 1:length(id)) 

seedBagExperiment<-seedBagExperiment %>%
  dplyr::left_join(referenceTable,by="id")

viabilityExperiment<-viabilityExperiment %>%
  dplyr::left_join(referenceTable,by="id")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------
# -------------------------------------------------------------------



seedBagExperiment = seedBagExperiment %>%
  dplyr::mutate(year = as.factor(yearStart)) %>%
  dplyr::select(site,year,totalJan,seedStart,seedlingJan) %>%
  dplyr::rename(siteBags = site,
                yearBags = year)

viabilityExperiment = viabilityExperiment %>%
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

data <- tidybayes::compose_data(seedBagExperiment,viabilityExperiment)

data$n = dim(seedBagExperiment)[1]
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
                  mu0_g = rep(0,data$n_siteBags), sigma0_g = rep(.5,data$n_siteBags),
                  sigma_g = matrix(rep(.5,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
                  mu0_v = rep(0,data$n_siteBags), sigma0_v = rep(.5,data$n_siteBags),
                  sigma_v = matrix(rep(.5,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags)),
             list(mu0_1 = rep(-1,data$n_siteBags), sigma0_1 = rep(1,data$n_siteBags),
                  sigma_1 = matrix(rep(1,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
                  mu0_2 = rep(-1,data$n_siteBags), sigma0_2 = rep(1,data$n_siteBags),
                  sigma_2 = matrix(rep(1,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
                  mu0_g = rep(-1,data$n_siteBags), sigma0_g = rep(1,data$n_siteBags),
                  sigma_g = matrix(rep(1,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
                  mu0_v = rep(-1,data$n_siteBags), sigma0_v = rep(1,data$n_siteBags),
                  sigma_v = matrix(rep(1,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags)),
             list(mu0_1 = rep(1,data$n_siteBags), sigma0_1 = rep(1.25,data$n_siteBags),
                  sigma_1 = matrix(rep(1.25,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
                  mu0_2 = rep(1,data$n_siteBags), sigma0_2 = rep(1.25,data$n_siteBags),
                  sigma_2 = matrix(rep(1.25,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
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
parsToMonitor_g = c("theta_g","mu0_g","sigma0_g","mu_g","sigma_g","p_g")
parsToMonitor_v = c("theta_v","mu0_v","sigma0_v","mu_v","sigma_v","p_v")
parsToMonitor_deriv = c("nu_1","s1","s2")


# chain (n.iter)
samples.rjags = coda.samples(jm, 
                             variable.names = c(parsToMonitor_1,parsToMonitor_2,
                                                parsToMonitor_g,parsToMonitor_v,
                                                parsToMonitor_deriv), 
                             n.iter = n.iterations, thin = n.thin)


MCMCsummary(samples.rjags, params = c("mu0_1","sigma0_1"))
MCMCsummary(samples.rjags, params = c("mu0_2","sigma0_2"))
MCMCsummary(samples.rjags, params = c("mu0_g","sigma0_g"))
MCMCsummary(samples.rjags, params = c("mu0_v","sigma0_v"))

MCMCsummary(samples.rjags, params = c("mu_1","sigma_1"))
MCMCsummary(samples.rjags, params = c("mu_2","sigma_2"))

MCMCsummary(samples.rjags, params = c("p_1","p_2","p_g","p_v"))

MCMCsummary(samples.rjags, params = c("s1","s2"))

fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBurial/")
dir.create(file.path(fileDirectory), showWarnings = FALSE)
# 
saveRDS(samples.rjags,file=paste0(fileDirectory,"seedBurialSamples.rds"))

# # -------------------------------------------------------------------
# # -------------------------------------------------------------------
# # Visualize
# # -------------------------------------------------------------------
# # -------------------------------------------------------------------

# par(mfrow=c(2,1))
# plot(density(MCMCchains(samples.rjags,params="p_1")[,1]),ylim=c(0,10))
# lines(density(MCMCchains(samples.rjags,params="p_2")[,1]),lty='dashed')
# lines(density(MCMCchains(samples.rjags,params="p_pop")[,1]),col='red')
# 
# plot(density(MCMCchains(samples.rjags,params="p_1")[,2]),ylim=c(0,10))
# lines(density(MCMCchains(samples.rjags,params="p_2")[,2]),lty='dashed')
# lines(density(MCMCchains(samples.rjags,params="p_pop")[,2]),col='red')
# 
# 
# lines(density(boot::inv.logit(MCMCchains(samples.rjags2,params="mu")[,c(1)])),lty='dashed')
# lines(density(boot::inv.logit(MCMCchains(samples.rjags2,params="mu")[,c(3)])),lty='dashed')
# lines(density(boot::inv.logit(MCMCchains(samples.rjags2,params="mu")[,c(5)])),lty='dashed')

plot(density((boot::inv.logit(MCMCchains(samples.rjags,params="s1")[,1]))),
     ylim=c(0,10),type='n',xlim=c(0,1))
for(i in 1:2){
  lines(density(MCMCchains(samples.rjags,params="s1")[,i]),col=i)
  lines(density(MCMCchains(samples.rjags,params="s2")[,i]),col=i,lty="dashed")
}


# fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBurial/")
# dir.create(file.path(fileDirectory), showWarnings = FALSE)
# 
# saveRDS(samples.rjags,file=paste0(fileDirectory,"samples.rjags.rds"))

data$freq <- data$totalJan/data$seedStart

summary.LogitCentered<-apply(MCMCchains(samples.rjags,params="theta_1"),2,quantile, probs = c(.05,.5,.95))

summaryList<-list(summary.LogitCentered)

#pdf(paste0(dirFigures,"shrinkage.pdf"), width=8, height=6)
est=seedBagExperiment %>%
  dplyr::group_by(siteBags) %>% 
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
       col=data$siteBags)
  segments(x0=data$freq,
           y0= t(summaryList[[i]])[,1],
           y1= t(summaryList[[i]])[,3],
           col=data$siteBags)
  abline(a=0,b=1)
  abline(h=est$prop,col=as.factor(est$siteBags))
}
#dev.off()

plot(density((boot::inv.logit(MCMCchains(samples.rjags,params="s1")[,1]))),
     ylim=c(0,10),type='n',xlim=c(0,1))
for(i in 1:2){
  lines(density(MCMCchains(samples.rjags,params="s1")[,i]),col=i)
  lines(density(MCMCchains(samples.rjags,params="s2")[,i]),col=i,lty="dashed")
}


# fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBurial/")
# dir.create(file.path(fileDirectory), showWarnings = FALSE)
# 
# saveRDS(samples.rjags,file=paste0(fileDirectory,"samples.rjags.rds"))

data$freq <- data$seedlingJan/data$totalJan

summary.LogitCentered<-apply(MCMCchains(samples.rjags,params="theta_2"),2,quantile, probs = c(.05,.5,.95))

summaryList<-list(summary.LogitCentered)

#pdf(paste0(dirFigures,"shrinkage.pdf"), width=8, height=6)
est=seedBagExperiment %>%
  dplyr::group_by(siteBags) %>% 
  dplyr::summarise(totalJan=sum(totalJan),seedlingJan=sum(seedlingJan)) %>%
  dplyr::mutate(prop = seedlingJan/totalJan)
par(mfrow=c(1,1))
for(i in 1:1){
  plot(x = data$freq,
       y = t(summaryList[[i]])[,2],
       xlab= "Proportion of germinants remaining in bag",
       ylab="MCMC estimate and 95% CI for probability",
       pch=16,cex=.5,
       xlim=c(0,1),ylim=c(0,1),
       col=data$siteBags)
  segments(x0=data$freq,
           y0= t(summaryList[[i]])[,1],
           y1= t(summaryList[[i]])[,3],
           col=data$siteBags)
  abline(a=0,b=1)
  abline(h=est$prop,col=as.factor(est$siteBags))
}
#dev.off()

# germination trials
data$freq <- data$germCount/data$germStart

summary.LogitCentered<-apply(MCMCchains(samples.rjags,params="theta_g"),2,quantile, probs = c(.05,.5,.95))

summaryList<-list(summary.LogitCentered)

#pdf(paste0(dirFigures,"shrinkage.pdf"), width=8, height=6)
est=viabilityExperiment %>%
  dplyr::group_by(siteViab) %>% 
  dplyr::summarise(germCount=sum(germCount),germStart=sum(germStart)) %>%
  dplyr::mutate(prop = germCount/germStart)
par(mfrow=c(1,1))
for(i in 1:1){
  plot(x = data$freq,
       y = t(summaryList[[i]])[,2],
       xlab= "Proportion of germinants in test",
       ylab="MCMC estimate and 95% CI for probability",
       pch=16,cex=.5,
       xlim=c(0,1),ylim=c(0,1),
       col=data$siteViab)
  segments(x0=data$freq,
           y0= t(summaryList[[i]])[,1],
           y1= t(summaryList[[i]])[,3],
           col=data$siteViab)
  abline(a=0,b=1)
  abline(h=est$prop,col=as.factor(est$siteViab))
}

# viability trials
data$freq <- data$viabStain/data$viabStart

summary.LogitCentered<-apply(MCMCchains(samples.rjags,params="theta_v"),2,quantile, probs = c(.05,.5,.95))

summaryList<-list(summary.LogitCentered)

#pdf(paste0(dirFigures,"shrinkage.pdf"), width=8, height=6)
est=viabilityExperiment %>%
  dplyr::group_by(siteViab) %>% 
  dplyr::summarise(viabStain=sum(viabStain,na.rm=TRUE),viabStart=sum(viabStart)) %>%
  dplyr::mutate(prop = viabStain/viabStart)
par(mfrow=c(1,1))
for(i in 1:1){
  plot(x = data$freq,
       y = t(summaryList[[i]])[,2],
       xlab= "Proportion of germinants in test",
       ylab="MCMC estimate and 95% CI for probability",
       pch=16,cex=.5,
       xlim=c(0,1),ylim=c(0,1),
       col=data$siteViab)
  segments(x0=data$freq,
           y0= t(summaryList[[i]])[,1],
           y1= t(summaryList[[i]])[,3],
           col=data$siteViab)
  abline(a=0,b=1)
  abline(h=est$prop,col=as.factor(est$siteViab))
}