#################################################################################
################################################################################
################################################################################
# Code for figures to compare the following modeling approaches for the seedling survivorship data
# 
# Scripts by Gregor Siegmund
# fausto.siegmund@gmail.com
# last updated 04-22-2020
#################################################################################
#################################################################################
#################################################################################
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)

library(MCMCvis)
library(tidybayes)
library(tidyverse)
library(magrittr)
library(bayesplot)
library(rethinking)


directory = "/Users/Gregor/Dropbox/dataLibrary/posteriors/"
simFiles <- paste0(directory,list.files(directory))

mcmcSamples <- readRDS(simFiles[[11]])

dirPars = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/parameterSummary/"
dir.create(file.path(dirPars), showWarnings = FALSE)

################################################################################
# Data directory
#################################################################################

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedlingSurvival/"
dataFiles <- paste0(directory,list.files(directory))

data <- readRDS(dataFiles[[2]])
censusSeedlingsFruitingPlants <- readRDS(dataFiles[[1]])

################################################################################
# Seedling survival
#################################################################################

siteIndex <- data.frame(siteIndex=unique(censusSeedlingsFruitingPlants$site),site=1:20)
yearIndex <- data.frame(yearIndex=as.numeric(unique(censusSeedlingsFruitingPlants$year)),
                        year=unique(unique(data$year)))

sigma0_p<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(sigma0_1[site]) %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(mu = mean(sigma0_1),
                   ci.lo95 = quantile(sigma0_1,probs=0.025), 
                   ci.hi95 = quantile(sigma0_1,probs=0.975),
                   med = median(sigma0_1), 
                   hpdi.lo95 = HPDI(sigma0_1,.95)[1], 
                   hpdi.hi95 = HPDI(sigma0_1,.95)[2]
  )

sigma0_p<-sigma0_p %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(site)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::select(site,mu,ci.lo95,ci.hi95,med,hpdi.lo95,hpdi.hi95)

ggplot(sigma0_p)+
  geom_linerange(aes(x=site,ymin=ci.lo95,ymax=ci.hi95),size=.25) +
  geom_point(aes(x=site,y=med))

# saveRDS(sigma_py,file=paste0(dirPars,"sigmaSummary.RDS"))


sigma0_py<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(sigma_1[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(mu = mean(sigma_1),
                   ci.lo95 = quantile(sigma_1,probs=0.025), 
                   ci.hi95 = quantile(sigma_1,probs=0.975),
                   med = median(sigma_1), 
                   hpdi.lo95 = HPDI(sigma_1,.95)[1], 
                   hpdi.hi95 = HPDI(sigma_1,.95)[2]
  )

sigma0_py<-sigma0_py %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(site,year)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::select(site,year,mu,ci.lo95,ci.hi95,med,hpdi.lo95,hpdi.hi95)

ggplot(sigma0_py)+
  geom_linerange(aes(x=year,ymin=ci.lo95,ymax=ci.hi95),size=.25) +
  geom_point(aes(x=year,y=med)) +
  facet_wrap(~site)

sampleSize = censusSeedlingsFruitingPlants %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(year = as.double(year)) 
  

trialNumber = censusSeedlingsFruitingPlants %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(trials = sum(seedlingNumber),
                   trialsVar = var(seedlingNumber)) %>%
  dplyr::mutate(year = as.double(year)) 

sigmaSamples <- sigma0_py %>%
  dplyr::left_join(sampleSize,by=c("site","year")) %>%
  dplyr::left_join(trialNumber,by=c("site","year"))


ggplot(data=sigmaSamples, aes(x=n,y=med)) +
  geom_point() +
  theme_bw() +
  xlab("Number of plots") + ylab("Median estimate of sigma")

ggplot(data=sigmaSamples, aes(x=trials,y=med)) +
  geom_point() +
  theme_bw() +
  xlab("Number of trials") + ylab("Median estimate of sigma") 

ggplot(data=sigmaSamples, aes(x=trialsVar,y=med)) +
  geom_point() +
  theme_bw() +
  xlab("Variance in trials") + ylab("Median estimate of sigma") +
  xlim(c(0,10))

sigmaSamples %>%
  dplyr::group_by(site) %>%
  dplyr::mutate(mean())

ggplot(sigmaSamples)+
  geom_linerange(aes(x=year,ymin=ci.lo95,ymax=ci.hi95),size=.25) +
  geom_point(aes(x=year,y=med,size=trials)) +
  facet_wrap(~site)
