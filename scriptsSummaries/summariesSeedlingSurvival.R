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

sigma_py<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p_1[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(mu = mean(p_1),
                   ci.lo95 = quantile(p_1,probs=0.025), 
                   ci.hi95 = quantile(p_1,probs=0.975),
                   med = median(p_1), 
                   hpdi.lo95 = HPDI(p_1,.95)[1], 
                   hpdi.hi95 = HPDI(p_1,.95)[2]
  )

sigma_py<-sigma_py %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(site,year)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::select(site,year,mu,ci.lo95,ci.hi95,med,hpdi.lo95,hpdi.hi95)

# saveRDS(sigma_py,file=paste0(dirPars,"sigmaSummary.RDS"))
write.csv(sigma_py,file=paste0(dirPars,"sigmaSummary.csv"))
