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
library(HDInterval)
library(coda)
library(rethinking)

dirFigures = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/"


directory = "/Users/Gregor/Dropbox/dataLibrary/posteriors/"
posteriorFiles <- paste0(directory,list.files(directory))

mcmcSamples <- readRDS(posteriorFiles[[7]])

dirPars = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/parameterSummary/"
dir.create(file.path(dirPars), showWarnings = FALSE)

################################################################################
# Read in raw data
#################################################################################

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/"
dataFiles <- paste0(directory,list.files(directory))

data <- readRDS(dataFiles[[5]])

################################################################################
# Create composite
#################################################################################

countFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countFruitsPerPlantAllPlots.RDS")

siteNames = unique(countFruitsPerPlantAllPlots$site)

# recover site indices and year indices
siteIndex <- data.frame(siteIndex=unique(countFruitsPerPlantAllPlots$site),site=unique(data$site))
yearIndex <- data.frame(yearIndex=unique(countFruitsPerPlantAllPlots$year),
                        year=unique(data$year)) 

tfeDF<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(mu_py_tfe[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(mu = mean(mu_py_tfe),
                   ci.lo95 = quantile(mu_py_tfe,probs=0.025), 
                   ci.hi95 = quantile(mu_py_tfe,probs=0.975),
                   med = median(mu_py_tfe), 
                   hpdi.lo95 = HPDI(mu_py_tfe,.95)[1], 
                   hpdi.hi95 = HPDI(mu_py_tfe,.95)[2]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  dplyr::select(-c(site,year)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex)


countSeedPerDamagedFruit <- readRDS(dataFiles[[2]])

siteNames = unique(countSeedPerDamagedFruit$site4) 

# recover site indices and year indices
siteIndex <- data.frame(siteIndex=unique(countSeedPerDamagedFruit$site4),site=unique(data$site4))
yearIndex <- data.frame(yearIndex=unique(countSeedPerDamagedFruit$year4),
                        year=unique(data$year4)) 

tfeCompDF<- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(mu_py_tfe_comp[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(mu = mean(mu_py_tfe_comp),
                   ci.lo95 = quantile(mu_py_tfe_comp,probs=0.025), 
                   ci.hi95 = quantile(mu_py_tfe_comp,probs=0.975),
                   med = median(mu_py_tfe_comp), 
                   hpdi.lo95 = HPDI(mu_py_tfe_comp,.95)[1], 
                   hpdi.hi95 = HPDI(mu_py_tfe_comp,.95)[2]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  dplyr::select(-c(site,year)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) 

# manually filter out S22 2013
tfeCompDF<- tfeCompDF %>% 
  dplyr::filter(site!="S22"|year!=2013) %>%
  dplyr::filter(site!="OSR"|year!=2013) %>%
  dplyr::filter(site!="LCE"|year!=2013) %>%
  dplyr::filter(site!="GCN"|year!=2013)

reference<-tfeCompDF %>% dplyr::select(site,year) %>% unique()
reference$year <- as.integer(reference$year)
tfeCompDF$year<-as.integer(tfeCompDF$year)

# undamaged
countFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/countUndamagedDamagedFruitsPerPlantAllPlots.rds")
siteIndex <- data.frame(siteIndex=unique(countFruitsPerPlantAllPlots$site),site=unique(data$site))
yearIndex <- data.frame(yearIndex=unique(countFruitsPerPlantAllPlots$year),
                        year=unique(data$year)) 

ufDF<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(mu_py_und[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(mu = mean(mu_py_und),
                   ci.lo95 = quantile(mu_py_und,probs=0.025), 
                   ci.hi95 = quantile(mu_py_und,probs=0.975),
                   med = median(mu_py_und), 
                   hpdi.lo95 = HPDI(mu_py_und,.95)[1], 
                   hpdi.hi95 = HPDI(mu_py_und,.95)[2]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  dplyr::select(-c(site,year)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex)

ufDF$year<-as.integer(ufDF$year)

referenceUF<-ufDF %>% dplyr::select(site,year) %>% unique()

ufOnly<-setdiff(referenceUF,reference) %>%
  dplyr::mutate(site.year=paste0(site,year))

ufDF_missing<-ufDF %>%
  dplyr::mutate(site.year=paste0(site,year)) %>%
  dplyr::filter(site.year %in%ufOnly$site.year) %>%
  dplyr::select(-site.year)

# bind all datasets
tfeFullDF<-bind_rows(tfeDF,tfeCompDF,ufDF_missing)

# fill in with NAs
referenceDF <- expand.grid(site=unique(tfeFullDF$site),year=unique(tfeFullDF$year))

tfeFullDF<-tfeFullDF %>%
  dplyr::full_join(referenceDF,by=c("site","year"))

################################################################################
# Summaries
#################################################################################

tfeFullDF<-tfeFullDF %>%
  dplyr::select(site,year,mu,ci.lo95,ci.hi95,med,hpdi.lo95,hpdi.hi95)

#saveRDS(tfeFullDF,file=paste0(dirPars,"fruitEquivalentsPerPlantSummary.RDS"))
write.csv(tfeFullDF,file=paste0(dirPars,"fruitEquivalentsPerPlantSummary.csv"))
