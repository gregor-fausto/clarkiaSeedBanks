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
posteriorFiles <- paste0(directory,list.files(directory))

mcmcSamples <- readRDS(posteriorFiles[[7]])

dirPars = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/parameterSummary/"
dir.create(file.path(dirPars), showWarnings = FALSE)

################################################################################
# Data directory
#################################################################################

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/"
dataFiles <- paste0(directory,list.files(directory))

data <- readRDS(dataFiles[[5]])

################################################################################
# Seedling survival
#################################################################################

countSeedPerFruit <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countSeedPerFruit.RDS")
countSeedPerFruit <- countSeedPerFruit %>% dplyr::filter(demography==1)
siteNames = unique(countSeedPerFruit$site)

# recover site indices and year indices
siteIndex <- data.frame(siteIndex=unique(countSeedPerFruit$site),site=unique(data$site3))
yearIndex <- data.frame(yearIndex=unique(countSeedPerFruit$year),
                        year=unique(data$year3)) 


seeds_py<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(mu_py_seeds[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(mu = mean(mu_py_seeds),
                   ci.lo95 = quantile(mu_py_seeds,probs=0.025), 
                   ci.hi95 = quantile(mu_py_seeds,probs=0.975),
                   med = median(mu_py_seeds), 
                   hpdi.lo95 = HPDI(mu_py_seeds,.95)[1], 
                   hpdi.hi95 = HPDI(mu_py_seeds,.95)[2]
  )

seeds_py<-seeds_py %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(site,year)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::select(site,year,mu,ci.lo95,ci.hi95,med,hpdi.lo95,hpdi.hi95)

#saveRDS(seeds_py,file=paste0(dirPars,"seedsSummary.RDS"))

write.csv(seeds_py,file=paste0(dirPars,"seedsSummary.csv"))