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

mcmcSamples <- readRDS(posteriorFiles[[9]])

dirPars = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/parameterSummary/"
dir.create(file.path(dirPars), showWarnings = FALSE)

################################################################################
# Germination
#################################################################################

directory = "/Users/Gregor/Dropbox/dataLibrary/workflow/tidyData/"
simFiles <- paste0(directory,list.files(directory))

data <- readRDS(simFiles[[1]])
siteNames = unique(data$siteBags)

################################################################################
# Germination 1
#################################################################################

g1_p <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(g1[site]) %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(mu = mean(g1),
                   ci.lo95 = quantile(g1,probs=0.025), 
                   ci.hi95 = quantile(g1,probs=0.975),
                   med = median(g1), 
                   hpdi.lo95 = HPDI(g1,.95)[1], 
                   hpdi.hi95 = HPDI(g1,.95)[2]
  )

g1_p<-cbind(siteNames=siteNames,g1_p) %>%
  dplyr::select(-site) %>%
  dplyr::rename(site = siteNames)



################################################################################
# Winter seed survival (s1)
#################################################################################

s1_p <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(s1[site]) %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(mu = mean(s1),
                   ci.lo95 = quantile(s1,probs=0.025), 
                   ci.hi95 = quantile(s1,probs=0.975),
                   med = median(s1), 
                   hpdi.lo95 = HPDI(s1,.95)[1], 
                   hpdi.hi95 = HPDI(s1,.95)[2]
  )

s1_p<-cbind(siteNames=siteNames,s1_p) %>%
  dplyr::select(-site) %>%
  dplyr::rename(site = siteNames)


# ################################################################################
# # Summer seed survival (s2)
# #################################################################################

s2_p <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(s2[site]) %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(mu = mean(s2),
                   ci.lo95 = quantile(s2,probs=0.025), 
                   ci.hi95 = quantile(s2,probs=0.975),
                   med = median(s2), 
                   hpdi.lo95 = HPDI(s2,.95)[1], 
                   hpdi.hi95 = HPDI(s2,.95)[2]
  )

s2_p<-cbind(siteNames=siteNames,s2_p) %>%
  dplyr::select(-site) %>%
  dplyr::rename(site = siteNames)



# ################################################################################
# # Winter seed survival, year 2 (s3)
# #################################################################################

s3_p <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(s3[site]) %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(mu = mean(s3),
                   ci.lo95 = quantile(s3,probs=0.025), 
                   ci.hi95 = quantile(s3,probs=0.975),
                   med = median(s3), 
                   hpdi.lo95 = HPDI(s3,.95)[1], 
                   hpdi.hi95 = HPDI(s3,.95)[2]
  )

s3_p<-cbind(siteNames=siteNames,s3_p) %>%
  dplyr::select(-site) %>%
  dplyr::rename(site = siteNames)


# saveRDS(g1_p,file=paste0(dirPars,"g1Summary.RDS"))
# saveRDS(s1_p,file=paste0(dirPars,"s1Summary.RDS"))
# saveRDS(s2_p,file=paste0(dirPars,"s2Summary.RDS"))
# saveRDS(s3_p,file=paste0(dirPars,"s3Summary.RDS"))

write.csv(g1_p,file=paste0(dirPars,"g1Summary.csv"))
write.csv(s1_p,file=paste0(dirPars,"s1Summary.csv"))
write.csv(s2_p,file=paste0(dirPars,"s2Summary.csv"))
write.csv(s3_p,file=paste0(dirPars,"s3Summary.csv"))
