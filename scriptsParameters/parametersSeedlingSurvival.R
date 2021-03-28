#################################################################################
################################################################################
################################################################################
# Code for figures to compare the following modeling approaches for the seedling survivorship data
# 
# Scripts by Gregor Siegmund
# fausto.siegmund@gmail.com
# last updated 03-07-2021
#################################################################################
#################################################################################
#################################################################################
# rm(list=ls(all=TRUE)) # clear R environment
rm(list=setdiff(ls(all=TRUE),c("fileDirectory","outputDirectory"))) # if using in source(script)
options(stringsAsFactors = FALSE)

library(MCMCvis)
library(tidybayes)
library(tidyverse)
library(magrittr)
library(bayesplot)
library(rethinking)

mcmcSampleDirectory <- paste0(fileDirectory,list.files(fileDirectory))
mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("seedlingSurvivalSamples.rds",mcmcSampleDirectory)]])
data <- readRDS(mcmcSampleDirectory[[grep("seedlingSurvivalData.rds",mcmcSampleDirectory)]])

################################################################################
# Data directory
#################################################################################

mu=MCMCchains(mcmcSamples,params="mu")
#mu.p = apply(mu,2,boot::inv.logit)

saveRDS(mu,"/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/sigma.RDS")