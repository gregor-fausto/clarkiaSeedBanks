#################################################################################
################################################################################
################################################################################
# Code for seeds per fruit
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
mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("seedsPerFruitSamples.rds",mcmcSampleDirectory)]])
data <- readRDS(mcmcSampleDirectory[[grep("seedsPerFruitData.rds",mcmcSampleDirectory)]])

################################################################################
# Seeds per undamaged fruit
#################################################################################
mu_seeds=exp(MCMCchains(mcmcSamples, params=c("mu.log_seeds")))

saveRDS(mu_seeds,"/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/seedsPerFruit.RDS")