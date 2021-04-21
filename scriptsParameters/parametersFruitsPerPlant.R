#################################################################################
################################################################################
################################################################################
# Code for figures of total fruit equivalents per plant
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
mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("fruitsSamples.rds",mcmcSampleDirectory)]])
data <- readRDS(mcmcSampleDirectory[[grep("fruitsData.rds",mcmcSampleDirectory)]])

mcmcSamplesSeeds <- readRDS(mcmcSampleDirectory[[grep("seedsPerFruitSamples.rds",mcmcSampleDirectory)]])
dataSeeds <- readRDS(mcmcSampleDirectory[[grep("seedsPerFruitData.rds",mcmcSampleDirectory)]])

################################################################################
# Total fruit equivalents
#################################################################################
mu_tfe=exp(MCMCchains(mcmcSamples, params=c("mu.log_tfe")))
saveRDS(mu_tfe,"/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/observedTFE.RDS")

################################################################################
# Create composite
#################################################################################

mu_und=exp(MCMCchains(mcmcSamples, params=c("mu.log_und")))
mu_dam=exp(MCMCchains(mcmcSamples, params=c("mu.log_dam")))

mu_seeds=exp(MCMCchains(mcmcSamplesSeeds, params=c("mu.log_seeds")))
mu_dam_seeds=exp(MCMCchains(mcmcSamplesSeeds, params=c("mu.log_dam_seeds")))

index = c()
for(i in 8:15){
  index=c(index,grep(paste0(",",i,"\\]"),colnames(mu_seeds)))
}
mu_seeds[,index] %>% dim
mu_dam_seeds %>% dim

ratio.seeds=mu_dam_seeds/mu_seeds[,index] 
mu_composite=mu_und + mu_dam*ratio.seeds
saveRDS(mu_composite,"/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/compositeTFE.RDS")


################################################################################
# No pool
#################################################################################

mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("fruitsSamples-noPool.rds",mcmcSampleDirectory)]])
mcmcSamplesSeeds <- readRDS(mcmcSampleDirectory[[grep("seedsPerFruitSamples-noPool.rds",mcmcSampleDirectory)]])

################################################################################
# Total fruit equivalents
#################################################################################
mu_tfe=exp(MCMCchains(mcmcSamples, params=c("mu.log_tfe")))
saveRDS(mu_tfe,"/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/observedTFE-noPool.RDS")

################################################################################
# Create composite
#################################################################################

mu_und=exp(MCMCchains(mcmcSamples, params=c("mu.log_und")))
mu_dam=exp(MCMCchains(mcmcSamples, params=c("mu.log_dam")))

mu_seeds=exp(MCMCchains(mcmcSamplesSeeds, params=c("mu.log_seeds")))
mu_dam_seeds=exp(MCMCchains(mcmcSamplesSeeds, params=c("mu.log_dam_seeds")))

index = c()
for(i in 8:15){
  index=c(index,grep(paste0(",",i,"\\]"),colnames(mu_seeds)))
}
mu_seeds[,index] %>% dim
mu_dam_seeds %>% dim

ratio.seeds=mu_dam_seeds/mu_seeds[,index] 
mu_composite=mu_und + mu_dam*ratio.seeds
saveRDS(mu_composite,"/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/compositeTFE-noPool.RDS")
