# -------------------------------------------------------------------
# Analysis of fitness models
# Outputs reproductive success estimates
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE,max.print=100000)
# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(rjags) # jags interface
library(MCMCvis)
library(tidyverse)
library(reshape2)
library(HDInterval)
library(bayesplot)
library(coda)
library(rethinking)

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/"
sampleFiles <- paste0(directory,list.files(directory))

# -------------------------------------------------------------------
# Seedling survival to fruiting
# -------------------------------------------------------------------
sigmaMCMCsamples = readRDS(sampleFiles[grep("sigma.RDS",sampleFiles)])
sigmaMCMCsamples.noPool = readRDS(sampleFiles[grep("sigma-noPool.RDS",sampleFiles)])

# -------------------------------------------------------------------
# Fruits per plant
# -------------------------------------------------------------------

observedTfeMCMCsamples <- readRDS(sampleFiles[grep("observedTFE.RDS",sampleFiles)])
compositeTfeMCMCsamples <- readRDS(sampleFiles[grep("compositeTFE.RDS",sampleFiles)])
observedTfeMCMCsamples.noPool <- readRDS(sampleFiles[grep("observedTFE-noPool.RDS",sampleFiles)])
compositeTfeMCMCsamples.noPool <- readRDS(sampleFiles[grep("compositeTFE-noPool.RDS",sampleFiles)])

# -------------------------------------------------------------------
# Seeds per fruit
# -------------------------------------------------------------------

observedSeedsMCMCsamples <- readRDS(sampleFiles[grep("seedsPerFruit.RDS",sampleFiles)])
observedSeedsMCMCsamples.noPool <- readRDS(sampleFiles[grep("seedsPerFruit-noPool.RDS",sampleFiles)])

################################################################################
# Create data frames for analysis
#################################################################################

sigmaMCMCsamples = sigmaMCMCsamples
tfeMCMCsamples = cbind(observedTfeMCMCsamples,compositeTfeMCMCsamples)

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/"
saveRDS(boot::inv.logit(sigmaMCMCsamples),file="~/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/sigma-analysis.RDS")
saveRDS(tfeMCMCsamples,file="~/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/tfe-analysis.RDS")
saveRDS(observedSeedsMCMCsamples,file="~/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/phi-analysis.RDS")

# check to make sure data frames are the same size
# this is not automatic at the moment
dim(sigmaMCMCsamples)==dim(tfeMCMCsamples)
dim(sigmaMCMCsamples)==dim(observedSeedsMCMCsamples)
dim(tfeMCMCsamples)==dim(observedSeedsMCMCsamples)

# -------------------------------------------------------------------
# Calculate reproductive success
# -------------------------------------------------------------------
rsMCMCsamples = boot::inv.logit(sigmaMCMCsamples)*tfeMCMCsamples*observedSeedsMCMCsamples
saveRDS(rsMCMCsamples,file="~/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/rsPosterior.RDS")

seedsPerPlant = tfeMCMCsamples*observedSeedsMCMCsamples
saveRDS(seedsPerPlant,file="~/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/seedsPerPlantPosterior.RDS")



################################################################################
# Create data frames for analysis: no pool
#################################################################################

sigmaMCMCsamples.noPool = sigmaMCMCsamples.noPool
tfeMCMCsamples.noPool = cbind(observedTfeMCMCsamples.noPool,compositeTfeMCMCsamples.noPool)

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/"
saveRDS(boot::inv.logit(sigmaMCMCsamples.noPool),file="~/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/sigma-analysis-noPool.RDS")
saveRDS(tfeMCMCsamples.noPool,file="~/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/tfe-analysis-noPool.RDS")
saveRDS(observedSeedsMCMCsamples.noPool,file="~/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/phi-analysis-noPool.RDS")

# check to make sure data frames are the same size
# this is not automatic at the moment
dim(sigmaMCMCsamples.noPool)==dim(tfeMCMCsamples.noPool)
dim(sigmaMCMCsamples.noPool)==dim(observedSeedsMCMCsamples.noPool)
dim(tfeMCMCsamples.noPool)==dim(observedSeedsMCMCsamples.noPool)

# -------------------------------------------------------------------
# Calculate reproductive success
# -------------------------------------------------------------------
rsMCMCsamples.noPool = boot::inv.logit(sigmaMCMCsamples.noPool)*tfeMCMCsamples.noPool*observedSeedsMCMCsamples.noPool
saveRDS(rsMCMCsamples.noPool,file="~/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/rsPosterior-noPool.RDS")

seedsPerPlant.noPool = tfeMCMCsamples.noPool*observedSeedsMCMCsamples.noPool
saveRDS(seedsPerPlant.noPool,file="~/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/seedsPerPlantPosterior-noPool.RDS")


