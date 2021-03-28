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

# -------------------------------------------------------------------
# Fruits per plant
# -------------------------------------------------------------------

observedTfeMCMCsamples <- readRDS(sampleFiles[grep("observedTFE.RDS",sampleFiles)])
compositeTfeMCMCsamples <- readRDS(sampleFiles[grep("compositeTFE.RDS",sampleFiles)])

# -------------------------------------------------------------------
# Seeds per fruit
# -------------------------------------------------------------------

observedSeedsMCMCsamples <- readRDS(sampleFiles[grep("seedsPerFruit.RDS",sampleFiles)])

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

