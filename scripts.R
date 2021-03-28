# -------------------------------------------------------------------
# Run scripts for model fitting
# -------------------------------------------------------------------
n.adapt = 3000
n.update = 5000
n.iterations = 15000
n.thin = 1

dataDirectory = "/Users/Gregor/Desktop/postProcessingData-2021/"
modelDirectory = "/Users/Gregor/Desktop/jagsScripts/"
fileDirectory = "/Users/Gregor/Desktop/mcmcSamples/"
scriptDirectory = "/Users/Gregor/Desktop/scripts/"

# run models with partial pooling
source(paste0(scriptDirectory,"modelScripts-viability-parallel.R"))
Sys.sleep(60*10)
source(paste0(scriptDirectory,"modelScripts-seeds-parallel.R"))
Sys.sleep(60*10)
source(paste0(scriptDirectory,"modelScripts-seedlingSurvival-parallel.R"))
Sys.sleep(60*10)
source(paste0(scriptDirectory,"modelScripts-seedsperfruit-parallel.R"))
Sys.sleep(60*10)
source(paste0(scriptDirectory,"modelScripts-fruits-parallel.R"))

# run models with no pooling
source(paste0(scriptDirectory,"modelScripts-seedlingSurvival-parallel-noPool.R"))
Sys.sleep(60*10)
source(paste0(scriptDirectory,"modelScripts-seedsperfruit-parallel-noPool.R"))
Sys.sleep(60*10)
source(paste0(scriptDirectory,"modelScripts-fruits-parallel-noPool.R"))

# -------------------------------------------------------------------
# Run scripts for evaluating convergence 
# -------------------------------------------------------------------
# scriptConvergenceDirectory = "/Users/Gregor/Desktop/scriptsModelConvergence/"
# fileDirectory = "/Users/Gregor/Desktop/mcmcSamples/"
# outputDirectory = "/Users/Gregor/Desktop/convergence/"

scriptConvergenceDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsModelConvergence/"
fileDirectory = "/Users/Gregor/Dropbox/dataLibrary/mcmcSamplesThinned/"
outputDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/convergence/"

# check convergence for belowground models
source(paste0(scriptConvergenceDirectory,"modelConvergenceViability.R"))
source(paste0(scriptConvergenceDirectory,"modelConvergenceBelowground.R"))

# check convergence for partially pooled aboveground models
source(paste0(scriptConvergenceDirectory,"modelConvergenceSeedlingSurvival.R"))
source(paste0(scriptConvergenceDirectory,"modelConvergenceFruits.R"))
source(paste0(scriptConvergenceDirectory,"modelConvergenceSeedsperfruit.R"))

# check convergence for no pooling aboveground models
# source(paste0(scriptConvergenceDirectory,"modelConvergenceSeedlingSurvival-noPool.R"))
# source(paste0(scriptConvergenceDirectory,"modelConvergenceFruits-noPool.R"))
# source(paste0(scriptConvergenceDirectory,"modelConvergenceSeedsperfruit-noPool.R"))

# pause script here to check that convergence has been achieved

# -------------------------------------------------------------------
# Run scripts for model checks
# -------------------------------------------------------------------
scriptCheckDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecks/"
fileDirectory = "/Users/Gregor/Dropbox/dataLibrary/mcmcSamplesThinned/"
outputDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/modelChecks-2020/"

source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/modelChecksBelowground.R")
source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/modelChecksSurvival.R")
source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/modelChecksFruitsPerPlant.R")
source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/modelChecksSeedsPerFruit.R")

# following scripts not yet written
# source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/modelChecksViability.R")

# -------------------------------------------------------------------
# Run scripts to output summary tables and summary plots
# -------------------------------------------------------------------

# need to consider what the right scripts are for here

# -------------------------------------------------------------------
# Run scripts to recover structured model parameters
# -------------------------------------------------------------------

source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsSummaries/parametersBelowground-1.R")
source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsSummaries/parametersFruitsPerPlant.R")
source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsSummaries/parametersSeedsPerFruit.R")

# -------------------------------------------------------------------
# Run scripts to analyze data for paper
# -------------------------------------------------------------------

source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/scripts/calculateReproductiveSuccess.R")
source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/scripts/estimateCorrelationGerminationSurvival-Reanalysis.R")
source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/scripts/estimateCorrelationGerminationReproductiveSuccess-Reanalysis.R")
source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/scripts/simulatePopulationGrowthRateDIModel-Reanalysis.R")
