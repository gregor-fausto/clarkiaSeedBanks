# -------------------------------------------------------------------
# Run scripts for model fitting
# -------------------------------------------------------------------
n.adapt = 3000
n.update = 5000
n.iterations = 15000
n.thin = 1

dataDirectory = "/Users/Gregor/Desktop/data/"
modelDirectory = "/Users/Gregor/Desktop/jagsScripts/"
fileDirectory = "/Users/Gregor/Desktop/mcmcSamples/"
scriptDirectory = "/Users/Gregor/Desktop/scripts/"

source(paste0(scriptDirectory,"modelScripts-seeds-parallel.R"))
source(paste0(scriptDirectory,"modelScripts-viability-parallel.R"))
source(paste0(scriptDirectory,"modelScripts-seedlingSurvival-parallel"))
source(paste0(scriptDirectory,"modelScripts-fecundity-parallel.R"))

# -------------------------------------------------------------------
# Run scripts for evaluating convergence and model checking
# -------------------------------------------------------------------

source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/modelChecksBelowground.R")
source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/modelChecksSurvival.R")

# following scripts not yet written
# source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/modelChecksFruitsPerPlant.R")
# source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/modelChecksSeedsPerFruit.R")

# pause script here to check that convergence has been achieved

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
