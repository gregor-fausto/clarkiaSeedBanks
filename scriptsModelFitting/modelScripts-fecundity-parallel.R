# -------------------------------------------------------------------
# Models for fitting data on fruits per plant and seeds per fruit
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
# rm(list=setdiff(ls(all=TRUE),c(dataDirectory,modelDirectory,fileDirectory,n.adapt,n.update,n.iterations,n.thin))) # if using in source(script)
options(stringsAsFactors = FALSE)
# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(rjags) 
library(tidybayes)
library(tidyverse)
library(parallel)

# -------------------------------------------------------------------
# Set directories
# -------------------------------------------------------------------
dataDirectory = "/Users/Gregor/Dropbox/dataLibrary/postProcessingData/"
modelDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsModelFitting/jagsScripts/"
fileDirectory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/parallel/"

# -------------------------------------------------------------------
# Import and organize data
# -------------------------------------------------------------------
countFruitsPerPlantAllPlots <- readRDS(paste0(dataDirectory,"countFruitsPerPlantAllPlots.RDS"))

countFruitsPerPlantAllPlots <- countFruitsPerPlantAllPlots %>%
  dplyr::rename(y_tfe = countFruitNumberPerPlant) %>%
  dplyr::select(site,year,y_tfe)

countFruitsPerPlantAllPlots$year <- as.character(countFruitsPerPlantAllPlots$year)

countUndamagedDamagedFruitsPerPlantAllPlots <- readRDS(paste0(dataDirectory,"countUndamagedDamagedFruitsPerPlantAllPlots.RDS"))

countUndamagedDamagedFruitsPerPlantAllPlots <- countUndamagedDamagedFruitsPerPlantAllPlots %>%
  dplyr::rename(y_und = countUndamagedFruitNumberPerPlant) %>%
  dplyr::rename(y_dam = countDamagedFruitNumberPerPlant) %>%
  dplyr::rename(site2 = site) %>%
  dplyr::rename(year2 = year) %>%
  dplyr::select(site2,year2,y_und,y_dam) 

countUndamagedDamagedFruitsPerPlantAllPlots$year2 <- as.character(countUndamagedDamagedFruitsPerPlantAllPlots$year2)

countSeedPerFruit <-readRDS(paste0(dataDirectory,"countSeedPerFruit.RDS"))

countSeedPerUndamagedFruit <- countSeedPerFruit %>%
  dplyr::filter(demography==1) %>%
  dplyr::filter(damaged==0) %>%
  dplyr::rename(site3 = site) %>%
  dplyr::rename(year3 = year) %>%
  dplyr::select(site3,year3,sdno)

countSeedPerUndamagedFruit$year3 <- as.character(countSeedPerUndamagedFruit$year3)

countSeedPerDamagedFruit <- countSeedPerFruit %>%
  dplyr::filter(demography==1) %>%
  dplyr::filter(damaged==1) %>%
  dplyr::rename(site4 = site) %>%
  dplyr::rename(year4 = year) %>%
  dplyr::rename(sdno_dam = sdno) %>%
  dplyr::select(site4,year4,sdno_dam)

countSeedPerDamagedFruit$year4 <- as.character(countSeedPerDamagedFruit$year4)

# -------------------------------------------------------------------
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------
data <- tidybayes::compose_data(countFruitsPerPlantAllPlots,
                                countUndamagedDamagedFruitsPerPlantAllPlots,
                                countSeedPerUndamagedFruit,
                                countSeedPerDamagedFruit)

data$n = dim(countFruitsPerPlantAllPlots)[1]
data$n2 = dim(countUndamagedDamagedFruitsPerPlantAllPlots)[1]
data$n3 = dim(countSeedPerUndamagedFruit)[1]
data$n4 = dim(countSeedPerDamagedFruit)[1]

detach("package:tidyverse", unload=TRUE)

# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# -------------------------------------------------------------------
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 3000
n.update = 5000
n.iterations = 10000
n.thin = 1

# -------------------------------------------------------------------
# Describe functions
# -------------------------------------------------------------------

initsMu0 <- function(samps = data$n_site2){
  rgamma(n = samps, shape = 1, rate = 1)
}

initsSigma0 <- function(samps = data$n_site){
  extraDistr::rhnorm(n = samps, sigma = 1)
}

initsSigma <- function(rows = data$n_site, cols = data$n_year){
  matrix(extraDistr::rhnorm(n = rows*cols, sigma = 1), rows, cols)
}

# # set inits for JAGS
#  inits <- list()
# for(i in 1:3){
#   inits[[i]] <- list(initsMu0(), initsMu0(), initsMu0(), initsMu0(), initsMu0(),
#                      initsSigma0(), initsSigma0(), initsSigma0(), initsSigma0(), initsSigma0(),
#                      initsSigma(rows = data$n_site2,cols=data$n_year),
#                      initsSigma(rows = data$n_site2,cols=data$n_year2),
#                      initsSigma(rows = data$n_site2,cols=data$n_year2),
#                      initsSigma(rows = data$n_site2,cols=data$n_year3),
#                      initsSigma(rows = data$n_site2,cols=data$n_year4) )
# 
#   names(inits[[i]]) = c(paste(rep("nu",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
#                         paste(rep("sigma0",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
#                         paste(rep("sigma",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"))
# 
# }
# 
# # Call to JAGS
# 
# # tuning (n.adapt)
# jm = jags.model("~/Dropbox/clarkiaSeedBanks/scriptsModelFitting/jagsScripts/jags-fecundity.R",
#                 data = data, inits = inits,
#                 n.chains = length(inits), n.adapt = n.adapt)

# -------------------------------------------------------------------
# Set up JAGS for parallel run
# -------------------------------------------------------------------
library(stringr)

cl <- makeCluster(3)

myWorkers <- NA
for(i in 1:3) myWorkers[i] <- stringr::word(capture.output(cl[[i]]), -1)

initsFun = function(){temp = list(initsMu0(), initsMu0(), initsMu0(), initsMu0(), initsMu0(),
                                  initsSigma0(), initsSigma0(), initsSigma0(), initsSigma0(), initsSigma0(),
                                  initsSigma(rows = data$n_site2,cols=data$n_year),
                                  initsSigma(rows = data$n_site2,cols=data$n_year2),
                                  initsSigma(rows = data$n_site2,cols=data$n_year2),
                                  initsSigma(rows = data$n_site2,cols=data$n_year3),
                                  initsSigma(rows = data$n_site2,cols=data$n_year4),
                                  "base::Mersenne-Twister", runif(1, 1, 2000))

names(temp) = c(paste(rep("nu",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                paste(rep("sigma0",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                paste(rep("sigma",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                ".RNG.name",".RNG.seed")
return(temp)
}

inits = list(initsFun(),initsFun(),initsFun())

parsToMonitor = c(paste(rep("nu",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                  paste(rep("sigma0",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                  paste(rep("mu",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
                  paste(rep("sigma",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"))
parsToCheck = c("y_tfe.sim","chi2.tfe.obs","chi2.tfe.sim",
                       "y_und.sim","chi2.und.obs","chi2.und.sim",
                       "y_dam.sim","chi2.dam.obs","chi2.dam.sim",
                       "y_sd.sim","chi2.sd.obs","chi2.sd.sim",
                       "y_sd_dam.sim","chi2.sd_dam.obs","chi2.sd_dam.sim")

parallel::clusterExport(cl, c("myWorkers","data", "inits", "n.adapt", "n.update",
                              "n.iterations", "parsToMonitor", "parsToCheck"))

out <- clusterEvalQ(cl, {
  library(rjags)
  jm = jags.model(file=paste0(modelDirectory,"jags-fecundity.R"),
                  data = data, n.chains = 1,
                  n.adapt = n.adapt, 
                  inits = inits[[which(myWorkers==Sys.getpid())]])
  update(jm, n.iter = n.update)
  zm = coda.samples(jm, variable.names = c(parsToMonitor,parsToCheck),
                    n.iter = n.iterations, thin = 1)
  return(as.mcmc(zm))
})
stopCluster(cl)
samples.rjags = mcmc.list(out)

dir.create(file.path(fileDirectory), showWarnings = FALSE)
saveRDS(samples.rjags,file=paste0(fileDirectory,"fecunditySamples.rds"))
saveRDS(data,file=paste0(fileDirectory,"fecundityData.rds"))
