# -------------------------------------------------------------------
# Models for seedling survival to fruiting
# -------------------------------------------------------------------
# rm(list=ls(all=TRUE)) # clear R environment
rm(list=setdiff(ls(all=TRUE),c("dataDirectory","modelDirectory","fileDirectory","n.adapt","n.update","n.iterations","n.thin"))) # if using in source(script)
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
# dataDirectory = "/Users/Gregor/Dropbox/dataLibrary/postProcessingData/"
# modelDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsModelFitting/jagsScripts/"
# fileDirectory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/parallel/"

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# setwd and read data files
censusSeedlingsFruitingPlants <- readRDS(paste0(dataDirectory,"censusSeedlingsFruitingPlants.RDS"))

# -------------------------------------------------------------------
# Recode and filter data
# -------------------------------------------------------------------
censusSeedlingsFruitingPlants <- censusSeedlingsFruitingPlants %>% 
  # more fruiting plants than seedlings
  # recode s.t. number of seedlings equals number of fruiting plants
  dplyr::mutate(seedlingNumber=ifelse(fruitplNumber>seedlingNumber,fruitplNumber,seedlingNumber)) %>% 
  # NA for seedlings
  # filter these out; missing counts
  dplyr::filter(!is.na(seedlingNumber)) %>% 
  # NA for fruiting plants
  #  filter these out; missing response
  dplyr::filter(!is.na(fruitplNumber)) 
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------
# -------------------------------------------------------------------
censusSeedlingsFruitingPlants$year <- as.character(censusSeedlingsFruitingPlants$year)

data <- tidybayes::compose_data(censusSeedlingsFruitingPlants)

detach("package:tidyverse", unload=TRUE)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# n.adapt = 3000
# n.update = 5000
# n.iterations = 15000
# n.thin = 1

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Describe functions
# -------------------------------------------------------------------
# -------------------------------------------------------------------

initsMu <- function(rows = data$n_site, cols = data$n_year){
  matrix(rnorm(n = rows*cols, mean = 0, sd = 1), rows, cols)
}

initsSigma <- function(rows = data$n_site, cols = data$n_year){
  matrix(extraDistr::rhnorm(n = rows*cols, sigma = 1), rows, cols)
}

# -------------------------------------------------------------------
# Set up JAGS for parallel run
# -------------------------------------------------------------------
library(stringr)

cl <- makeCluster(3)

myWorkers <- NA
for(i in 1:3) myWorkers[i] <- stringr::word(capture.output(cl[[i]]), -1)

initsFun = function(){temp = list(initsMu(),  initsSigma(), 
                                  "base::Mersenne-Twister", runif(1, 1, 2000))

names(temp) = c("mu",
                "sigma",
                ".RNG.name",".RNG.seed")
return(temp)
}

inits = list(initsFun(),initsFun(),initsFun())

parsToMonitor = c("mu","sigma")
parsToCheck = c("fruitplNumber_sim","chi2.obs","chi2.sim")

parallel::clusterExport(cl, c("myWorkers","data", "inits", "n.adapt", "n.update",
                              "n.iterations", "parsToMonitor", "parsToCheck","modelDirectory"))

out <- clusterEvalQ(cl, {
  library(rjags)
  jm = jags.model(file=paste0(modelDirectory,"jags-seedlingSurvival-noPool.R"),
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
saveRDS(samples.rjags,file=paste0(fileDirectory,"seedlingSurvivalSamples-noPool.rds"))
saveRDS(data,file=paste0(fileDirectory,"seedlingSurvivalData-noPool.rds"))
