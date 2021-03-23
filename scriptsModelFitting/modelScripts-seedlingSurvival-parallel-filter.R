# ------------------------------------------------------------------- 
# Models for seedling survival to fruiting
# -------------------------------------------------------------------
rm(list = ls(all = TRUE))  # clear R environment
options(stringsAsFactors = FALSE)
# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(rjags)
library(tidybayes)
library(tidyverse)
library(parallel)

# -------------------------------------------------------------------
# ------------------------------------------------------------------- 
#Import and organize data
# -------------------------------------------------------------------
# ------------------------------------------------------------------- 
# setwd and read data files
censusSeedlingsFruitingPlants <- readRDS("~/Dropbox/dataLibrary/postProcessingData/censusSeedlingsFruitingPlants.RDS")

# ------------------------------------------------------------------- 
# Issue with number of fruiting plants greater than number of seedlings
write.table(censusSeedlingsFruitingPlants %>% dplyr::filter(fruitplNumber > seedlingNumber), 
  "~/Dropbox/clarkiaSeedBanks/dataToCheck/seedlingSurvival.txt", sep = "\t", row.names = FALSE)
write(c("Issue with rows of data where number of fruiting plants is greater than number of seedlings"), 
  file = "~/Dropbox/clarkiaSeedBanks/dataToCheck/seedlingSurvival-metadata.txt")

## Issue with rows of data where number of seedlings is NA but number of fruiting
## plants is numeric; these are true missing data
write.table(censusSeedlingsFruitingPlants %>% dplyr::filter(is.na(seedlingNumber)), 
  "~/Dropbox/clarkiaSeedBanks/dataToCheck/seedlingSurvival.txt", sep = "\t", row.names = FALSE, 
  append = TRUE, col.names = FALSE)
write(c("Rows of data where number of seedlings is NA but number of fruiting plants is numeric; these are true missing data from 2006"), 
  file = "~/Dropbox/clarkiaSeedBanks/dataToCheck/seedlingSurvival-metadata.txt", 
  append = TRUE)

## Issue with rows of data where number of seedlings is numeric but number of
## fruiting plants is NA
write.table(censusSeedlingsFruitingPlants %>% dplyr::filter(is.na(fruitplNumber)), 
  "~/Dropbox/clarkiaSeedBanks/dataToCheck/seedlingSurvival.txt", sep = "\t", row.names = FALSE, 
  append = TRUE, col.names = FALSE)
write(c("Issue with rows of data where number of seedlings is numeric but number of fruiting plants is NA"), 
  file = "~/Dropbox/clarkiaSeedBanks/dataToCheck/seedlingSurvival-metadata.txt", 
  append = TRUE)


# filter out all data with issues
censusSeedlingsFruitingPlants <- censusSeedlingsFruitingPlants %>% # more fruiting plants than seedlings
dplyr::filter(!(fruitplNumber > seedlingNumber)) %>% # NA for seedlings
dplyr::filter(!is.na(seedlingNumber)) %>% # NA for fruiting plants
dplyr::filter(!is.na(fruitplNumber))
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# ------------------------------------------------------------------- 
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------
# -------------------------------------------------------------------
censusSeedlingsFruitingPlants$year <- as.character(censusSeedlingsFruitingPlants$year)

data <- tidybayes::compose_data(censusSeedlingsFruitingPlants)

detach("package:tidyverse", unload = TRUE)

# -------------------------------------------------------------------
# ------------------------------------------------------------------- 
# Set JAGS parameters and random seed
# -------------------------------------------------------------------
# ------------------------------------------------------------------- 
# scalars that specify the number of iterations in the chain for adaptation number of
# iterations for burn-in number of samples in the final chain
n.adapt = 3000
n.update = 5000
n.iterations = 10000
n.thin = 1

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Describe functions 
# -------------------------------------------------------------------
# -------------------------------------------------------------------

initsMu0 <- function(samps = data$n_site) {
  rnorm(n = samps, mean = 0, sd = 1)
}

initsSigma0 <- function(samps = data$n_site) {
  extraDistr::rhnorm(n = samps, sigma = 1)
}

initsSigma <- function(rows = data$n_site, cols = data$n_year) {
  matrix(extraDistr::rhnorm(n = rows * cols, sigma = 1), rows, cols)
}

# # set inits for JAGS inits <- list() for(i in 1:3){ inits[[i]] <-
# list(initsMu0(), initsSigma0(), initsSigma() ) names(inits[[i]]) =
# c('mu0','sigma0','sigma') } # # Call to JAGS # # # tuning (n.adapt) jm =
# jags.model('~/Dropbox/clarkiaSeedBanks/scriptsModelFitting/jagsScripts/jags-seedlingSurvival.R',
# data = data, inits = inits, n.chains = length(inits), n.adapt = n.adapt)

# ------------------------------------------------------------------- Set up JAGS
# for parallel run
# -------------------------------------------------------------------
library(stringr)

cl <- makeCluster(3)

myWorkers <- NA
for (i in 1:3) myWorkers[i] <- stringr::word(capture.output(cl[[i]]), -1)

initsFun = function() {
  temp = list(initsMu0(), initsSigma0(), initsSigma(), "base::Mersenne-Twister", 
    runif(1, 1, 2000))
  
  names(temp) = c("mu0", "sigma0", "sigma", ".RNG.name", ".RNG.seed")
  return(temp)
}

inits = list(initsFun(), initsFun(), initsFun())

parsToMonitor = c("mu0", "sigma0", "mu", "sigma")
parsToCheck = c("fruitplNumber_sim", "chi2.obs", "chi2.sim")

parallel::clusterExport(cl, c("myWorkers", "data", "inits", "n.adapt", "n.update", 
  "n.iterations", "parsToMonitor", "parsToCheck"))

out <- clusterEvalQ(cl, {
  library(rjags)
  jm = jags.model(file = "~/Dropbox/clarkiaSeedBanks/scriptsModelFitting/jagsScripts/jags-seedlingSurvival.R", 
    data = data, n.chains = 1, n.adapt = n.adapt, inits = inits[[which(myWorkers == 
      Sys.getpid())]])
  update(jm, n.iter = n.update)
  zm = coda.samples(jm, variable.names = c(parsToMonitor, parsToCheck), n.iter = n.iterations, 
    thin = 1)
  return(as.mcmc(zm))
})
stopCluster(cl)
samples.rjags = mcmc.list(out)

fileDirectory <- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/parallel/")
dir.create(file.path(fileDirectory), showWarnings = FALSE)

saveRDS(samples.rjags, file = paste0(fileDirectory, "seedlingSurvivalSamples-filter.rds"))
saveRDS(data, file = paste0(fileDirectory, "seedlingData-filter.rds"))
