# -------------------------------------------------------------------
# Checking individual models
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)
# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(rjags) # jags interface
library(MCMCvis)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(HDInterval)
library(bayesplot)
library(stringr)

set.seed(10)
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
setwd("~/Dropbox/clarkiaSeedBanks/dataForAnalysis")
load(file="sigmaDF.RData") 

setwd("~/Dropbox/projects/clarkiaScripts/data/cleanData")
load(file="fecDF.RData")

setwd("~/Dropbox/projects/clarkiaScripts/data/cleanData")
load(file="phiIndDF.RData")

### MULTIPLE SITE, MULTIPLE YEAR MODEL (POOLED)

siteFilter = as.character(unique(sigmaDF$site))[1:20]
yearFilter = 2007:2012

# need to exlude 2 observations if using all data points without latent states
sigmaDat <- sigmaDF %>% 
  dplyr::filter(year %in% yearFilter & site %in% siteFilter) %>%
  dplyr::select(c(year,site,transect,position,noSeedlings,noFruitingPlants)) %>%
  dplyr::mutate( p = noFruitingPlants/noSeedlings ) %>%
  dplyr::mutate( test = ifelse(p>1, 1, 0)) %>%
  dplyr::filter( test < 1 | is.na(test)) %>% 
  dplyr::filter( !is.na(noSeedlings)) %>%
  dplyr::select(c(site, year, transect, position, noSeedlings, noFruitingPlants)) %>%
  unique

# remove plots without plants (NA)
fecDat <- fecDF %>% 
  dplyr::filter(year %in% yearFilter & site %in% siteFilter) %>%
  dplyr::select(c(year,site,transect,position,totalFruitEquivalents)) %>%
  dplyr::filter( !is.na(totalFruitEquivalents) )

# remove plots without plants (NA)
phiDat <- phiIndDF %>% 
  dplyr::filter(year %in% yearFilter & site %in% siteFilter) %>%
  dplyr::select(c(year,site,response))

# index 
sigmaDat$site.index <- as.integer(as.factor(sigmaDat$site))
sigmaDat$year.index <- as.integer(as.factor(sigmaDat$year))

fecDat$site.index <- as.integer(as.factor(fecDat$site))
fecDat$year.index <- as.integer(as.factor(fecDat$year))

phiDat$site.index <- as.integer(as.factor(phiDat$site))
phiDat$year.index <- as.integer(as.factor(phiDat$year))

if(length(unique(sigmaDat$year.index)) != length(unique(fecDat$year.index))) stop("unequal years")
if(length(unique(sigmaDat$year.index)) != length(unique(phiDat$year.index))) stop("unequal years")

if(length(unique(sigmaDat$site.index)) != length(unique(fecDat$site.index))) stop("unequal sites")
if(length(unique(sigmaDat$site.index)) != length(unique(phiDat$site.index))) stop("unequal sites")


data = list(
  n = as.double(sigmaDat$noSeedlings),
  y = as.double(sigmaDat$noFruitingPlants),
  N = nrow(sigmaDat),
  siteSigma = as.double(sigmaDat$site.index),
  yearSigma = as.double(sigmaDat$year.index),
  y1 = as.double(fecDat$totalFruitEquivalents),
  N1 = nrow(fecDat),
  siteFec = as.double(fecDat$site.index),
  yearFec = as.double(fecDat$year.index),
  y2 = as.double(phiDat$response),
  N2 = nrow(phiDat),
  sitePhi = as.double(phiDat$site.index),
  yearPhi = as.double(phiDat$year.index),
  nyears = length(unique(sigmaDat$year.index)),
  nsites = length(unique(sigmaDat$site.index))
)

initsS = list(
  list( alphaS = rep(-10, data$nsites), 
        betaS = rep(-10, data$nyears), 
        gammaS = matrix(-10, nrow=data$nsites,ncol=data$nyears)),
  list( 
        alphaS = rep(0, data$nsites), 
        betaS = rep(0, data$nyears), 
        gammaS = matrix(0, nrow=data$nsites,ncol=data$nyears) ),
  list( 
        alphaS = rep(10, data$nsites), 
        betaS = rep(10, data$nyears), 
        gammaS = matrix(10, nrow=data$nsites,ncol=data$nyears))
)


initsF = list(
  list( 
        alphaF = rep(-10, data$nsites), 
        betaF = rep(-10, data$nyears), 
        gammaF = matrix(-10, nrow=data$nsites,ncol=data$nyears), 
        rF = matrix(10, nrow=data$nsites,ncol=data$nyears)),
  list( 
        alphaF = rep(0, data$nsites), 
        betaF = rep(0, data$nyears), 
        gammaF = matrix(0, nrow=data$nsites,ncol=data$nyears), 
        rF = matrix(10, nrow=data$nsites,ncol=data$nyears)),
  list( 
        alphaF = rep(10, data$nsites), 
        betaF = rep(10, data$nyears), 
        gammaF = matrix(10, nrow=data$nsites,ncol=data$nyears), 
        rF = matrix(10, nrow=data$nsites,ncol=data$nyears))
)


initsP = list(
  list( 
        alphaP = rep(-10, data$nsites),
        betaP = rep(-10, data$nyears),
        gammaP = matrix(-10, nrow=data$nsites,ncol=data$nyears),
        rP = matrix(10, nrow=data$nsites,ncol=data$nyears)),
  list(
        alphaP = rep(0, data$nsites),
        betaP = rep(0, data$nyears),
        gammaP = matrix(0, nrow=data$nsites,ncol=data$nyears),
        rP = matrix(10, nrow=data$nsites,ncol=data$nyears)),
  list( 
        alphaP = rep(10, data$nsites),
        betaP = rep(10, data$nyears),
        gammaP = matrix(10, nrow=data$nsites,ncol=data$nyears),
        rP = matrix(10, nrow=data$nsites,ncol=data$nyears))
)

# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 1000
n.update = 10000
n.iter = 10000

# Call to JAGS

# tuning (n.adapt)
jm = jags.model("~/Dropbox/modelsF2019/models/sigmaJAGS.R", data = data, inits = initsS,
                n.chains = length(initsS), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

intercepts = c("alphaS","alphaF","alphaP")
slopes = c("betaS","betaF","betaP")
variances = c("gammaS","gammaF","gammaP")
sims = c("y.sim","y1.sim","y2.sim")

# chain (n.iter)
zmP = coda.samples(jm, variable.names = c(intercepts,slopes,variances,sims), n.iter = n.iter)

save(zmP,file="~/Dropbox/modelsF2019/output/sigmaFit")

# tuning (n.adapt)
jm = jags.model("~/Dropbox/modelsF2019/models/fecJAGS.R", data = data, inits = initsF,
                n.chains = length(initsF), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

intercepts = c("alphaS","alphaF","alphaP")
slopes = c("betaS","betaF","betaP")
variances = c("gammaS","gammaF","gammaP")
sims = c("y.sim","y1.sim","y2.sim")

# chain (n.iter)
zmP = coda.samples(jm, variable.names = c(intercepts,slopes,variances,sims), n.iter = n.iter)

save(zmP,file="~/Dropbox/modelsF2019/output/fecFit")

# tuning (n.adapt)
jm = jags.model("~/Dropbox/modelsF2019/models/phiJAGS.R", data = data, inits = initsP,
                n.chains = length(initsP), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

intercepts = c("alphaS","alphaF","alphaP")
slopes = c("betaS","betaF","betaP")
variances = c("gammaS","gammaF","gammaP")
sims = c("y.sim","y1.sim","y2.sim")

# chain (n.iter)
zmP = coda.samples(jm, variable.names = c(intercepts,slopes,variances,sims), n.iter = n.iter)

save(zmP,file="~/Dropbox/modelsF2019/output/phiFit")


