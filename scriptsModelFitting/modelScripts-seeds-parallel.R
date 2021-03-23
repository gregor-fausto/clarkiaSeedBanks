# -------------------------------------------------------------------
# Models for seed bag burial experiment and s0
# -------------------------------------------------------------------
# log likelihood
# https://sourceforge.net/p/mcmc-jags/discussion/610036/thread/8211df61/
# https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.43.903&rep=rep1&type=pdf
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
dataDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/"
modelDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsModelFitting/jagsScripts/"
fileDirectory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/parallel/"

# -------------------------------------------------------------------
# Import and organize seed bag burial data
# -------------------------------------------------------------------

data = readRDS(paste0(dataDirectory,"seedBagMaster.RDS"))

# rename variables and oragnize
data = data %>%
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="-", remove=FALSE) %>%
  tidyr::unite(col='siteBag', c(site,bagNo), sep="-", remove=FALSE) %>%
  dplyr::mutate(siteBag = as.factor(siteBag))   %>%
  dplyr::mutate(year = as.factor(yearStart),
                age = as.factor(age)) %>%
  dplyr::select(siteBag,site,year,age,totalJan,seedStart,seedlingJan,intactOct) %>%
  dplyr::rename(siteBags = site,
                yearBags = year,
                ageBags = age)

# create a reference dataframe
refIndex <- data.frame(name=rep(c("totalJan","intactOct"),3),
                       ageBags=as.factor(c(1,1,2,2,3,3)),
                       betaIndex=c(1,2,3,4,5,6),
                       months=c(3,12,15,24,27,36),
                       g_1 = c(0,1,1,1,1,1),
                       g_2 = c(0,0,0,1,1,1),
                       g_3 = c(0,0,0,0,0,1)) %>%
  dplyr::mutate(months=months/36)

# create a data frame for intact seed data
# pivot to storing data in long format indexed by response type
intactSeedData=data %>%
  dplyr::select(-seedlingJan) %>%
  tidyr::pivot_longer(cols=c(totalJan,intactOct)) %>%
  dplyr::left_join(refIndex,by=c("name","ageBags")) %>%
  dplyr::rename(response = name, y = value, siteSurvival = siteBags) %>%
  dplyr::mutate(gIndexSurvival = ageBags,
                yearSurvival = yearBags)

# create a reference data frame for event histories
gCompIndex <- data.frame(yearSurvival = as.factor(c(rep(2006,6),rep(2007,4),rep(2008,2))),
                         ageBags=as.factor(c(1,1,2,2,3,3,1,1,2,2,1,1)),
                         compIndex=c(1:12),
                         response=rep(c("totalJan","intactOct"),6)) 

# join intact seed data frame with event history data frame
intactSeedData<-intactSeedData %>% 
  dplyr::left_join(gCompIndex,by=c('yearSurvival','ageBags',"response"))

# create a data frame for germination data
germinationData=data %>%
  dplyr::select(-c(seedStart,intactOct)) %>%
  dplyr::mutate(gIndex=ageBags) %>%
  dplyr::rename(siteGermination = siteBags,
                yearGermination = yearBags) %>%
  dplyr::mutate(seedlingJan2=seedlingJan)

# indexing to avoid issue of dealing with ragged arrays in JAGS
germinationAge = c(1,2,3,1,2,1)
germinationYear = c(2006,2006,2006,2007,2007,2008)
germinationIndex = 1:6

df.index=data.frame(ageBags=as.factor(germinationAge),yearGermination=as.factor(germinationYear),germinationIndex=germinationIndex)
germinationData <- germinationData %>%
  # dplyr::mutate(yearGermination=as.integer(yearGermination),ageBags=as.numeric(ageBags)) %>%
  dplyr::left_join(df.index,by=c("ageBags","yearGermination")) %>%
  dplyr::mutate(yearGermination=as.factor(yearGermination),ageBags=as.factor(ageBags))
# -------------------------------------------------------------------
# Import and organize aboveground plot data on fecundity
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Set directories
# -------------------------------------------------------------------
dataDirectory = "/Users/Gregor/Dropbox/dataLibrary/postProcessingData/"

# counts of fruits per plot
censusSeedlingsFruitingPlants <- readRDS(paste0(dataDirectory,"censusSeedlingsFruitingPlants.RDS"))

censusSeedlings <- censusSeedlingsFruitingPlants %>%
  dplyr::filter(year%in%c(2008,2009)) %>%
  dplyr::mutate(yearRef=year-1) %>%
  dplyr::rename(t1=year) %>%
  dplyr::select(-fruitplNumber)

# counts of fruits per plant
countFruitsPerPlantTransects <- readRDS(paste0(dataDirectory,"countFruitsPerPlantTransects.RDS"))

countFruitsPerPlotTransects <- countFruitsPerPlantTransects %>%
  dplyr::filter(year%in%c(2007,2008)) %>%
  dplyr::group_by(site,year,transect,position,) %>%
  dplyr::summarise(totalFruitsPerPlot = sum(countFruitsPerPlant))

# counts of seeds per fruit
countSeedPerFruit <- readRDS(paste0(dataDirectory,"countSeedPerFruit.RDS"))

countSeedPerTotalFruitEquivalent <- countSeedPerFruit %>%
  dplyr::filter(demography==1) %>%
  dplyr::filter(year %in% c(2007,2008)) %>%
  dplyr::filter(damaged==0) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(meanSeedPerTotalFruitEquivalent = round(mean(sdno)))

seedRainPlot<-countFruitsPerPlotTransects %>%
  dplyr::left_join(countSeedPerTotalFruitEquivalent,by=c("site","year")) %>%
  dplyr::mutate(fec.est = totalFruitsPerPlot*meanSeedPerTotalFruitEquivalent) %>%
  dplyr::mutate(t=year) %>% 
  dplyr::rename(yearRef=year) %>%
  dplyr::mutate(site=as.character(site),transect=as.character(transect),position=as.character(position)) %>%
  dplyr::left_join(censusSeedlings,by=c("site","transect","position","yearRef")) %>%
  dplyr::mutate(p = seedlingNumber/fec.est)

#max(seedRainPlot$p)

seedRainTransect = seedRainPlot %>%
  dplyr::group_by(site,yearRef,transect) %>%
  dplyr::summarise(fec.est = sum(fec.est), seedlingNumber = sum(seedlingNumber)) %>%
  dplyr::mutate(p = seedlingNumber/fec.est)

seedRainTransectData = seedRainTransect %>% 
  dplyr::rename(sitePlot = site) %>%
  dplyr::rename(yearPlot = yearRef) %>%
  dplyr::rename(fec = fec.est) %>%
  dplyr::rename(plotSeedlings = seedlingNumber) %>%
  dplyr::select(-p) %>%
  dplyr::mutate(yearPlot = as.factor(yearPlot))

# indexing to avoid ragged array

df.index=data.frame(ageBags=as.factor(germinationAge),yearGermination=as.factor(germinationYear),germinationIndex=germinationIndex)
df.index=df.index %>%
  dplyr::filter(ageBags==1 &yearGermination==2007|yearGermination==2008) %>%
  dplyr::rename(yearPlot=yearGermination) %>%
  dplyr::rename(fecIndex = germinationIndex) %>%
  dplyr::select(yearPlot,fecIndex)
seedRainTransectData <- seedRainTransectData %>%
  dplyr::left_join(df.index,by=c("yearPlot")) 
# -------------------------------------------------------------------
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------


data <- tidybayes::compose_data(intactSeedData,germinationData,seedRainTransectData)

data$n1 = dim(intactSeedData)[1]
data$n2 = dim(germinationData)[1]
data$n3 = dim(seedRainTransectData)[1]

refDf<-data.frame(year=data$yearGermination,gIndex=data$gIndex) %>%
  unique %>%
  dplyr::mutate(germinationReference=1:6)

data$yearRefGerm = refDf$year
data$indexRefGerm = refDf$gIndex
data$refGerm = refDf$germinationReference

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

initsA <- function(samps = data$n_siteSurvival){
  rgamma(n = samps, shape = 2, rate = 2)
}

initsMu0 <- function(samps = data$n_siteSurvival){
  rnorm(n = samps, mean = 0, sd = 1)
}

initsMu0Mat <- function(rows = data$n_siteSurvival, cols = data$n_yearGermination){
  matrix(rnorm(n = rows*cols, mean = 0, sd = 1), rows, cols)
}

initsSigma0Weak <- function(samps = data$n_siteSurvival){
  extraDistr::rhnorm(n = samps, sigma = 1)
}

initsSigma0LessWeak <- function(samps = data$n_siteSurvival){
  extraDistr::rhnorm(n = samps, sigma = 2)
}

initsSigma0LessWeakMat <- function(rows = data$n_siteSurvival, cols = data$n_yearSurvival){
  matrix(extraDistr::rhnorm(n = rows*cols, sigma = 2), rows, cols)
}

initsSigmaWeak <- function(rows = data$n_siteSurvival, cols = data$n_yearSurvival){
  matrix(extraDistr::rhnorm(n = rows*cols, sigma = 1), rows, cols)
}

initsSigmaWeak2 <- function(rows = data$n_siteSurvival, cols = 6){
  matrix(extraDistr::rhnorm(n = rows*cols, sigma = 1), rows, cols)
}

initsSigmaWeak3 <- function(rows = data$n_siteSurvival, cols = 2){
  matrix(extraDistr::rhnorm(n = rows*cols, sigma = 1), rows, cols)
}

initsSigmaLessWeak <- function(rows = data$n_siteSurvival, cols = data$n_yearPlot){
  matrix(extraDistr::rhnorm(n = rows*cols, sigma = 2), rows, cols)
}

# set inits for JAGS
inits <- list()
for(i in 1:3){
  inits[[i]] <- list(initsA(), initsMu0(), initsSigma0Weak(), initsSigmaWeak(),
                     initsMu0Mat(), initsSigma0LessWeakMat(), initsSigmaWeak2(),
                     initsMu0(), initsSigma0LessWeak(), initsSigmaWeak3())

  names(inits[[i]]) = c("a","mu0_s","sigma0_s","sigma_s",
                        "mu0_g", "sigma0_g", "sigma_g",
                        "mu0_s0", "sigma0_s0", "sigma_s0")

}

# # Call to JAGS
#
# # tuning (n.adapt)
# jm = jags.model("~/Dropbox/clarkiaSeedBanks/scriptsModelFitting/jagsScripts/jags-seeds-parallel.R",
#                 data = data, inits = inits,
#                 n.chains = length(inits), n.adapt = n.adapt)

# need to reparameterize the sigma_g to get this to run in parallel to not deal with ragged array
# see viability data

# -------------------------------------------------------------------
# Set up JAGS for parallel run
# -------------------------------------------------------------------
library(stringr)

cl <- makeCluster(3)

myWorkers <- NA
for(i in 1:3) myWorkers[i] <- stringr::word(capture.output(cl[[i]]), -1)

initsFun = function(){temp = list(initsA(), initsMu0(), initsSigma0Weak(), initsSigmaWeak(),
                                  initsMu0Mat(), initsSigma0LessWeakMat(), initsSigmaWeak2(),
                                  initsMu0(), initsSigma0LessWeak(), initsSigmaWeak3(), 
                                  "base::Mersenne-Twister", runif(1, 1, 2000))

names(temp) = c("a","mu0_s","sigma0_s","sigma_s",
                "mu0_g", "sigma0_g", "sigma_g",
                "mu0_s0", "sigma0_s0", "sigma_s0",
                ".RNG.name",".RNG.seed")
return(temp)
}

inits = list(initsFun(),initsFun(),initsFun())

parsToMonitor = c("a","mu0_s","sigma0_s","sigma_s",
                  "mu0_g", "sigma0_g", "sigma_g",
                  "mu0_s0", "sigma0_s0", "sigma_s0")
parsToCheck = c("y_sim","chi2.yobs","chi2.ysim",
                "seedlingJan_sim","chi2.obs","chi2.sim",
                "plotSeedlings_sim","chi2.plot.obs","chi2.plot.sim")

parallel::clusterExport(cl, c("myWorkers","data", "inits", "n.adapt", "n.update",
                              "n.iterations", "parsToMonitor", "parsToCheck"))

out <- clusterEvalQ(cl, {
  library(rjags)
  jm = jags.model(file=paste0(modelDirectory,"jags-seeds.R"),
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
saveRDS(samples.rjags,file=paste0(fileDirectory,"seedSamples.rds"))
saveRDS(data,file=paste0(fileDirectory,"seedData.rds"))