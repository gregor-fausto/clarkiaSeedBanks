#################################################################################
################################################################################
################################################################################
# Code for comparing the following modeling approaches for the seedling survivorship data
# 1. binomial likelihood
# 2. binomial likelihood with beta prior, complete pooling
# 3. binomial likelihood with beta prior, partial pooling, parameterize mode
# 4. binomial likelihood with beta prior, partial pooling, parameterize mean
# 4. binomial likelihood with logit parameterization
#
# Scripts by Gregor Siegmund
# fausto.siegmund@gmail.com
# last updated 02-29-2020
#################################################################################
#################################################################################
#################################################################################

################################################################################
# Load packages
################################################################################
library(janitor)
library(tidyverse)
library(tidybayes)
library(rjags)
library(magrittr) # for %<>% with recover_types
library(dplyr)
library(MCMCvis)

#################################################################################
# get paths for JAGS script and figure directories
#################################################################################
dirJagsScripts = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/appendixOneScripts/")

################################################################################
# Load survivorship data
################################################################################
survivorshipData0615 <-
  readxl::read_excel(path = "~/Dropbox/dataLibrary/Clarkia-LTREB/20_demography_sites/Survivorship & Fecundity_06-15.xls")

survData <- survivorshipData0615 %>%
  janitor::clean_names(case = "lower_camel")

survVariables <- names(survData)
seedlingNumNames <-
  tidyselect::vars_select(survVariables, contains("seedlingNumber"))
fruitingPlantNumNames <-
  tidyselect::vars_select(survVariables, contains("fruitplNumber6"))

survDataSubset <- survData %>%
  dplyr::select(site,
                transect,
                position,
                seedlingNumNames,
                fruitingPlantNumNames)

seedlingData <- survDataSubset %>%
  dplyr::select(site, transect, position, seedlingNumNames) %>%
  tidyr::pivot_longer(cols = seedlingNumNames,
                      names_to = "year",
                      values_to = "seedlingNumber") %>%
  tidyr::separate(year, into = c("discard", "year")) %>%
  dplyr::mutate(year = as.numeric(paste0(20, year))) %>%
  dplyr::select(-discard)

fruitingplantData <- survDataSubset %>%
  dplyr::select(site, transect, position, fruitingPlantNumNames) %>%
  tidyr::pivot_longer(cols = fruitingPlantNumNames,
                      names_to = "year",
                      values_to = "fruitingPlantNumber") %>%
  tidyr::separate(year, into = c("discard", "year")) %>%
  dplyr::mutate(year = as.numeric(paste0(20, year))) %>%
  dplyr::select(-discard)

survDataReformatted <- seedlingData %>%
  dplyr::left_join(fruitingplantData, by = c("site", "transect", "position", "year")) %>%
  dplyr::rename(plot = position)

################################################################################
# Load site data
################################################################################
siteData <- read.csv(file="/Users/Gregor/Dropbox/dataLibrary/Datafile_Site_Environment_corr.csv",header=TRUE)

siteEastingData<-siteData %>%
  janitor::clean_names(case="lower_camel") %>%
  dplyr::select(site,easting) %>% unique %>%
  dplyr::filter(site %in% unique(survData$site))

################################################################################
# Prepare data for analysis
################################################################################
# filter dataset to rows in which
# the seedling number > 0 AND there are no NAs in the fruiting plant count
survDataAnalysis <- survDataReformatted %>%
    dplyr::filter(seedlingNumber>0 & !is.na(fruitingPlantNumber))

survDataAnalysis<-survDataAnalysis %>%
  dplyr::mutate(seedlingNumber = ifelse(seedlingNumber<fruitingPlantNumber,fruitingPlantNumber,seedlingNumber))

################################################################################
################################################################################
# Maximum likelihood estimate of survivorship, binomial likelihood
# by minimizing the negative log likelihood
################################################################################
################################################################################

# split the data frame into a list of lists by site and year
survDataAnalysisList <- split(survDataAnalysis,
                              list(survDataAnalysis$site,survDataAnalysis$year))

# create list (length 10, for 10 years) of lists (each length 20, for 20 sites)
survDataAnalysisListList <- list(survDataAnalysisList[1:20],
               survDataAnalysisList[21:40],
               survDataAnalysisList[41:60],
               survDataAnalysisList[61:80],
               survDataAnalysisList[81:100],
               survDataAnalysisList[101:120],
               survDataAnalysisList[121:140],
               survDataAnalysisList[141:160],
               survDataAnalysisList[161:180],
               survDataAnalysisList[181:200])

# empty matrix for maximum likelihood estimates
mles <- matrix(NA,20,10)
# empty matrix for number of trials
nTrials <- matrix(NA,20,10)

# obtain maximum likelihood estimates by minimizing the negative log likelihood for each year and site
for(i in 1:10){
  tmp.list <- survDataAnalysisListList[[i]]
  for(k in 1:20){
    tmp <- tmp.list[[k]]
    y = tmp$fruitingPlantNumber
    n = tmp$seedlingNumber
    negLogLik <- function(p,x=tmp$fruitingPlantNumber,size=tmp$seedlingNumber) -sum(dbinom(x=x, size=size, prob=p,log=TRUE))
    opt<-optimise(negLogLik,interval=c(0,1),maximum=FALSE)
    mles[k,i] <- opt$minimum
    nTrials[k,i] <- sum(tmp$seedlingNumber)
  }
}

# place maximum likelihood estimates in data frame
mleDataWide <- data.frame(site=unique(survDataAnalysis$site), mles )

names(mleDataWide)=c("site",2006:2015)

mleDataLong <- mleDataWide %>%
  tidyr::pivot_longer(cols=`2006`:`2015`,names_to="year",values_to="pHat")
mleDataLong$year<-as.numeric(mleDataLong$year)

saveRDS(mleDataLong,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/survivorshipMaxLikEst.rds")

################################################################################
################################################################################
# Prepare data for fit with JAGS using tidybayes
################################################################################
################################################################################

survDataAnalysisBayes <- survDataAnalysis %>% 
  dplyr::select(-c(transect,plot)) %>%
  dplyr::mutate(year = as.factor(year))

survDataAnalysisBayesJags <- tidybayes::compose_data(survDataAnalysisBayes)

saveRDS(survDataAnalysisBayes, file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/survivorshipDataBayes.rds")
saveRDS(survDataAnalysis,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/survivorshipDataAnalysis.rds")


################################################################################
################################################################################
# JAGS fit for model with binomial likelihood, beta prior, complete pooling
################################################################################
################################################################################

################################################################################
# Initial values
################################################################################

nSites = survDataAnalysisBayesJags$n_site
nYears = survDataAnalysisBayesJags$n_year

inits = list(list(p = matrix(rep(.1,nSites*nYears),nrow=nSites,ncol=nYears)),
             list(p = matrix(rep(.5,nSites*nYears),nrow=nSites,ncol=nYears)),
             list(p = matrix(rep(.9,nSites*nYears),nrow=nSites,ncol=nYears))) 

################################################################################
# Set parameters to monitor
################################################################################

parsToMonitor = c("p")

################################################################################
# Set JAGS parameters
################################################################################

# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 500
n.update = 5000
n.iterations = 1000
n.thin = 1

################################################################################
# Fit model 
################################################################################

# set random seed
set.seed(2)

# tuning (n.adapt)
jm = jags.model(paste0(dirJagsScripts,"survivorshipModel-binLik-betaPrior-completePooling.R"), 
                 data = survDataAnalysisBayesJags, inits = inits,
                 n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iterations = n.update)

# chain (n.iter)
samplesBinLinkBetaPriorComplete = coda.samples(jm, variable.names = c(parsToMonitor), n.iter = n.iterations, thin = n.thin)

samplesBinLinkBetaPriorComplete %<>% recover_types(survDataAnalysisBayes)

saveRDS(samplesBinLinkBetaPriorComplete,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/survivorshipBinLinkBetaPriorComplete.rds")

################################################################################
################################################################################
# JAGS fit for model with binomial likelihood, beta prior, partial pooling
# parameterized via mode
################################################################################
################################################################################

################################################################################
# Initial values
################################################################################

nSites = survDataAnalysisBayesJags$n_site
nYears = survDataAnalysisBayesJags$n_year

kappaInit = 100;

inits = list( list( #theta = matrix(rep(.1,nSites*nYears),nrow=nSites,ncol=nYears),
                          omega=rep(.5,nSites) ,
                          kappaMinusTwo=rep(98,nSites) ),
              list( #theta = matrix(rep(.5,nSites*nYears),nrow=nSites,ncol=nYears),
                    omega=rep(.5,nSites) ,
                    kappaMinusTwo=rep(98,nSites) ),
              list( #theta = matrix(rep(.9,nSites*nYears),nrow=nSites,ncol=nYears),
                    omega=rep(.5,nSites) ,
                    kappaMinusTwo=rep(98,nSites) )
              )

################################################################################
# Set parameters to monitor
################################################################################

parsToMonitor = c("theta","omega","kappa")

################################################################################
# Set JAGS parameters
################################################################################

# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 500
n.update = 5000
n.iterations = 10000
n.thin = 1

################################################################################
# Fit model 
################################################################################

# set random seed
set.seed(2)

# tuning (n.adapt)
jm = jags.model(paste0(dirJagsScripts,"survivorshipModel-binLik-betaPriorMode-partialPooling.R"), 
                 data = survDataAnalysisBayesJags, inits = inits,
                 n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iterations = n.update)

# chain (n.iter)
samplesBinLinkBetaPriorPartialMode = coda.samples(jm, variable.names = c(parsToMonitor), 
                              n.iter = n.iterations, thin = n.thin)

samplesBinLinkBetaPriorPartialMode %<>% recover_types(survDataAnalysisBayes)

saveRDS(samplesBinLinkBetaPriorPartialMode,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/survivorshipBinLinkBetaPriorPartialMode.rds")

################################################################################
################################################################################
# JAGS fit for model with binomial likelihood, beta prior, partial pooling
# parameterized via mean
################################################################################
################################################################################

################################################################################
# Initial values
################################################################################

nSites = survDataAnalysisBayesJags$n_site
nYears = survDataAnalysisBayesJags$n_year

inits = list( list( phi=rep(.1,nSites) ,
                    kappa=rep(1.1,nSites) ),
              list( phi=rep(.5,nSites) ,
                    kappa=rep(1.5,nSites) ),
              list( phi=rep(.9,nSites) ,
                    kappa=rep(2,nSites) )
)

################################################################################
# Set parameters to monitor
################################################################################

parsToMonitor = c("theta","phi","kappa")

################################################################################
# Set JAGS parameters
################################################################################

# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 500
n.update = 5000
n.iterations = 10000
n.thin = 1

################################################################################
# Fit model 
################################################################################

# set random seed
set.seed(2)

# tuning (n.adapt)
jm = jags.model(paste0(dirJagsScripts,"survivorshipModel-binLik-betaPriorMean-partialPooling.R"), 
                 data = survDataAnalysisBayesJags, inits = inits,
                 n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iterations = n.update)

# chain (n.iter)
samplesBinLinkBetaPriorPartialMean = coda.samples(jm, variable.names = c(parsToMonitor), 
                                              n.iter = n.iterations, thin = n.thin)

samplesBinLinkBetaPriorPartialMean %<>% recover_types(survDataAnalysisBayes)

saveRDS(samplesBinLinkBetaPriorPartialMean,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/survivorshipBinLinkBetaPriorPartialMean.rds")

################################################################################
################################################################################
# JAGS fit for model with binomial likelihood, logit parameterization, partial pooling
################################################################################
################################################################################

################################################################################
# Initial values
################################################################################
 
nSites = survDataAnalysisBayesJags$n_site
nYears = survDataAnalysisBayesJags$n_year
 
inits = list(list(mu.alpha = rep(rnorm(1),nSites), 
                  sigma.site = rep(rlnorm(1),nSites)),
             list(mu.alpha = rep(rnorm(1),nSites), 
                  sigma.site = rep(rlnorm(1),nSites)),
             list(mu.alpha = rep(rnorm(1),nSites),
                  sigma.site = rep(rlnorm(1),nSites))
) 


################################################################################
# Set parameters to monitor
################################################################################

parsToMonitor = c("mu.alpha","sigma.site","theta","alpha.i")

################################################################################
# Set JAGS parameters
################################################################################

# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 500
n.update = 5000
n.iterations = 1000
n.thin = 1

################################################################################
# Fit model 
################################################################################

# set random seed
set.seed(2)

# tuning (n.adapt)
jm = jags.model(paste0(dirJagsScripts,"survivorshipModel-binLik-logitLink-partialPooling.R"), 
                 data = survDataAnalysisBayesJags, inits = inits,
                 n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iterations = n.update)

# chain (n.iter)
samplesBinLinkLogitLinkPartial = coda.samples(jm, variable.names = c(parsToMonitor), 
                              n.iter = n.iterations, thin = n.thin)

samplesBinLinkLogitLinkPartial %<>% recover_types(survDataAnalysisBayes)

saveRDS(samplesBinLinkLogitLinkPartial,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/survivorshipBinLinkLogitLink.rds")
