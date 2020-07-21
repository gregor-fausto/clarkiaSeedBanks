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
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(HDInterval)
library(bayesplot)

# -------------------------------------------------------------------
# Import data
# -------------------------------------------------------------------
sigmaSummary<-read.csv("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/parameterSummary/sigmaSummary.csv")[,-1]
fecSummary<-read.csv("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/parameterSummary/fruitEquivalentsPerPlantSummary.csv")[,-1]
phiSummary<-read.csv("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/parameterSummary/seedsSummary.csv")[,-1]

# -------------------------------------------------------------------
# Reproductive success: medians
# -------------------------------------------------------------------
nsite <- 20
nyear <- 6

sigmaSummary <- sigmaSummary %>%
  dplyr::select(site,year,med) %>%
  dplyr::rename(med_sigma = med)

fecSummary <- fecSummary %>%
  dplyr::select(site,year,med) %>%
  dplyr::rename(med_fec = med)

phiSummary <- phiSummary %>%
  dplyr::select(site,year,med) %>%
  dplyr::rename(med_phi = med)

abovegroundMedians <- sigmaSummary %>% 
  dplyr::full_join(fecSummary,by=c("site","year")) %>%
  dplyr::full_join(phiSummary,by=c("site","year")) %>%
  dplyr::filter(year<2019)

rsMedianEstimates <- abovegroundMedians %>% 
  dplyr::mutate(rs = med_sigma*med_fec*med_phi) %>%
  dplyr::filter(year<2019)

# save estimates of RS with sites/years where there is missing data excluded
saveRDS(rsMedianEstimates,file="~/Dropbox/clarkiaSeedBanks/products/dataFiles/rsMedianEstimates.RDS")




# -------------------------------------------------------------------
# Reproductive success: full posterior
# -------------------------------------------------------------------

# recover p0_1 using tidybayes
# write for loop filtering to site
# write for loop filtering to year
# for each site-year, multiply all draws/iterations
# save to a 

directory = "/Users/Gregor/Dropbox/dataLibrary/posteriors/"
simFiles <- paste0(directory,list.files(directory))

sigma_mcmcSamples <- readRDS("~/Dropbox/dataLibrary/posteriors/seedSurvivalSamples.RDS")

mcmcSigma <- sigma_mcmcSamples %>%
  tidybayes::spread_draws(p_1[site,year]) 

fitness_mcmcSamples <- readRDS("~/Dropbox/dataLibrary/posteriors/fitnessSamplesTesting3.RDS")

mcmcPhi <-fitness_mcmcSamples %>%
  tidybayes::spread_draws(mu_py_seeds[site,year])

################################################################################
# Create composite
#################################################################################

data <- readRDS("~/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/data.rds")


countFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countFruitsPerPlantAllPlots.RDS")
siteIndex <- data.frame(siteIndex=unique(countFruitsPerPlantAllPlots$site),site=unique(data$site))
yearIndex <- data.frame(yearIndex=unique(countFruitsPerPlantAllPlots$year),
                        year=unique(data$year)) 

tfeDF<-fitness_mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(mu_py_tfe[site,year]) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  #dplyr::select(-c(site,year)) %>%
  # dplyr::rename(site = siteIndex) %>%
  # dplyr::rename(year = yearIndex) %>%
  dplyr::mutate(mu_py_fruits=mu_py_tfe)

countSeedPerDamagedFruit <- readRDS("~/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/countSeedPerDamagedFruit.rds")
siteIndex <- data.frame(siteIndex=unique(countSeedPerDamagedFruit$site4),site=unique(data$site4))
yearIndex <- data.frame(yearIndex=unique(countSeedPerDamagedFruit$year4),
                        year=unique(data$year4)) 

tfeCompDF<- fitness_mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(mu_py_tfe_comp[site,year])  %>%
  dplyr::ungroup() %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  #dplyr::select(-c(site,year)) %>%
  # dplyr::rename(site = siteIndex) %>%
  # dplyr::rename(year = yearIndex) %>%
  dplyr::mutate(mu_py_fruits=mu_py_tfe_comp)
  
countUndamagedDamagedFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/countUndamagedDamagedFruitsPerPlantAllPlots.rds")
siteIndex <- data.frame(siteIndex=unique(countUndamagedDamagedFruitsPerPlantAllPlots$site2),site=unique(data$site2))
yearIndex <- data.frame(yearIndex=unique(countUndamagedDamagedFruitsPerPlantAllPlots$year2),
                        year=unique(data$year2)) 

ufDF<-fitness_mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(mu_py_und[site,year])   %>%
  dplyr::ungroup() %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  #dplyr::select(-c(site,year)) %>%
  # dplyr::rename(site = siteIndex) %>%
  # dplyr::rename(year = yearIndex)  %>%
  dplyr::mutate(mu_py_fruits=mu_py_und)

ufRef<-ufDF %>%
  dplyr::select(siteIndex,yearIndex) %>%
  unique() %>% 
  dplyr::mutate(site.year=paste0(siteIndex,".",yearIndex))
tfeCompRef<-tfeCompDF %>%
  dplyr::select(siteIndex,yearIndex) %>%
  unique() %>% 
  dplyr::mutate(site.year=paste0(siteIndex,".",yearIndex))
tfeRef<-tfeDF %>%
  dplyr::select(siteIndex,yearIndex) %>%
  unique() %>% 
  dplyr::mutate(site.year=paste0(siteIndex,".",yearIndex))

# find index with appropriate samples
findIndex<-function(i=siteNumber,j=yearNumber){
  site<-siteIndex$siteIndex[i]
  year<-(2006:2018)[j]
   paste0(site,".",year)
}

# find corresponding matrix
findMat<-function(i=siteNumber,j=yearNumber){
if (findIndex(i,j) %in% tfeRef$site.year) {
  return(tfeDF)
} else if (findIndex(i,j) %in% tfeCompRef$site.year) { 
  return(tfeCompDF)
} else {
  return(ufDF)
}
}

# get samples corresponding to particular site and years
f<-function(chains=zc,site=i,year=j){
  tmp <- chains %>% dplyr::filter(site==i&year==j)
  return(tmp)
}

rsEstimates<-list()
#iter<-max(dim(mcmcSigma)[1],dim(mcmcFec)[1],dim(mcmcPhi)[1])
iter=1000
rs<-matrix(NA,nrow=iter,ncol=13)

# change if number of years changes
for(i in 1:20){
  for(j in 1:13){
    sigma=f(chains=mcmcSigma,site=i,year=j)$p_1
    fec=f(chains=findMat(i,j),site=i,year=j)$mu_py_fruits
    phi=f(chains=mcmcPhi,site=i,year=j)$mu_py_seeds
    
    # ignore missing data for now?
    sigma2 = if(length(sigma)>0) {sample(sigma,iter)} else {rep(NA,iter)}
    fec2= if(length(fec)>0) {sample(fec,iter)} else {rep(NA,iter)}
    phi2= if(length(phi)>0) {sample(phi,iter)} else {rep(NA,iter)}
    
    rs[,j] <- sigma2*fec2*phi2
  }
  rsEstimates[[i]] <- rs
}

saveRDS(rsEstimates,file="~/Dropbox/clarkiaSeedBanks/products/dataFiles/rsVarFullPosterior.RDS")
