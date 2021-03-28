# -------------------------------------------------------------------
# Models for germination in seed bag burial experiment
# -------------------------------------------------------------------
# log likelihood
# https://sourceforge.net/p/mcmc-jags/discussion/610036/thread/8211df61/
# https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.43.903&rep=rep1&type=pdf
# -------------------------------------------------------------------
 rm(list=ls(all=TRUE)) # clear R environment
# rm(list=setdiff(ls(all=TRUE),c("dataDirectory","modelDirectory","fileDirectory","n.adapt","n.update","n.iterations","n.thin"))) # if using in source(script)
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
 dataDirectory = "~/Dropbox/dataLibrary/postProcessingData/"
 modelDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/"

# -------------------------------------------------------------------
# Import seed bag burial data
# -------------------------------------------------------------------
seedBagsData = readRDS(paste0(dataDirectory,"seedBagsData.rds"))
 
# -------------------------------------------------------------------
# Manually relabel a bag number and filter out 1 row of data
# -------------------------------------------------------------------
# relabel bag
seedBagsData[seedBagsData$site=="DLW" & seedBagsData$bagNo==43 & seedBagsData$age==1,]$bagNo = 42

# remove row of data with 108 seeds in January
seedBagsData <- seedBagsData %>% dplyr::filter(totalJan<=100)

# -------------------------------------------------------------------
# Prep data for JAGS
# -------------------------------------------------------------------
# Create a variable for the number of seeds at the start of the trial
seedBagsData$seedStart<-as.double(100)

# Rename data frame
data = seedBagsData

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

# create a reference data frame for event histories
gCompIndex <- data.frame(yearSurvival = as.factor(c(rep(2006,6),rep(2007,4),rep(2008,2))),
                         ageBags=as.factor(c(1,1,2,2,3,3,1,1,2,2,1,1)),
                         compIndex=c(1:12),
                         response=rep(c("totalJan","intactOct"),6)) 

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
  dplyr::left_join(df.index,by=c("ageBags","yearGermination")) %>%
  dplyr::mutate(yearGermination=as.factor(yearGermination),ageBags=as.factor(ageBags))

# remove data with missing totalJan from germination dataset
germinationData <- germinationData %>%
  dplyr::filter(!is.na(totalJan))

germinationData = germinationData %>% dplyr::filter(ageBags==1)

# -------------------------------------------------------------------
# Read in climate data
# -------------------------------------------------------------------
climate=readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/climateData.RDS")
climate=climate %>% dplyr::ungroup()
climate <- climate %>% dplyr::filter(intenseDemography==1)
climate <- climate %>% dplyr::filter(year %in% 2005:2007&season=="winter") %>% dplyr::mutate(year = year+1)

climate = germinationData %>%
  dplyr::rename(site = siteGermination , year = yearGermination) %>%
  dplyr::mutate(site=as.character(site),year=as.numeric(as.character(year))) %>%
  dplyr::left_join(climate,by=c("site","year"))

muTemp = mean(climate$t,na.rm=TRUE)
sdTemp = sd(climate$t,na.rm=TRUE)
muPrecip = mean(climate$p,na.rm=TRUE)
sdPrecip = sd(climate$p,na.rm=TRUE)

climate$t.std = (climate$t - muTemp)/sdTemp
climate$p.std = (climate$p - muPrecip)/sdPrecip

# -------------------------------------------------------------------
# Predictions
# -------------------------------------------------------------------
# for plotting probability of germination as function of temperature
temp_pred = seq(0, 15, .25) 
#Standardized for making predictions of probability of occupancy at new
temp_pred_std = (temp_pred - muTemp) / sdTemp

# for plotting probability of germination as function of temperature
precip_pred = seq(0, 200, 5) 
#Standardized for making predictions of probability of occupancy at new
precip_pred_std = (precip_pred - muPrecip) / sdPrecip

# -------------------------------------------------------------------
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------

data <- tidybayes::compose_data(germinationData)

data$n = dim(germinationData)[1]

refDf<-data.frame(year=data$yearGermination,gIndex=data$gIndex) %>%
  unique %>%
  dplyr::mutate(germinationReference=1:6)

data$yearRefGerm = refDf$year
data$indexRefGerm = refDf$gIndex
data$refGerm = refDf$germinationReference
data$t.std = climate$t.std
data$p.std = climate$p.std
data$temp_pred_std = temp_pred_std
data$precip_pred_std = precip_pred_std

detach("package:tidyverse", unload=TRUE)

# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
# -------------------------------------------------------------------
 
n.adapt = 3000
n.update = 5000
n.iterations = 5000
n.thin = 1

# -------------------------------------------------------------------
# Describe functions
# -------------------------------------------------------------------

initsMu0 <- function(samps = 1){
  rnorm(n = samps, mean = 0, sd = 1)
}

initsSigma0Weak <- function(samps = 1){
  extraDistr::rhnorm(n = samps, sigma = 1)
}

initsSigmaWeakMat <- function(samps = data$n_siteGermination){
  extraDistr::rhnorm(n = samps, sigma = 1)
}

initsBeta <- function(samps = 3){
  rnorm(n = samps, mean = 0, sd = 1)
}

# -------------------------------------------------------------------
# Set inits for JAGS
# -------------------------------------------------------------------
inits <- list()
for(i in 1:3){
  inits[[i]] <- list(initsMu0(), initsSigma0Weak(), initsSigmaWeakMat(),
                     initsBeta(), initsMu0(samps=3))

  names(inits[[i]]) = c("mu0_g", "sigma0_g", "sigma_g", "beta", "mu0_y")

}

# -------------------------------------------------------------------
# Call to JAGS
# -------------------------------------------------------------------
jm = jags.model("~/Dropbox/clarkiaSeedBanks/modelBuild/jags-germinationInteractions.R",
                data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

update(jm, n.iter = n.update)

parsToMonitor = c("mu0_g", "sigma0_g", "sigma_g","beta","mu_g",
                  "psiPredict", "phiPredict","mu0_y")

zm = coda.samples(jm, variable.names = c(parsToMonitor),
                  n.iter = n.iterations, thin = 1)

# -------------------------------------------------------------------
# Check convergence
# -------------------------------------------------------------------
MCMCvis::MCMCsummary(zm,c("mu0_g","sigma0_g"))
MCMCvis::MCMCsummary(zm,c("sigma_g"))
MCMCvis::MCMCsummary(zm,c("mu_g"))
MCMCvis::MCMCsummary(zm,c("beta"))


# -------------------------------------------------------------------
# Evaluate interactions
# -------------------------------------------------------------------
psi = MCMCvis::MCMCchains(zm,params="psiPredict")

par(mfrow=c(1,1),mar = c(2,2,2,2))
plot(NA,NA,type='l',xlim=c(0,15),ylim=c(0,1))

for(i in 1:20){
  index=grep(paste0("\\[",i,","),colnames(psi))
  med.pred = apply(psi[,index],2,median)
  temp_pred_std_natural=temp_pred_std*sdTemp+muTemp
  
  lines(temp_pred_std_natural,med.pred,type='l')
}

phi = MCMCvis::MCMCchains(zm,params="phiPredict")

par(mfrow=c(1,1),mar = c(2,2,2,2))
plot(NA,NA,type='l',xlim=c(0,200),ylim=c(0,1))

for(i in 1:20){
  index=grep(paste0("\\[",i,","),colnames(phi))
  med.pred = apply(phi[,index],2,median)
  precip_pred_std_natural=precip_pred_std*sdPrecip+muPrecip
  lines(precip_pred_std_natural,med.pred,type='l')
}

mu_g = MCMCchains(zm,"mu_g")
sigma_g = MCMCchains(zm,"sigma_g")

i=1
index=grep(paste0("\\[",i,"\\]"),colnames(mu_g))
parm1=mu_g[,index]
index=grep(paste0("\\[",i,"\\]"),colnames(sigma_g))
parm2=sigma_g[,index]
alpha=rnorm(dim(mu_g)[1],parm1,parm2)
index=grep(paste0("\\[",i,","),colnames(phi))
alpha+

for(i in 1:n_siteGermination){
  
  alpha_pred[i] ~ dnorm(mu_g[siteGermination[i]],tau_g[siteGermination[i]])
  
  for (j in 1:length(temp_pred_std)) {
    # gets the germination probability across temps at mean precipitation (0, centered) 
    logit(psiPredict[i,j]) <- alpha_pred[i] + beta[1] * temp_pred_std[j]
  }
  
  for (j in 1:length(precip_pred_std)) {
    # gets the germination probability across temps at mean precipitation (0, centered) 
    logit(phiPredict[i,j]) <- alpha_pred[i] + beta[2] * precip_pred_std[j]
  }
  
}


# HPDI <-  MCMCpstr(zm, params=c("psiPredict"), func = function(x) hdi(x, .95))
# el_maxHPDI =  MCMCpstr(zm, params=c("elevationMax"), func = function(x) hdi(x,.95)) 
# 
# par(mfrow = c(1, 2))
# plot(elevation_pred, med_psi$psiPredict, type = "l", ylim = c(0,1), xlab = "Elevation (m)", 
#      lwd = 2, ylab ="Probability of occupancy")
# lines(elevation_pred, HPDI$psiPredict[,1], type = "l", lty = "dashed", lwd = 2)
# lines(elevation_pred, HPDI$psiPredict[,2], type = "l", lty = "dashed", lwd = 2)
# abline(v = mean(MCMCchains(zm,"elevationMax")), lwd = 2)
# text(1300, .9, "Optimum elevation", cex = .8)
