
# log likelihood
# https://sourceforge.net/p/mcmc-jags/discussion/610036/thread/8211df61/

# 
# https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.43.903&rep=rep1&type=pdf

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
library(magrittr)
library(tidybayes)
library(parallel)
library(stringr)
library(loo)

set.seed(10)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Read data
# -------------------------------------------------------------------
# -------------------------------------------------------------------

data = readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/seedBagMaster.RDS")

data = data %>%
  #dplyr::filter(site=="EC") %>%
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="-", remove=FALSE) %>%
  tidyr::unite(col='siteBag', c(site,bagNo), sep="-", remove=FALSE) %>%
  dplyr::mutate(siteBag = as.factor(siteBag))   %>%
  dplyr::mutate(year = as.factor(yearStart),
                age = as.factor(age)) %>%
  dplyr::select(siteBag,site,year,age,totalJan,seedStart,seedlingJan,intactOct) %>%
  dplyr::rename(siteBags = site,
                yearBags = year,
                ageBags = age)

refIndex <- data.frame(name=rep(c("totalJan","intactOct"),3),
                       ageBags=as.factor(c(1,1,2,2,3,3)),
                       betaIndex=c(1,2,3,4,5,6),
                       months=c(3,12,15,24,27,36),
                       g_1 = c(0,1,1,1,1,1),
                       g_2 = c(0,0,0,1,1,1),
                       g_3 = c(0,0,0,0,0,1)) %>%
  dplyr::mutate(months=months/36)

survivalData=data %>%
  dplyr::select(-seedlingJan) %>%
  tidyr::pivot_longer(cols=c(totalJan,intactOct)) %>%
  dplyr::left_join(refIndex,by=c("name","ageBags")) %>%
  dplyr::rename(response = name, y = value, siteSurvival = siteBags) %>%
  dplyr::mutate(gIndexSurvival = ageBags,
                yearSurvival = yearBags)

#survivalData=survivalData[!duplicated(survivalData$months),]
# survivalData=survivalData[!duplicated(interaction(survivalData$yearSurvival,
#                                                   survivalData$gIndexSurvival)),]

gCompIndex <- data.frame(yearSurvival = as.factor(c(rep(2006,6),rep(2007,4),rep(2008,2))),
                         ageBags=as.factor(c(1,1,2,2,3,3,1,1,2,2,1,1)),
                         compIndex=c(1:12),
                         response=rep(c("totalJan","intactOct"),6)) 

survivalData<-survivalData %>% 
  dplyr::left_join(gCompIndex,by=c('yearSurvival','ageBags',"response"))

#survivalData=survivalData[!duplicated(survivalData$compIndex),]


germinationData=data %>%
  dplyr::select(-c(seedStart,intactOct)) %>%
  dplyr::mutate(gIndex=ageBags) %>%
  dplyr::rename(siteGermination = siteBags,
                yearGermination = yearBags) %>%
  dplyr::mutate(seedlingJan2=seedlingJan)


data <- tidybayes::compose_data(survivalData,germinationData)

data$n1 = dim(survivalData)[1]
data$n2 = dim(germinationData)[1]

refDf<-data.frame(year=data$yearGermination,gIndex=data$gIndex) %>%
  unique %>%
  dplyr::mutate(germinationReference=1:6)

data$yearRefGerm = refDf$year
data$indexRefGerm = refDf$gIndex
data$refGerm = refDf$germinationReference

d=survivalData %>%
  dplyr::mutate(index=1:dim(survivalData)[1]) %>%
  dplyr::group_by(months,yearBags) %>%
  dplyr::top_n(n=1)  %>%
  dplyr::ungroup() %>%
  dplyr::rename(seedStart_pred=seedStart,
                siteSurvival_pred=siteSurvival,
                yearSurvival_pred=yearSurvival,
                months_pred=months,
                compIndex_pred=compIndex) %>%
  dplyr::select(seedStart_pred,siteSurvival_pred,yearSurvival_pred,months_pred,compIndex_pred)

d
d<-tidybayes::compose_data(d)
names(d)[8]="n_pred"
data=c(data,d)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.chain = 3
n.adapt = 3000
n.update = 5000
n.iterations = 10000
n.thin = 1

dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/priorChecking/jagsScripts/")

# set inits for JAGS
inits <- list()

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Weibull residence time, age-dependent germination
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# # tuning (n.adapt)
rjags::load.module("glm", quiet = TRUE)

jmWrAg = jags.model(paste0(dir,"binomialLikelihood-logLink-wrDecay-4.R"), data = data, n.adapt = n.adapt,
                    n.chains=3)

# burn-in (n.update)
update(jmWrAg, n.iter = n.update)

samples.rjags = coda.samples(jmWrAg,
                             variable.names =
                               # GERMINATION
                               c("mu0_g", "mu_g",
                                 "sigma0_g","sigma_g",
                                 #"tau_g",#"tau0_g",
                                # "g",
                                # "mu0_pred","tau0_pred","mu_pred","tau_pred",
                                # "sigma0_pred","sigma_pred","g_pred",
                                "logLik_g",
                                 #PRIOR PREDICTIVE
                                 "y_pred","y_prior.pred",
                                 # SURVIVAL
                                 "mu0_s","mu_s",
                                "sigma0_s","sigma_s",
                                "a", 
                                #"inv.b", "theta_c","mu",
                                # "beta","beta_mod","alpha_s","theta_s",
                                 # POSTERIOR PREDICTIVE
                                # "y_sim","y_surv","mu_survival",
                                 "logLik_y"),
                             n.iter = n.iterations, thin = n.thin)

fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/decayModel/")
dir.create(file.path(fileDirectory), showWarnings = FALSE)
# 
samples.rjags=samples.rjags
saveRDS(samples.rjags,file=paste0(fileDirectory,"seedBurialSamples.rds"))
saveRDS(data,file=paste0(fileDirectory,"data.rds"))
#saveRDS(seedBagsData,file=paste0(fileDirectory,"seedBagExperiment.rds"))
#saveRDS(viabilityRawData,file=paste0(fileDirectory,"viabilityExperiment.rds"))


