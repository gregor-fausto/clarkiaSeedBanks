# list of JAGS scripts used in this file
# seedBagsCompletePoolingViabilityPartialPoolingBurialJAGS.R
# seedBagsPartialPoolingViabilityPartialPoolingBurialJAGS.R

# -------------------------------------------------------------------
# Models for joint estimates of year 1 belowground rates
# Seed survival, germination, and viability
# Models use log-odds parameterization
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

set.seed(10)
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize seed bag data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# setwd and read data files
load("/Users/Gregor/Dropbox/clarkiaSeedBanks/library/dataFromWorkflowFile/seedBagsData.rda")
df <- seedBags

# -------------------------------------------------------------------
# Clean and organize seed bag data
# -------------------------------------------------------------------
df$seedStart<-100

df$seedStart <- as.integer(df$seedStart)

df$totalJan <- ifelse(df$totalJan>100,100,df$totalJan)

# `MC III 5 13 1 is a problem (91 intact, 17 seedlings)`

# df <- df %>% dplyr::rename(bagNo=bag)

## FOR NOW REMOVE MISSING DATA AND PROBLEMS
# the main issue will be bags that were not recovered in January
# but then recovered in October; there is no way of getting an estimate on how many 
# seeds 'started' those trials
df %>% dplyr::filter(is.na(intactJan<intactOct))

df<-subset(df,!is.na(df$totalJan))
df<-subset(df,!is.na(df$intactOct))
df<-subset(df,!(intactJan<intactOct))

seedBagExperiment <- df
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize viability trial data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# setwd and read data files
load("/Users/Gregor/Dropbox/clarkiaSeedBanks/library/dataFromWorkflowFile/viabilityRawData.rds")
df <- viabilityRawData

df <- df %>% 
  dplyr::select(-c(germPerc,germNot,viabPerc,viabPerc2,condTest))

#df <- df %>% dplyr::rename(bagNo=bag)
df$bag <-as.integer(as.numeric(df$bagNo))

viabilityExperiment <- df

# another check
seedBurial <- seedBagExperiment %>% dplyr::select(site,bagNo,round,age) %>%
  dplyr::mutate(seedBurial = 1)
viabilityTrial <- viabilityExperiment %>% dplyr::select(site,bag,round,age) %>%
  dplyr::mutate(viabilityTrial = 1)

joinedDF<-seedBurial %>% 
  dplyr::full_join(viabilityTrial) 

# dplyr::group_by(site,round,age) %>%
#   dplyr::count(seedBurial,viabilityTrial)

# one row is coded differently so that 
# viabStart=NA and viabStain=NA
# all others have viabStart=NA and viabStain=NA
# recode
viabilityExperiment[is.na(viabilityExperiment$viabStart),]$viabStart = 0

## data check
viabilityExperiment %>% dplyr::filter(germStart - germCount - viabStart<0) 

# filter out rows with problems
# these need to be corrected
viabilityExperiment<-viabilityExperiment %>% 
  dplyr::filter(germStart - germCount - viabStart >= 0)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------
# -------------------------------------------------------------------

## filter the dataset for testing purposes
filterData<-function(x) {
  x %>%
    dplyr::filter(age==1) %>%
    dplyr::filter(site=="BG")
}

seedBagExperiment<-filterData(seedBagExperiment)
viabilityExperiment<-filterData(viabilityExperiment)

# assign variable that combines site and bag; unique id for each bag
# for each dataset, create that unique identifier again and then
# use that to link it to the reference identifier created above
seedBagExperiment<-seedBagExperiment %>%
  tidyr::unite(col='idBag', c(site,bagNo,round,age), sep="", remove=FALSE) %>%
  tidyr::unite(col='idSiteRound', c(site,round), sep="", remove=FALSE)

viabilityExperiment<-viabilityExperiment %>%
  tidyr::unite(col='idBag', c(site,bagNo,round,age), sep="", remove=FALSE) %>%
  tidyr::unite(col='idSiteRound', c(site,round), sep="", remove=FALSE) 

# this line creates a unique id for the subsetted data that is then 
# used to index each of the 2 datasets
# and provides the reference set of bags that were included in the experiment

#https://community.rstudio.com/t/how-to-add-a-counter-to-each-group-in-dplyr/12986/2

referenceTableBag<-data.frame(idBag=union(seedBagExperiment$idBag, viabilityExperiment$idBag)) %>%
  dplyr::mutate(indexBag = 1:length(idBag))

referenceTableSite<-data.frame(site=union(seedBagExperiment$site, viabilityExperiment$site)) %>%
  dplyr::mutate(indexSite = 1:length(site)) 

# referenceTableSiteRound<-data.frame(idSiteRound=union(seedBagExperiment$idSiteRound, viabilityExperiment$idSiteRound)) %>%
#   dplyr::mutate(indexSiteRound = 1:length(idSiteRound)) 

seedBagExperiment<-seedBagExperiment %>%
  dplyr::left_join(referenceTableBag,by="idBag") %>%
  dplyr::left_join(referenceTableSite,by="site")

viabilityExperiment<-viabilityExperiment %>%
  dplyr::left_join(referenceTableBag,by="idBag") %>%
  dplyr::left_join(referenceTableSite,by="site")

# exclude these bags; not in the seed bag dataset
# should check these
d<-referenceTableBag$idBag[!(referenceTableBag$idBag %in% seedBagExperiment$idBag)]
referenceTableBag<-referenceTableBag %>% dplyr::filter(!(idBag%in%d))

viabilityExperiment<-viabilityExperiment %>% dplyr::filter(!(idBag%in%d))


# seedBagExperiment <- seedBagExperiment %>%
#   dplyr::group_by(site,round) %>%
#   dplyr::mutate(indexSiteRound = dplyr::row_number()) %>%
#   dplyr::group_by(site) %>%
#   dplyr::mutate(indexSite = dplyr::row_number()) %>%
#   dplyr::select(idBag,site,idSiteRound,indexBag,indexSite,indexSiteRound) %>% View
# 
# viabilityExperiment <- viabilityExperiment %>%
#   dplyr::group_by(site,yearStart) %>%
#   dplyr::mutate(indexSiteYear = dplyr::row_number()) %>%
#   dplyr::group_by(site) %>%
#   dplyr::mutate(indexSite = dplyr::row_number())

# relevel variable
# not necessary CAN DELETE
# seedBagExperiment$siteBag<-forcats::fct_relevel(seedBagExperiment$siteBag, as.vector(unique(ve$siteBag)))
# viabilityExperiment$siteBag<-forcats::fct_relevel(viabilityExperiment$siteBag, as.vector(unique(ve$siteBag)))

d<-unique(data.frame(site=as.double(as.numeric(as.factor(seedBagExperiment$site))),
                     siteyear=as.double(as.numeric(as.factor(seedBagExperiment$idSiteRound)))))

# pass data to list for JAGS
data = list(
  # nbags comes from a reference table
  # that indexes all the bags across both experiments
  nbags = max(referenceTableBag$indexBag), 
  nsites = max(referenceTableSite$indexSite),
  site_year = as.double(d$site),
  
  # Seed burial experiment, year one
  y_total = as.double(seedBagExperiment$totalJan),
  n_buried = as.double(seedBagExperiment$seedStart),
  
  N_burial = nrow(seedBagExperiment),
  bag_burial = as.double(seedBagExperiment$indexBag),
  
  site = as.double(seedBagExperiment$indexSite),
  
  siteyear = as.double(as.numeric(as.factor(seedBagExperiment$idSiteRound))),
  nsiteyears = length(unique(seedBagExperiment$idSiteRound))
  
)

#save(data,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBagIndexSiteModelData.rds")
#save(seedBagExperiment,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBagExperimentIndexSiteData.rds")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 30000
n.update = 5000
n.iter = 10000
n.thin = 1

set.seed(10)
dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/jagsScripts/")

# # -------------------------------------------------------------------
# # -------------------------------------------------------------------
# New attempt
# # -------------------------------------------------------------------
# # -------------------------------------------------------------------

# set inits for JAGS
inits = list(list(mu.alpha = rep(rnorm(0),data$nsiteyears), sigma.site = rep(rlnorm(1),data$nsites)),
             list(mu.alpha = rep(rnorm(0),data$nsiteyears), sigma.site = rep(rlnorm(1),data$nsites)),
             list(mu.alpha = rep(rnorm(0),data$nsiteyears), sigma.site = rep(rlnorm(1),data$nsites))
)

# Call to JAGS

# tuning (n.adapt)
jm = jags.model(paste0(dir,"testJags.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

parsToMonitor = c("sigma.year","mu.alpha","beta.i")
# chain (n.iter)
zc = coda.samples(jm, variable.names = c(parsToMonitor), n.iter = n.iter, thin = n.thin)

#save(zc,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBagsEvansModelFit.rds")

MCMCsummary(zc,params=c("sigma.year","mu.alpha","beta.i"))
MCMCtrace(zc,params=c("mu.alpha","sigma.year"))

d<-matrix(rep(NA,3*length(MCMCchains(zc,params="mu.alpha"))),ncol=3)
for(i in 1:3){
  d[,i]<-MCMCchains(zc,params="mu.alpha")+MCMCchains(zc,params="beta.i")[,i]
}

out<-data.frame(cbind(
(seedBagExperiment %>% dplyr::group_by(site,round) %>% dplyr::summarise(mean(totalJan/seedStart))),
bayes=apply(boot::inv.logit(MCMCchains(zc,params="mu.alpha")),2,mean)))
out

