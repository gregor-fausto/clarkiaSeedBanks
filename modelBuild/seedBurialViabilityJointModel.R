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


l1<-seedBagExperiment %>%
  dplyr::select(site,bagNo,round,age) %>%
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="", remove=FALSE) %>%
  unique()

l2<-viabilityExperiment %>%
  dplyr::select(site,bagNo,round,age) %>%
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="", remove=FALSE) %>%
  unique()

referenceTable<-data.frame(id=union(l2$id,l1$id)) %>%
  dplyr::mutate(idNo = 1:length(id)) 

# problem: need a reference list of bags and sites that joins the
# viability dataset and the seed bag experiment

filterData<-function(x) {
  x %>%
    dplyr::filter(age==1) %>%
    dplyr::filter(site=="BG")
}

seedBagExperiment<-filterData(seedBagExperiment)
viabilityExperiment<-filterData(viabilityExperiment)

# assign variable that combines site and bag; unique id for each bag
seedBagExperiment<-seedBagExperiment %>%
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="", remove=FALSE) %>%
  tidyr::unite(col='siteBag', c(site,bagNo), sep="", remove=FALSE) %>%
  dplyr::mutate(siteBag = as.factor(siteBag)) %>%
  dplyr::left_join(referenceTable,by="id")

viabilityExperiment<-viabilityExperiment %>%
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="", remove=FALSE) %>%
  tidyr::unite(col='siteBag', c(site,bagNo), sep="", remove=FALSE) %>%
  dplyr::mutate(siteBag = as.factor(siteBag)) %>%
  dplyr::left_join(referenceTable,by="id")


# relevel variable
# not necessary CAN DELETE
# seedBagExperiment$siteBag<-forcats::fct_relevel(seedBagExperiment$siteBag, as.vector(unique(ve$siteBag)))
# viabilityExperiment$siteBag<-forcats::fct_relevel(viabilityExperiment$siteBag, as.vector(unique(ve$siteBag)))

# pass data to list for JAGS
data = list(
  # Germination and Viability Trials
  yg = as.double(viabilityExperiment$germCount),
  ng = as.double(viabilityExperiment$germStart),
  yv = as.double(viabilityExperiment$viabStain),
  nv = as.double(viabilityExperiment$viabStart),
  N = nrow(viabilityExperiment),
  bag = as.double(viabilityExperiment$idNo),
  nbags = length(unique(viabilityExperiment$idNo)), 
  
  # Seed burial experiment, year one
  y_seedlings = as.double(seedBagExperiment$seedlingJan),
  y_total = as.double(seedBagExperiment$totalJan),
  y_october = as.double(seedBagExperiment$intactOct),
  n_buried = as.double(seedBagExperiment$seedStart),
  N_burial = nrow(seedBagExperiment),
  
  bag_burial = as.double(seedBagExperiment$idNo),
  nbags_burial = length(unique(seedBagExperiment$idNo))
)

save(data,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBagsModelData.rds")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 1000
n.update = 10000
n.iter = 10000
n.thin = 10

set.seed(10)
dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/jagsScripts/")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Complete pooling
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# set inits for JAGS
inits = list(list(pv = .1,pg = .1, pi = .1, ps = .1), 
             list(pv = .5,pg = .5, pi = .5, ps = .5), 
             list(pv = .9,pg = .9, pi = .9, ps = .9))

# Call to JAGS

# tuning (n.adapt)
jm = jags.model(paste0(dir,"seedBagsCompletePoolingJAGS.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

parsToMonitor = c("pv","pg","pi","ps","viability")
sims = c("ygSim","yvSim","ySeedlingsSim","yTotalSim")
# chain (n.iter)
zc_pool = coda.samples(jm, variable.names = c(parsToMonitor,sims), n.iter = n.iter, thin = n.thin)

MCMCsummary(zc_pool, params = c("pv","pg","pi","ps","viability"))
save(zc_pool,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBagsCompletePoolingFit.rds")
MCMCsummary(zc_pool, params = c("ygSim","yvSim","ySeedlingsSim","yTotalSim"))

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# No pooling for both datasets - not possible for joining the datasets;
# Viability has to at least be pooled to the bag level;
# No pooling for the seed burials
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# set inits for JAGS
inits = list(list(pv = rep(.1,data$nbags_burial),pg = rep(.1,data$nbags_burial),
                  pi = rep(.1,data$N_burial), ps = rep(.1,data$N_burial)), 
             list(pv = rep(.5,data$nbags_burial),pg = rep(.5,data$nbags_burial),
                  pi = rep(.5,data$N_burial), ps = rep(.5,data$N_burial)), 
             list(pv = rep(.9,data$nbags_burial),pg = rep(.9,data$nbags_burial),
                  pi = rep(.9,data$N_burial), ps = rep(.9,data$N_burial)))

# Call to JAGS

# tuning (n.adapt)
jm = jags.model(paste0(dir,"seedBagsNoPoolingJAGS.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

parsToMonitor = c("pv","pg","pi","ps","viability")
sims = c("ygSim","yvSim","ySeedlingsSim","yTotalSim")
# chain (n.iter)
zc_nopool = coda.samples(jm, variable.names = c(parsToMonitor,sims), n.iter = n.iter, thin = n.thin)

MCMCsummary(zc_nopool, params = c("pv","pg","pi","ps","viability"))
save(zc_nopool,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBagsNoPoolingFit.rds")
MCMCsummary(zc_nopool, params = c("ygSim","yvSim","ySeedlingsSim","yTotalSim"))

MCMCsummary(zc_nopool, params = c("viability"))

