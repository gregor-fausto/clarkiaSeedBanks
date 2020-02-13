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

df <- df %>% dplyr::rename(bagNo=bag)

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

df <- df %>% dplyr::rename(bagNo=bag)
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
# Multiple-site model (3 year)
# Log-odds parameterization for viability, trials summed across bags
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
  yg = as.double(viabilityExperiment$germCount),
  ng = as.double(viabilityExperiment$germStart),
  yv = as.double(viabilityExperiment$viabStain),
  nv = as.double(viabilityExperiment$viabStart),
  N = nrow(viabilityExperiment),
  bag = as.double(viabilityExperiment$siteBag),
  nbags = length(unique(viabilityExperiment$siteBag))
)

# pass data to list for JAGS
# data = list()
#    yg = as.double(dat$seedlingJan),
#   yt = as.double(dat$totalJan),
#   yo = as.double(dat$intactOct),
#   yv = as.double(dat$y_new),
#   n = as.double(dat$seedStart),
#   nv = as.double(dat$n_new),
#   N = nrow(dat),
#   site = as.double(dat$site.index),
#   nsites = length(unique(dat$site.index)),
#   year = as.double(dat$year.index),
#   nyears = length(unique(dat$year.index)),
#   
#   yt2 = as.double(datTwo$totalJan),
#   n2 = as.double(datTwo$seedStart),
#   N2 = nrow(datTwo),
#   site2 = as.double(datTwo$site.index),
#   nsites2 = length(unique(datTwo$site.index)),
#   year2 = as.double(datTwo$year.index),
#   nyears2 = length(unique(datTwo$year.index))
# )

# right now each bag has its own prior;
# want to give nsite priors 
# inits = list(
#   list( sigma.b = matrix(rep(50,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears), 
#         mu.b = matrix(rep(0,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears), 
#         #alpha = rep(-10, data$N),
#         # sigmaS1 = rep(.1, data$nsites), 
#         # sigmaG1 = rep(.1, data$nsites), 
#         # sigmaS2 = rep(.1, data$nsites),
#         alphaS1 = rep(-10, data$nsites), 
#         alphaG1 = rep(-10, data$nsites), 
#         alphaS2 = rep(-10, data$nsites),
#         alphaS3 = rep(0, data$nsites)
#         # betaS1 = matrix(-10, nrow=data$nsites,ncol=data$nyears), 
#         # betaG1 = matrix(-10, nrow=data$nsites,ncol=data$nyears), 
#         # betaS2 = matrix(-10, nrow=data$nsites,ncol=data$nyears)
#   ),
#   list( sigma.b = matrix(rep(10,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears), 
#         mu.b = matrix(rep(0,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears),
#         #alpha = rep(0, data$N),
#         # sigmaS1 = rep(.5, data$nsites), 
#         # sigmaG1 = rep(.5, data$nsites), 
#         # sigmaS2 = rep(.5, data$nsites),
#         alphaS1 = rep(0, data$nsites), 
#         alphaG1 = rep(0, data$nsites), 
#         alphaS2 = rep(0, data$nsites),
#         alphaS3 = rep(0, data$nsites)
#         # betaS1 = matrix(0, nrow=data$nsites,ncol=data$nyears), 
#         # betaG1 = matrix(0, nrow=data$nsites,ncol=data$nyears), 
#         # betaS2 = matrix(0, nrow=data$nsites,ncol=data$nyears)
#   ),
#   list( sigma.b = matrix(rep(20,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears), 
#         mu.b = matrix(rep(0,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears), 
#         #alpha = rep(10, data$N),
#         #sigmaS1 = rep(.9, data$nsites), 
#         #sigmaG1 = rep(.9, data$nsites), 
#         #sigmaS2 = rep(.9, data$nsites),
#         alphaS1 = rep(10, data$nsites), 
#         alphaG1 = rep(10, data$nsites), 
#         alphaS2 = rep(10, data$nsites),
#         alphaS3 = rep(0, data$nsites)
#         # betaS1 = matrix(10, nrow=data$nsites,ncol=data$nyears), 
#         # betaG1 = matrix(10, nrow=data$nsites,ncol=data$nyears), 
#         # betaS2 = matrix(10, nrow=data$nsites,ncol=data$nyears)
#   )
# )


# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 1000
n.update = 10000
n.iter = 10000

# Call to JAGS

dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/jagsScripts/")

# tuning (n.adapt)
jm = jags.model(paste0(dir,"seed_modelJAGS.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

intercepts = c("alphaS1","alphaG1","alphaS2","alphaS3")
#slopes = c("betaS1","betaG1","betaS2")
#variances = c("sigmaS1","sigmaG1","sigmaS2")
viab = c("mu.b","sigma.b")
#sims = c("yv.sim","yt.sim","yg.sim","yo.sim","yt2.sim")



# chain (n.iter)
zc = coda.samples(jm, variable.names = c(intercepts,viab), n.iter = n.iter, thin=10)

save(zc,file="/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/output/seedbagfit.rds")
save(data,file="/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/output/seedbagdata.rds")
