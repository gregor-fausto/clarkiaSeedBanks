# list of JAGS scripts used in this file
# seedBagsCompletePoolingViabilityPartialPoolingBurialJAGS.R
# seedBagsPartialPoolingViabilityPartialPoolingBurialJAGS.R

# -------------------------------------------------------------------
# Models for joint estimates of year 1 below ground rates
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
library(magrittr)
library(tidybayes)
library(parallel)
library(stringr)

set.seed(10)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Read data
# -------------------------------------------------------------------
# -------------------------------------------------------------------

data = readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/seedBagMaster.RDS")

data = data %>%
  dplyr::filter(site=="BR") %>%
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="-", remove=FALSE) %>%
  tidyr::unite(col='siteBag', c(site,bagNo), sep="-", remove=FALSE) %>%
  dplyr::mutate(siteBag = as.factor(siteBag))   %>%
  dplyr::mutate(year = as.factor(yearStart),
                age = as.factor(age)) %>%
  dplyr::select(siteBag,site,year,age,totalJan,seedStart,seedlingJan,intactOct) %>%
  dplyr::rename(siteBags = site,
                yearBags = year,
                ageBags = age)

# refMat <- data.frame(name=rep(c("seedStart","seedlingJan","totalJan","intactOct"),3),
#                      months=c(0,3,3,12,0,15,15,24,27,27,27,36),
#                      yearBags=as.factor(c(rep(2006,4),rep(2007,4),rep(2008,4))))

octData = data %>%
  dplyr::select(siteBags,yearBags,ageBags,intactOct) %>%
  dplyr::mutate(
    months = ifelse(ageBags==1 , 12, NA),
    #
    months = ifelse(ageBags==2 , 24, months),
    #
    months = ifelse(ageBags==3 , 36, months)
  ) %>%
  dplyr::mutate(seedStart=100)

data<-data %>%
  dplyr::select(siteBags,yearBags,ageBags,totalJan,seedlingJan) %>%
  dplyr::mutate(
                months = ifelse(ageBags==1 , 3, NA),
                #
                months = ifelse(ageBags==2 , 15, months),
                #
                months = ifelse(ageBags==3 , 27, months)
  ) %>%
  dplyr::mutate(seedStart=100)

data1 <- data %>%
  dplyr::filter(ageBags=="1")

data2 <- data %>%
  dplyr::filter(ageBags=="2")

names(data2) = paste(names(data2),"3",sep="_")

octData1 <- octData %>%
  dplyr::filter(ageBags=="1")

names(octData1) = paste(names(octData1),"2",sep="_")

octData2 <- octData %>%
  dplyr::filter(ageBags=="2")

names(octData2) = paste(names(octData2),"4",sep="_")

data <- tidybayes::compose_data(data1,octData1,data2,octData2)

data$n1 = dim(data1)[1]
data$n2 = dim(data2)[1]

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

dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/seedModelComparison/")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Complete pooling of germination and viability trials
# Partial pooling of seed burial experiment (site level)
# -------------------------------------------------------------------
# -------------------------------------------------------------------


# set inits for JAGS
inits <- list()

# # tuning (n.adapt)
jm = jags.model(paste0(dir,"jagsCompexpConstant.R"), data = data, n.adapt = n.adapt)
              #  inits = inits,
              #  n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

# chain (n.iter)
# samples.rjags = coda.samples(jm,
#                              variable.names = c("beta","k"),
#                              n.iter = n.iterations, thin = n.thin)
# 
MCMCsummary(samples.rjags)

library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(4, "PuOr"))(4)

plot(rep(3,length=data$n1),data$totalJan/data$seedStart,xlim=c(0,24),ylim=c(0,1))
lines(seq(0,24,by=0.1),exp(-mean(MCMCchains(samples.rjags,"k"))*seq(0,24,by=0.1)))
abline(h=mean(data$totalJan/data$seedStart),lty=2,col=colors[1])

points(rep(12,length=data$n1),data$intactOct_2/data$seedStart_2,xlim=c(0,12),ylim=c(0,1))
lines(seq(0,24,by=0.1),exp(-mean(MCMCchains(samples.rjags,"k"))*seq(0,24,by=0.1)))
abline(h=mean(data$intactOct_2/data$seedStart_2),lty=2,col=colors[2])

points(rep(15,length=data$n2),data$totalJan_3/data$seedStart_3)
lines(seq(0,24,by=0.1),exp(-mean(MCMCchains(samples.rjags,"k"))*seq(0,24,by=0.1)))
abline(h=mean(data$totalJan_3/data$seedStart_3),lty=2,col=colors[3])

points(rep(24,length=data$n2),data$intactOct_4/data$seedStart_4)
lines(seq(0,24,by=0.1),exp(-mean(MCMCchains(samples.rjags,"k"))*seq(0,24,by=0.1)))
abline(h=mean(data$intactOct_4/data$seedStart_4),lty=2,col=colors[4])
# fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/posteriors/modelChecking/")
# dir.create(file.path(fileDirectory), showWarnings = FALSE)
#
# saveRDS(samples.rjags,file=paste0(fileDirectory,"seedSamplesAllAges.RDS"))
#

exp(-mean(MCMCchains(samples.rjags,"k"))*seq(0,24,by=0.1))

a=exp(-mean(MCMCchains(samples.rjags,"k"))*3)
b=exp(-mean(MCMCchains(samples.rjags,"k"))*12)
c=exp(-mean(MCMCchains(samples.rjags,"k"))*15)
d=exp(-mean(MCMCchains(samples.rjags,"k"))*24)

a/1
b/a
c/b
d/c


