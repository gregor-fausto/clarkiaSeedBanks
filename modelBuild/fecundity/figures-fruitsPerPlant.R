#################################################################################
################################################################################
################################################################################
# Code for figures to compare the following modeling approaches for the seedling survivorship data
#
# Scripts by Gregor Siegmund
# fausto.siegmund@gmail.com
# last updated 04-22-2020
#################################################################################
#################################################################################
#################################################################################
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)

library(MCMCvis)
library(tidybayes)
library(tidyverse)
library(magrittr)

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/fruitsPerPlantAllPlots/"
simFiles <- paste0(directory,list.files(directory))
dirFigures = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/")

mcmcSamples <- readRDS(simFiles[[3]])

################################################################################
# Checks
#################################################################################

parsToMonitor = c("mu0","sigma0","gamma","r")

# MCMCtrace(mcmcSamples,params="mu0")
# MCMCtrace(mcmcSamples,params="sigma0")

library(bayesplot)

################################################################################
# Data
#################################################################################

data <- readRDS(simFiles[[2]])
countFruitsPerPlantTransects <- readRDS(simFiles[[1]])

position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,northing,elevation) %>%
  dplyr::mutate(easting=easting/1000,northing=northing/1000)

################################################################################
# Site-level
#################################################################################

siteNames <- unique(countFruitsPerPlantTransects$site)

mcmcSummary<-mcmcSamples %>%
  
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(mu0[site]) %>%
  dplyr::group_by(site) %>%
  dplyr::mutate(mu0=exp(mu0)) %>%
  dplyr::summarise(med = median(mu0), 
                   ci.lo = quantile(mu0,probs=0.025), 
                   ci.hi = quantile(mu0,probs=0.975),
                   ci.lo2 = quantile(mu0,probs=0.25), 
                   ci.hi2 = quantile(mu0,probs=0.75)
  )
mcmcSummary<-cbind(mcmcSummary[,-1],site=siteNames)

vr = position %>% 
  dplyr::left_join(mcmcSummary,by="site")

pdf(paste0(dirFigures,"spatial-fruitsPerPlant-allplots.pdf"), width=8, height=6)

par(mfrow=c(1,1))
plot(vr$easting,vr$med,ylim=c(0,20),
     ylab="Fruits per plant",
     xlab="Easting (km)", 
     cex.lab = 1.5, cex.axis = 1.5,
     pch=16,type='n')

points(vr$easting,vr$med,ylim=c(0.0,1),
       pch=16)

segments(x0=vr$easting, y0=vr$ci.lo, y1=vr$ci.hi)
segments(x0=vr$easting, y0=vr$ci.lo2, y1=vr$ci.hi2,lwd=2)

dev.off()

# p_1<-MCMCchains(mcmcSamples,params="p_1")
# 
# plot(density(p_1),type='n',ylim=c(0,20))
# for(i in 1:20){
#   lines(density(p_1[,i]),lwd=.5)
# }
# 
# ggplot(data=vr) +
#   geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +
#   geom_linerange(aes(x=easting,ymin=ci.lo2,ymax=ci.hi2),size=1) +
#   geom_point(aes(x=easting,y=med)) +
#   ylim(c(0,1)) +
#   theme_bw()
# 
# reorder(Names, Proportion, mean)
# 
# ggplot(data=vr) +
#   geom_linerange(aes(x=reorder(site,desc(med)),ymin=ci.lo,ymax=ci.hi),size=.25) +
#   geom_linerange(aes(x=reorder(site,desc(med)),ymin=ci.lo2,ymax=ci.hi2),size=1) +
#   geom_point(aes(x=reorder(site,desc(med)),y=med)) +
#   coord_flip() +
#   ylim(c(0,1)) +
#   theme_bw()

################################################################################
# Site-year
#################################################################################

siteIndex <- data.frame(siteIndex=unique(countFruitsPerPlant$site),site=1)
yearIndex <- data.frame(year=1:6,yearIndex=2007:2012)

mcmcSummary<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(gamma[site,year]) %>%
  dplyr::mutate(lambda = exp(gamma)) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(med = median(lambda), 
                   ci.lo = quantile(lambda,probs=0.025), 
                   ci.hi = quantile(lambda,probs=0.975),
                   ci.lo2 = quantile(lambda,probs=0.25), 
                   ci.hi2 = quantile(lambda,probs=0.75)
  )
mcmcSummary<-mcmcSummary %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(site,year)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)

# fruitsPerPlantSummary <- mcmcSummary %>%
#   dplyr::rename(ci.lo95 = ci.lo,
#                 ci.hi95 = ci.hi,
#                 ci.lo50 = ci.lo2,
#                 ci.hi50 = ci.hi2)
# 
# fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/parameterSummary/")
# dir.create(file.path(fileDirectory), showWarnings = FALSE)
# # 
# saveRDS(fruitsPerPlantSummary,file=paste0(fileDirectory,"fruitsPerPlantSummary.RDS"))

vr = position %>% 
  dplyr::left_join(mcmcSummary,by="site")

#pdf(paste0(dirFigures,"spatial-sigma.pdf"), width=8, height=6)

par(mfrow=c(1,1))
plot(vr$easting,vr$med,ylim=c(0.0,1),
     ylab="Seedling survival probability [P(sigma)]",
     xlab="Easting (km)", 
     cex.lab = 1.5, cex.axis = 1.5,
     pch=16,type='n')

points(vr$easting,vr$med,ylim=c(0.0,1),
       pch=16)

segments(x0=vr$easting, y0=vr$ci.lo, y1=vr$ci.hi)
segments(x0=vr$easting, y0=vr$ci.lo2, y1=vr$ci.hi2,lwd=2)

#dev.off()

ggplot(data=vr) +
  geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +
  geom_linerange(aes(x=easting,ymin=ci.lo2,ymax=ci.hi2),size=1) +
  geom_point(aes(x=easting,y=med)) +
  facet_wrap(~year,nrow=2) +
  ylim(c(0,1)) +
  theme_bw()

ggplot(data=vr) +
  geom_linerange(aes(x=year,ymin=ci.lo,ymax=ci.hi),size=.25) +
  geom_linerange(aes(x=year,ymin=ci.lo2,ymax=ci.hi2),size=1) +
  geom_point(aes(x=year,y=med)) +
  facet_wrap(~site,nrow=4) +
  ylim(c(0,1)) +
  theme_bw()

