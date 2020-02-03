# -------------------------------------------------------------------
# Checking individual models
# to look for missing data values
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
library(bayesplot)
library(stringr)

set.seed(10)
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
setwd("~/Dropbox/projects/clarkiaScripts/data/cleanData")
load(file="sigmaDF.RData") 

setwd("~/Dropbox/projects/clarkiaScripts/data/cleanData")
load(file="fecDF.RData")

setwd("~/Dropbox/projects/clarkiaScripts/data/cleanData")
load(file="phiIndDF.RData")

### MULTIPLE SITE, MULTIPLE YEAR MODEL (POOLED)

siteFilter = as.character(unique(sigmaDF$site))[1:20]
yearFilter = 2007:2012

# need to exlude 2 observations if using all data points without latent states
sigmaDat <- sigmaDF %>% 
  dplyr::filter(year %in% yearFilter & site %in% siteFilter) %>%
  dplyr::select(c(year,site,transect,position,noSeedlings,noFruitingPlants)) %>%
  dplyr::mutate( p = noFruitingPlants/noSeedlings ) %>%
  dplyr::mutate( test = ifelse(p>1, 1, 0)) %>%
  dplyr::filter( test < 1 | is.na(test)) %>% 
  dplyr::filter( !is.na(noSeedlings)) %>%
  dplyr::select(c(site, year, transect, position, noSeedlings, noFruitingPlants)) %>%
  unique

# remove plots without plants (NA)
fecDat <- fecDF %>% 
  dplyr::filter(year %in% yearFilter & site %in% siteFilter) %>%
  dplyr::select(c(year,site,transect,position,totalFruitEquivalents)) %>%
  dplyr::filter( !is.na(totalFruitEquivalents) )

# remove plots without plants (NA)
phiDat <- phiIndDF %>% 
  dplyr::filter(year %in% yearFilter & site %in% siteFilter) %>%
  dplyr::select(c(year,site,response))

# index 
sigmaDat$site.index <- as.integer(as.factor(sigmaDat$site))
sigmaDat$year.index <- as.integer(as.factor(sigmaDat$year))

fecDat$site.index <- as.integer(as.factor(fecDat$site))
fecDat$year.index <- as.integer(as.factor(fecDat$year))

phiDat$site.index <- as.integer(as.factor(phiDat$site))
phiDat$year.index <- as.integer(as.factor(phiDat$year))

s<-sigmaDat %>% 
  tidyr::complete(site.index, nesting(year.index)) %>%
  dplyr::filter(is.na(noFruitingPlants)) %>%
  dplyr::select(site.index,year.index)

f<-fecDat %>% 
  tidyr::complete(site.index, nesting(year.index)) %>%
  dplyr::filter(is.na(totalFruitEquivalents)) %>%
  dplyr::select(site.index,year.index)

p<-phiDat %>% 
  tidyr::complete(site.index, nesting(year.index)) %>%
  dplyr::filter(is.na(response)) %>% 
  dplyr::select(site.index,year.index)

missing_dat<-unique(rbind(s,f,p))

save(missing_dat,file="~/Dropbox/modelsF2019/output/missingness")

