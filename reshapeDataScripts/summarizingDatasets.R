# -------------------------------------------------------------------
# script to tally number of observations for each parameter
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
library(xtable)

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

df<-subset(df,!is.na(df$totalJan))
df<-subset(df,!is.na(df$intactOct))
df<-subset(df,!(intactJan<intactOct))

seedBagExperiment <- df

# seed burial summary
seedBagExperimentSummary = seedBagExperiment %>%
#  dplyr::select(site,yearData,age,totalJan,seedStart) %>%
 # dplyr::filter(age==1) %>%
  dplyr::group_by(site,yearData,age) %>%
  dplyr::summarise(count=n()) %>%
  tidyr::pivot_wider(names_from = c(yearData,age), values_from = count) %>%
  dplyr::select(site, `2007_1`,`2008_1`,`2009_1`,`2008_2`,`2009_2`,`2009_3`)

print(xtable(seedBagExperimentSummary, type = "latex", align="lccccccc"),
      include.rownames=FALSE)

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

# seed burial summary
viabilityExperimentSummary = viabilityExperiment %>%
  dplyr::select(site,round,age,bag) %>% 
  unique %>%
    dplyr::group_by(site,round,age) %>% 
  dplyr::summarise(count=n()) %>%
  tidyr::pivot_wider(names_from = c(round,age), values_from = count) %>%
  dplyr::select(site, `1_1`,`2_1`,`3_1`,`1_2`,`2_2`,`1_3`)

print(xtable(viabilityExperimentSummary, type = "latex", align="lccccccc"),
      include.rownames=FALSE)


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize seedling survival data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
setwd("~/Dropbox/projects/clarkiaScripts/data/cleanData")
load(file="sigmaDF.RData") 

### MULTIPLE SITE, MULTIPLE YEAR MODEL (POOLED)

siteFilter = as.character(unique(sigmaDF$site))[1:20]
yearFilter = 2006:2015

# need to exclude 2 observations if using all data points without latent states
sigmaSummary <- sigmaDF %>% 
  dplyr::filter(year %in% yearFilter & site %in% siteFilter) %>%
  dplyr::select(c(year,site,transect,position,noSeedlings,noFruitingPlants)) %>% 
  dplyr::mutate( p = noFruitingPlants/noSeedlings ) %>%
  dplyr::mutate( test = ifelse(p>1, 1, 0)) %>% 
  dplyr::filter( !is.na(noSeedlings)) %>%
  dplyr::select(site,year,transect,position) %>% 
  unique %>%
  dplyr::group_by(site,year) %>% 
  dplyr::summarise(count=n()) %>%
  tidyr::pivot_wider(names_from = c(year), values_from = count)

print(xtable(sigmaSummary, type = "latex", align="lccccccccccc"),
      include.rownames=FALSE,NA.string="NA")  

# need to exclude 2 observations if using all data points without latent states
sigmaUndercountSummary <- sigmaDF %>% 
  dplyr::filter(year %in% yearFilter & site %in% siteFilter) %>%
  dplyr::select(c(year,site,transect,position,noSeedlings,noFruitingPlants)) %>% 
  dplyr::mutate( p = noFruitingPlants/noSeedlings ) %>%
  dplyr::mutate( test = ifelse(p>1, 1, 0)) %>% 
  dplyr::filter( !is.na(noSeedlings)) %>%
  dplyr::mutate( n = 1) %>% 
  dplyr::group_by(site,year) %>% 
    dplyr::summarise(sum1= sum(n[p<=1],na.rm=TRUE),
                     sum2 = sum(n[p>1], na.rm = TRUE)) %>% 
  dplyr::mutate(prop = sum2/(sum1+sum2)) %>% 
  dplyr::select(site,year,prop) %>%
  tidyr::pivot_wider(names_from = c(year), values_from = prop)

# manually edit table to set 0.00 -> 0
print(xtable(sigmaUndercountSummary, type = "latex", align="lccccccccccc",digits=2),
      include.rownames=FALSE,NA.string="NA")  
  
plot.sum =  sigmaDF %>% 
  dplyr::filter(year %in% yearFilter & site %in% siteFilter) %>%
  dplyr::select(c(year,site,transect,position,noSeedlings,noFruitingPlants)) %>% 
  dplyr::mutate( p = noFruitingPlants/noSeedlings ) %>%
  dplyr::mutate( test = ifelse(p>1, 1, 0)) %>% 
  dplyr::filter( !is.na(noSeedlings)) %>%
  dplyr::mutate( n = 1) %>% 
  dplyr::group_by(site,year) %>% 
  dplyr::summarise(sum1= sum(n[p<=1],na.rm=TRUE),
                   sum2 = sum(n[p>1], na.rm = TRUE)) %>% 
  dplyr::mutate(prop = sum2/(sum1+sum2), tot = sum1+sum2) 

ggplot(plot.sum) +
  geom_point(aes(x=year,y=tot)) +
  geom_point(aes(x=year,y=sum2),col="red") +
  facet_wrap(~site) + theme_bw()

p <- ggplot(plot.sum) +
  geom_line(aes(x=year,y=tot)) +
  geom_line(aes(x=year,y=sum2),col="red") +
  facet_wrap(~site) + theme_bw() + 
  scale_x_continuous(breaks = c(2006,2010, 2014))+
  xlab("Year") + ylab("Number of counts")

ggplot(plot.sum) +
  geom_line(aes(x=year,y=prop)) +
  facet_wrap(~site) + theme_bw()

dirFigures = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/")
ggsave(filename=paste0(dirFigures,"underCounting.pdf"), plot=p)


# mdat<-as.matrix(ifelse(sigmaUndercountSummary>0,2,0))
# mdat <- matrix(c(rep(0,3), 0, 5, 0, rep(0, (7*3)), 0, 5, 0), 
#                nrow = 10, ncol=3, byrow=TRUE)
# xtable(df,digits=mdat)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize fecundity data
# -------------------------------------------------------------------
# -------------------------------------------------------------------

setwd("~/Dropbox/projects/clarkiaScripts/data/cleanData")
load(file="fecDF.RData")

tail(fecDF)

siteFilter = as.character(unique(sigmaDF$site))[1:20]
yearFilter = 2006:2012

# remove plots without plants (NA)
fecDF %>% 
  dplyr::filter(year %in% yearFilter & site %in% siteFilter) %>%
  dplyr::select(c(year,site,transect,position,totalFruitEquivalents)) %>% 
  dplyr::filter( !is.na(totalFruitEquivalents) ) %>% 
  dplyr::select(site,year,transect,position) %>% 
  dplyr::group_by(site,year) %>% 
  dplyr::summarise(count=n()) %>%
  tidyr::pivot_wider(names_from = c(year), values_from = count)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import seeds per fruit data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
setwd("~/Dropbox/projects/clarkiaScripts/data/cleanData")
load(file="phiIndDF.RData")

tail(phiIndDF)

siteFilter = as.character(unique(sigmaDF$site))[1:20]
yearFilter = 2006:2015
# remove plots without plants (NA)
phiSummary <- phiIndDF %>% 
  dplyr::filter(variable=="noSeedsPerFruit"|variable=="noSeedsPerUndFruit") %>% 
  dplyr::filter(year %in% yearFilter & site %in% siteFilter) %>%
  dplyr::select(c(year,site,variable,response)) %>%
   dplyr::group_by(site,year) %>% 
   dplyr::summarise(count=n()) %>%
   tidyr::pivot_wider(names_from = c(year), values_from = count)
 
print(xtable(phiSummary, type = "latex", align="lccccccccccc"),
      include.rownames=FALSE,NA.string="NA")  

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize seedling survival data
# -------------------------------------------------------------------
# -------------------------------------------------------------------



setwd("~/Dropbox/projects/clarkiaScripts/data/cleanData")
load(file="phiIndDF.RData")
