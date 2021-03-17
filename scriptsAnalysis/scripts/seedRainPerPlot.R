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

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Seedling survival to fruiting
# -------------------------------------------------------------------
# -------------------------------------------------------------------
censusSeedlingsFruitingPlants <- readRDS("~/Dropbox/dataLibrary/postProcessingData/censusSeedlingsFruitingPlants.RDS")

noGerms=censusSeedlingsFruitingPlants %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(sdlg = sum(seedlingNumber),n.plot=n()) %>%
  dplyr::filter(sdlg==0)

fruitPl=censusSeedlingsFruitingPlants %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(fruitPl = mean(fruitplNumber))
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Seeds per plant
# -------------------------------------------------------------------
# -------------------------------------------------------------------
seedsPerPlant <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/seedsPerPlantPosterior.RDS")

medianEst <- function(x){return(apply(x,2,median))}
seedsPerPlant.med = medianEst(seedsPerPlant)

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/"
dataFiles <- paste0(directory,list.files(directory))
countSeedPerUndamagedFruit <- readRDS(dataFiles[[grep("countSeedPerUndamagedFruit",dataFiles)]])
# recover site indices and year indices
siteIndex <- data.frame(siteIndex=unique(countSeedPerUndamagedFruit$site3),site=1:20)
yearIndex <- data.frame(yearIndex=unique(countSeedPerUndamagedFruit$year3),
                        year=1:13) 
siteNames = unique(countSeedPerUndamagedFruit$site3)
df=data.frame(expand.grid(site=siteNames,year=2006:2018),fec=seedsPerPlant.med)
df

noGerms$yeartplus1 = noGerms$year
noGerms$year=noGerms$year-1
df$yeart = df$year

df.sum=noGerms %>% dplyr::left_join(df,by=c("year","site")) %>%
  dplyr::left_join(fruitPl,by=c('site','year')) #%>%
#df.sum$fec = round(df.sum$fec)
                             
s0 <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/s0-pop.RDS")
s1 <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/s1-pop.RDS")
g1 <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/g1-pop.RDS")

medianEst <- function(x){return(apply(x,2,median))}

s0.med = medianEst(s0)
g1.med = medianEst(g1)
s1.med = medianEst(s1)

df.sum=df.sum %>%
  dplyr::left_join(data.frame(site=siteNames,s0=s0.med,g1=g1.med,s1=s1.med),by="site")

test = df.sum %>%
  dplyr::mutate(seedRain = ceiling(ceiling(fec)*fruitPl*s0*s1)) 
