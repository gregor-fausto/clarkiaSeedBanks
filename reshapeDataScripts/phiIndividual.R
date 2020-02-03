# -------------------------------------------------------------------
# SEEDS PER FRUIT (PHI)
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)
# -------------------------------------------------------------------
# Set appropriate directories
# -------------------------------------------------------------------
setwd("~/Dropbox/projects/clarkiaScripts/data/reshapeData/")
# -------------------------------------------------------------------
# Load packages
# -------------------------------------------------------------------
library(dplyr) # manipulate data frames
library(tidyr) # manipulate data frames
library(reshape2) # manipulate data frames
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize data
# note; all these data come from permanet plots
# -------------------------------------------------------------------
# -------------------------------------------------------------------
seeds <- read.csv(file = "seeds.csv", header=TRUE)
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# visualize
# -------------------------------------------------------------------
# -------------------------------------------------------------------
par(mfrow=c(2,1))
hist(seeds$noSeedsPerFruit)
hist(seeds$noSeedsPerUndFruit,add=TRUE,col="red")
hist(seeds$noSeedsPerDamFruit,add=TRUE,col="blue")
hist(summarySeeds$meanSeedsPerFruitEq)
hist(summarySeeds$meanSeedsPerUndFruit,add=TRUE,col="red")
hist(summarySeeds$meanSeedsPerDamFruit,add=TRUE,col="blue")

plot(density(seeds$noSeedsPerFruit,na.rm=TRUE),ylim=c(0,0.08))
lines(density(seeds$noSeedsPerUndFruit,na.rm=TRUE),col="red")
lines(density(seeds$noSeedsPerDamFruit,na.rm=TRUE),col="blue")
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Gather
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# question whether or not to just include from permanent plot
df<-seeds %>% dplyr::select(-permanentPlot) %>% tidyr::gather(variable, response, -c(year,site))
df<-df[complete.cases(df), ]
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import predictor variables
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# site variables
climate <- read.csv(file = "allClimate.csv", header=TRUE) %>% 
  dplyr::filter(demography==1) %>%
  dplyr::select(-c(demography,T_summer,P_summer,T_springSummer,P_springSummer))
siteAbiotic <- read.csv(file = "siteAbiotic.csv", header=TRUE)
# make otherClarkia a factor (0 for no other Clarkia sp., 1 for other Clarkia sp.)
siteAbiotic$otherClarkia[siteAbiotic$otherClarkia>0] <- 1
siteAbiotic$otherClarkia<-as.factor(siteAbiotic$otherClarkia)
# rename abiotic variables
names(siteAbiotic)[6]="rock";names(siteAbiotic)[7]="aspect";
names(siteAbiotic)[8]="azimuth";names(siteAbiotic)[9]="slope";
allSitesAbiotic <- read.csv(file="allSitesAbiotic.csv",header=TRUE)
allSitesAbiotic <- dplyr::select(allSitesAbiotic,c(site,otherClarkiaSp))

allSitesAbiotic$S<-ifelse(grepl("S",allSitesAbiotic$otherClarkiaSp),1,0)
allSitesAbiotic$C<-ifelse(grepl("C",allSitesAbiotic$otherClarkiaSp),1,0)
allSitesAbiotic$U<-ifelse(grepl("U",allSitesAbiotic$otherClarkiaSp),1,0)
allSitesAbiotic<-allSitesAbiotic %>% dplyr::select(-otherClarkiaSp)

siteAbiotic<-siteAbiotic %>% dplyr::left_join(allSitesAbiotic,by="site")

# -------------------------------------------------------------------
# Create regression data frame 
# -------------------------------------------------------------------
regressionDf <- dplyr::select(df, c(site, year,  variable, response)) %>%
  dplyr::left_join(dplyr::select(siteAbiotic,
                   c(site,azimuth,slope,rock,otherClarkia,S,C,U)),
            by=c('site')) %>%
  dplyr::left_join(climate,
            by=c('site','year'))

# -------------------------------------------------------------------
# Clean up data for regression
# -------------------------------------------------------------------

# site, transect and position as factors
regressionDf$site <- as.factor(regressionDf$site)
regressionDf$rock <- as.factor(regressionDf$rock)
regressionDf$otherClarkia <- as.factor(regressionDf$otherClarkia)
regressionDf$S <- as.factor(regressionDf$S)
regressionDf$C <- as.factor(regressionDf$C)
regressionDf$U <- as.factor(regressionDf$U)

# -------------------------------------------------------------------
# Save data object
# -------------------------------------------------------------------
setwd("~/Dropbox/clarkiaSeedBanks/library/dataForAnalysis")
phiIndDF <- regressionDf
save(phiIndDF, file = "phiIndDF.RData")
