# -------------------------------------------------------------------
# Script to clean, organize, and reshape the data on 
# Counts for seedlings and number of fruiting plants
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)
# -------------------------------------------------------------------
# Set appropriate directories
# -------------------------------------------------------------------
setwd("~/Dropbox/projects/clarkiaScripts/reshapeData/")
# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(dplyr) # manipulate data frames
library(tidyr) # manipulate data frames
library(reshape2) # manipulate data frames
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
df <- read.csv(file = "survival.csv", header=TRUE)
df <- df %>% 
  dplyr::select(-c(noTotalFruitEqUndFruitsPerPl,noDamFruitsPerPl))
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Create vector for failure to fruit
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# create vector for count of seedlings that failed to fruit
df <- df %>%
  dplyr::mutate(failToFruit = noSeedlings - noFruitingPlants)

# in 2006, there were plots that had no seedlings but 
# had fruiting plants
# recode failure to fruit=NA (to 0)
# this is equivalent to all plants surviving to fruiting (100% survival)
df$failToFruit[is.na(df$failToFruit) & !is.na(df$noFruitingPlants)] = 0

# in several years, there were plots that had more fruiting plants
# than flowering plants
# recode failure to fruit<0 (to 0)
# this is equivalent to all plants surviving to fruiting
df$failToFruit[df$failToFruit<0] = 0

# recode noSeedlings==0 & noFruitingPlants==0 (to NA)
# this will remove rows from the dataset where there were no
# no seedlings AND no fruiting plants in the same plot
df$failToFruit[df$noSeedlings==0 &
                          df$noFruitingPlants==0] = NA
df$noFruitingPlants[df$noSeedlings==0 &
                          df$noFruitingPlants==0] = NA

# keep only rows according to criteria above
df <- df %>% 
  dplyr::filter(!is.na(noFruitingPlants) & !is.na(failToFruit))

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import predictor variables
# -------------------------------------------------------------------
# -------------------------------------------------------------------
climate <- read.csv(file = "allClimate.csv", header=TRUE) %>% 
  dplyr::filter(demography==1) %>%
  dplyr::select(-c(demography,T_summer,P_summer,T_springSummer,P_springSummer))
siteAbiotic <- read.csv(file = "siteAbiotic.csv", header=TRUE)
# rename abiotic variables
names(siteAbiotic)[6]="rock";names(siteAbiotic)[7]="aspect";
names(siteAbiotic)[8]="azimuth";names(siteAbiotic)[9]="slope";
# -------------------------------------------------------------------
# Create regression data frame 
# -------------------------------------------------------------------
sigmaDF <- df %>%
  dplyr::left_join(dplyr::select(siteAbiotic,
                   c(site,azimuth,slope,rock)),
            by=c('site')) %>%
  dplyr::left_join(climate,
                   by=c('site','year'))
# -------------------------------------------------------------------
# Clean up data for regression
# -------------------------------------------------------------------
# remove column with number of seedlings
# sigmaDF <- sigmaDF[,-5]

# create unique transect/site variables
sigmaDF$uniqueTransect<-paste(sigmaDF$site,
                              sigmaDF$transect,sep="")
sigmaDF$uniquePosition<-paste(sigmaDF$site,
                              sigmaDF$transect,
                              sigmaDF$position,sep="")

# site, transect and position as factors
sigmaDF$site <- as.factor(sigmaDF$site)
sigmaDF$rock <- as.factor(sigmaDF$rock)
sigmaDF$uniqueTransect <- as.factor(sigmaDF$uniqueTransect)
sigmaDF$uniquePosition <- as.factor(sigmaDF$uniquePosition)

# -------------------------------------------------------------------
# Save data object
# -------------------------------------------------------------------
setwd("~/Dropbox/clarkiaSeedBanks/library/dataForAnalysis")
save(sigmaDF, file = "sigmaDF.RData")
