# -------------------------------------------------------------------
# Script to clean, organize, and reshape the data on 
# Counts for seedlings and number of fruiting plants
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)
# -------------------------------------------------------------------
# Set appropriate directories
# -------------------------------------------------------------------
setwd("~/Dropbox/projects/clarkiaScripts/data/reshapeData/")
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

# site, transect and position as factors
sigmaDF$uniqueTransect <- as.factor(sigmaDF$uniqueTransect)
sigmaDF$uniquePosition <- as.factor(sigmaDF$uniquePosition)

# -------------------------------------------------------------------
# Save data object
# -------------------------------------------------------------------
setwd("~/Dropbox/clarkiaSeedBanks/library/dataForAnalysis")
save(sigmaDF, file = "sigmaDF.RData")
