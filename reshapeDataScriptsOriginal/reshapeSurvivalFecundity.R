###########################################################################
# load packages
###########################################################################
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(XLConnect)
library(xlsx)

###########################################################################
# Read in data
###########################################################################
setwd("/Users/Gregor/Dropbox/Clarkia-LTREB/20_demography_sites")
surv_fec<-read.xlsx(file="Surv_Fec_env_06-15.xls",sheetIndex = 1, header = TRUE)
surv_fec <- surv_fec %>%
  dplyr::rename(year=year_spring)

###########################################################################
# Survival and fecundity
###########################################################################
survival <- surv_fec %>%
  select(c(site,year,transect,position,seedling.,fruitpl.,tot_or_undamaged_frperpl.,damaged_frpl.))

survival <- survival %>%
  dplyr::rename(noSeedlings = seedling.) %>%
  dplyr::rename(noFruitingPlants = fruitpl.) %>%
  dplyr::rename(noTotalFruitEqUndFruitsPerPl = tot_or_undamaged_frperpl.) %>%
  dplyr::rename(noDamFruitsPerPl = damaged_frpl.)

###########################################################################
# Write data frames as .csv files
###################################################################### #####
setwd("/Users/Gregor/Dropbox/projects/clarkiaScripts/reshapeData/")
# survival
write.csv(survival, file = "survival.csv", row.names=FALSE)