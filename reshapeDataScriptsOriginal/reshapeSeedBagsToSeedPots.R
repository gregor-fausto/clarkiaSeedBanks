rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)
# -------------------------------------------------------------------
# Loading libraryd packages
# -------------------------------------------------------------------
library(xlsx)
library(dplyr) # manipulate data frames
library(tidyr) # manipulate data frames
library(reshape2) # manipulate data frames
###########################################################################
# Read in data
###########################################################################
setwd("/Users/Gregor/Dropbox/Clarkia-LTREB/20_demography_sites")
seedBags<-read.xlsx(file="seed bag data.xls",sheetIndex = 1, header = TRUE)

# rename seed bag variables
seedBags<-seedBags %>%
  dplyr::rename(bagNo=bag..) %>%
  dplyr::rename(round=Round) %>%
  dplyr::rename(startYear=Start_yr.Oct.) %>%
  dplyr::rename(age=Age) %>%
  dplyr::rename(viabilityEstYear=Yr_germ.Jan..viability.Oct.) %>%
  dplyr::select(-c(viability.))

# age 1 bags
s1g1Bags<-seedBags %>% 
  dplyr::filter(age==1)
###########################################################################
# Survival and fecundity
###########################################################################
seedPots <- seedPots %>%
  dplyr::select(c(site,block,row,col,seedStart,yearStart,yearOne,seedlingYearOne))

###########################################################################
# Write data frames as .csv files
###################################################################### #####
setwd("/Users/Gregor/Dropbox/projects/clarkiaScripts/reshapeData/")
# survival
write.csv(seedPots, file = "seedPots.csv", row.names=FALSE)