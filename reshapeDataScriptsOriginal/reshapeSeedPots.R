###########################################################################
# load packages
###########################################################################
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(xlsx)

###########################################################################
# Read in data
###########################################################################
setwd("/Users/Gregor/Dropbox/Clarkia-LTREB/20_demography_sites")
seedPots<-read.xlsx(file="pot_germination.xlsx",sheetIndex = 1, header = TRUE)
seedPots <- seedPots %>%
  dplyr::rename(site=Site) %>%
  dplyr::rename(block=Block) %>%
  dplyr::rename(row=Row) %>%
  dplyr::rename(col=Col) %>%
  dplyr::rename(seedStart=sd.) %>%
  dplyr::rename(yearStart=yr_pl) %>%
  dplyr::rename(yearEnd=yr_coll) %>%
  dplyr::rename(yearOne=census1) %>%
  dplyr::rename(seedlingYearOne=sdl._yr1) %>%
  dplyr::rename(yearTwo=census2) %>%
  dplyr::rename(seedlingYearTwo=sdl._yr2) %>%
  dplyr::rename(yearThree=census3) %>%
  dplyr::rename(seedlingYearThree=sdl._yr3)
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