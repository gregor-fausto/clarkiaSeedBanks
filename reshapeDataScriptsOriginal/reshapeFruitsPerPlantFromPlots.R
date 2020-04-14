###########################################################################
# Required packages
###########################################################################
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(XLConnect)
library(xlsx)

###########################################################################
# Fruits per ...x... (2)
###########################################################################
rm(list=ls(all=TRUE)) # clear R environment

setwd("/Users/Gregor/Dropbox/Clarkia-LTREB/20_demography_sites")

filenames<-list.files(pattern=paste("^fruit_"), recursive=TRUE)
worksheet_2006 <- filenames[1]
worksheet_2007.2012 <- filenames[2:7]
worksheet_2013.2015 <- filenames[8:10]

# within a directory
# read sheet 'x' from workbooks 'y*`
# this will create a large data frame with merged values
readSheet <- function(fileList,sheetNumber){
  bound <- do.call("rbind.fill", lapply(fileList, function(x) 
    read.xlsx(file = x, sheetIndex = sheetNumber,
              startRow = 1,
              header = TRUE, FILENAMEVAR = x)))
}

df.2006 <- read.xlsx(file = worksheet_2006,
                     sheetIndex = 1,
                     colIndex = 1:2,
                     header = TRUE)
df.2006$year <- 2006
df.2006 <- df.2006 %>%
  dplyr::rename(noFruitsPerPlant = Frno_per_pl) %>%
  dplyr::rename(site = SITE)

# this will only work if column headers are same in all files
df_2007.2012 <- readSheet(fileList = worksheet_2007.2012,
                          sheetNumber = 1)

## total fruit equivalents 2007-2012
df.2007.2012 <- df_2007.2012[,1:7] %>%
  dplyr::rename(noFruitingPlants = fruitpl.) %>%
  dplyr::rename(meanNoFruitsPerPlant = mean_.fruit.pl) %>%
  dplyr::rename(noFruitsPerPlant = X..fruits.ind_plant) %>%
  dplyr::select(-meanNoFruitsPerPlant)

## total fruit equivalents 2013-2015
# this will only work if column headers are same in all files
df_2013.2015 <- readSheet(fileList = worksheet_2013.2015,
                          sheetNumber = 1)

## total fruit equivalents 2007-2012
df.2013.2015 <- df_2013.2015[,1:9] %>%
  dplyr::rename(noFruitingPlants = fruitpl.) %>%
  dplyr::rename(meanNoUndFruitsPerPlant = mean_.und_fruit.pl) %>%
  dplyr::rename(meanNoDamFruitsPerPlant = mean_.dam_fruit.pl) %>%
  dplyr::rename(noUndFruitsPerPlant = X.und_fruits.ind_plant) %>%
  dplyr::rename(noDamFruitsPerPlant = X.dam_fruits.ind_plant) %>%
  dplyr::select(-meanNoUndFruitsPerPlant) %>%
  dplyr::select(-meanNoDamFruitsPerPlant)

# clean up
df.2006$noUndFruitsPerPlant <- as.numeric(as.character("NA"))
df.2006$noDamFruitsPerPlant <- as.numeric(as.character("NA"))
df.2006$transect <- as.factor(NA)
df.2006$position <- as.numeric(as.character("NA"))

df.2007.2012$noUndFruitsPerPlant <- as.numeric(as.character("NA"))
df.2007.2012$noDamFruitsPerPlant <- as.numeric(as.character("NA"))

df.2013.2015$noFruitsPerPlant <- as.numeric(as.character("NA"))

df_full <- df.2006 %>%
  full_join(df.2007.2012) %>%
  full_join(df.2013.2015) 

###########################################################################
# Write data frames as .csv files
###########################################################################
setwd("/Users/Gregor/Dropbox/projects/clarkiaScripts/reshapeData/")
# fruits per plant
write.csv(df_full, file = "fruitsPerPlantFromPlots.csv", row.names=FALSE)