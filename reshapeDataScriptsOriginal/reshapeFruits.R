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
# Fruits per ...x... (1)
###########################################################################
setwd("/Users/Gregor/Dropbox/Clarkia-LTREB/20_demography_sites")

bindSheets <- function(fileList,sheetNumber,numberOfColumns){
  bound <- do.call("rbind", lapply(fileList, function(x) 
    read.xlsx(file = x, sheetIndex = sheetNumber, 
              colIndex = (1:numberOfColumns), 
              header = TRUE, FILENAMEVAR = x)))
}

filenames<-list.files(pattern=paste("^fruit_data_"), recursive=TRUE)

# 2006
fruit_2006 <- read.xlsx(file="fruit_2006.xlsx",sheetIndex = 1, header = TRUE)

# 2007-2012
files <- filenames[1:6]
fruit_2007.2012 <- bindSheets(fileList = files,
                              sheetNumber = 1,
                              numberOfColumns = 11)

# 2013-2015
files <- c(filenames[7:9])
fruit_2013.2015 <- bindSheets(fileList = files,
                              sheetNumber = 1,
                              numberOfColumns = 13)

df_2006 <- fruit_2006
df_2006$year <- 2006
df_2006 <- df_2006[c("year", "SITE", "Frno_per_pl")]
df_2006 <- df_2006 %>%
  dplyr::rename(noFruitsPerPlant = Frno_per_pl) %>%
  dplyr::rename(site = SITE)

df_2007.2012 <- fruit_2007.2012 %>%
  dplyr::select(year:X..fruits.ind_plant) %>%
  dplyr::rename(noFruitingPlants = fruitpl.) %>%
  dplyr::rename(meanNoFruitsPerPlant = mean_.fruit.pl) %>%
  dplyr::rename(noFruitsPerPlant = X..fruits.ind_plant)

df_2013.2015 <- fruit_2013.2015 %>%
  dplyr::select(year:X.dam_fruits.ind_plant) %>%
  dplyr::rename(noFruitingPlants = fruitpl.) %>%
  dplyr::rename(meanNoUndFruitsPerPlant = mean_.und_fruit.pl) %>%
  dplyr::rename(meanNoDamFruitsPerPlant = mean_.dam_fruit.pl) %>%
  dplyr::rename(noUndFruitsPerPlant = X.und_fruits.ind_plant) %>%
  dplyr::rename(noDamFruitsPerPlant = X.dam_fruits.ind_plant)

df_fruits <- df_2006 %>%
  full_join(df_2007.2012) %>%
  full_join(df_2013.2015)

df_fruits <- df_fruits[c("year", "site", "transect", "position", 
                         "noFruitsPerPlant","noUndFruitsPerPlant","noDamFruitsPerPlant",
                         "noFruitingPlants",
                         "meanNoFruitsPerPlant","meanNoUndFruitsPerPlant","meanNoDamFruitsPerPlant")]

###########################################################################
# Write data frames as .csv files
###################################################################### #####
setwd("/Users/Gregor/Dropbox/projects/clarkiaScripts/reshapeData/")
# fruits 
write.csv(df_fruits, file = "fruits.csv", row.names=FALSE)