rm(list=ls(all=TRUE)) # clear R environment

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

setwd("/Users/Gregor/Dropbox/Clarkia-LTREB/20_demography_sites")

filenames<-list.files(pattern=paste("^fruit&seed"), recursive=TRUE)
worksheet_2006.2012 <- filenames[1:7]
worksheet_2013.2014 <- filenames[8:9]

# within a directory
# read sheet 'x' from workbooks 'y*`
# this will create a large data frame with merged values
readSheet <- function(fileList,sheetNumber){
  bound <- do.call("rbind.fill", lapply(fileList, function(x) 
    read.xlsx(file = x, sheetIndex = sheetNumber,
              startRow = 1,colIndex=1:2,
              header = TRUE, FILENAMEVAR = x)))
}

# this will only work if column headers are same in all files
df_2006.2012 <- readSheet(fileList = worksheet_2006.2012,
                          sheetNumber = 2) 

df_2006.2012<-cbind(dplyr::select(df_2006.2012,FILENAMEVAR),dplyr::select(df_2006.2012,-FILENAMEVAR))
yrs <- substring(df_2006.2012[,1],17,20)
df_2006.2012$year <- as.numeric(yrs)

df_2006.2012<-df_2006.2012 %>%
  dplyr::rename(noFruitsPerPlant = Fruit.number.per.plant) %>%
  dplyr::select(-FILENAMEVAR)

# within a directory
# read sheet 'x' from workbooks 'y*`
# this will create a large data frame with merged values
readSheet <- function(fileList,sheetNumber){
  bound <- do.call("rbind.fill", lapply(fileList, function(x) 
    read.xlsx(file = x, sheetIndex = sheetNumber,
              startRow = 1,colIndex=1:3,
              header = TRUE, FILENAMEVAR = x)))
}

# this will only work if column headers are same in all files
df_2013.2014 <- readSheet(fileList = worksheet_2013.2014,
                          sheetNumber = 2)

df_2013.2014<-cbind(dplyr::select(df_2013.2014,FILENAMEVAR),dplyr::select(df_2013.2014,-FILENAMEVAR))
yrs <- substring(df_2013.2014[,1],17,20)
df_2013.2014$year <- as.numeric(yrs)

df_2013.2014<-df_2013.2014 %>%
  dplyr::rename(noUndFruitsPerPlant = Fruit.number.per.plant) %>%
  dplyr::rename(noDamFruitsPerPlant = NA.)  %>%
  dplyr::select(-FILENAMEVAR)

df_full <- df_2006.2012 %>%
  dplyr::full_join(df_2013.2014[-(1),]) %>%
  dplyr::rename(site = Site)

###########################################################################
# Write data frames as .csv files
###########################################################################
setwd("/Users/Gregor/Dropbox/projects/clarkiaScripts/reshapeData/")
# fruits per plant
write.csv(df_full, file = "fruitsPerPlantAllData.csv", row.names=FALSE)