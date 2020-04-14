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
# Seeds
###########################################################################
setwd("/Users/Gregor/Dropbox/Clarkia-LTREB/20_demography_sites")

filenames<-list.files(pattern=paste("^seeds_"), recursive=TRUE)
worksheet_2006.2010 <- filenames[1:5]
worksheet_2011.2012 <- filenames[6:7]
worksheet_2013.2015 <- filenames[8:10]

# within a directory
# read sheet 'x' from workbooks 'y*`
# this will create a large data frame with merged values
readSheet <- function(fileList,sheetNumber,startingRow,endColumn){
  bound <- do.call("rbind.fill", lapply(fileList, function(x) 
    read.xlsx(file = x, sheetIndex = sheetNumber,
              startRow = startingRow, colIndex = 1:endColumn,
              header = TRUE, FILENAMEVAR = x)))
}

# reshape 2006-2010
df_2006.2010 <- readSheet(fileList = worksheet_2006.2010, 
                          sheetNumber=1, startingRow = 1,
                          endColumn = 2)

df.2006.2010<-cbind(select(df_2006.2010,FILENAMEVAR),select(df_2006.2010,-FILENAMEVAR))
yrs <- substring(df.2006.2010[,1],7,10)
df.2006.2010$FILENAMEVAR <- as.numeric(yrs)

df.2006.2010 <- df.2006.2010 %>%
  dplyr::rename(year = FILENAMEVAR) %>%
  dplyr::rename(noSeedsPerFruit = Seed.no.per.fruit) %>%
  dplyr::rename(site = Site) %>%
  dplyr::filter(site != "<NA>")
  
df.2006.2010$permanentPlot <- as.numeric(1)

# reshape 2011-2012
df_2011.2012 <- readSheet(fileList = worksheet_2011.2012, 
                          sheetNumber=1, startingRow = 1,
                          endColumn = 3)

df.2011.2012<-cbind(select(df_2011.2012,FILENAMEVAR),select(df_2011.2012,-FILENAMEVAR))
yrs <- substring(df.2011.2012[,1],7,10)
df.2011.2012$FILENAMEVAR <- as.numeric(yrs)

df.2011.2012 <- df.2011.2012 %>%
  dplyr::rename(year = FILENAMEVAR) %>%
  dplyr::rename(noSeedsPerFruit = Seed.no.per.fruit) %>%
  dplyr::rename(site = Site) %>%
  dplyr::filter(site != "<NA>")

df.2011.2012 <- df.2011.2012 %>%
  dplyr::filter(permanentPlot==1)

# reshape 2013-2015
df_2013.2015_undamaged <- readSheet(fileList = worksheet_2013.2015, 
                          sheetNumber=1, startingRow = 1,
                          endColumn = 3)

df_2013.2015_damaged <- readSheet(fileList = worksheet_2013.2015, 
                                    sheetNumber=2, startingRow = 1,
                                    endColumn = 2)

df.2013.2015_undamaged<-cbind(select(df_2013.2015_undamaged,FILENAMEVAR),
                              select(df_2013.2015_undamaged,-FILENAMEVAR))
yrs <- substring(df.2013.2015_undamaged[,1],7,10)
df.2013.2015_undamaged$FILENAMEVAR <- as.numeric(yrs)

df.2013.2015_undamaged <- df.2013.2015_undamaged %>%
  dplyr::rename(year = FILENAMEVAR) %>%
  dplyr::rename(noSeedsPerUndFruit = Seed.no.per.undamaged.fruit) %>%
  dplyr::rename(site = Site) %>%
  dplyr::filter(site != "<NA>")

df.2013.2015_undamaged <- df.2013.2015_undamaged %>%
  dplyr::filter(permanentPlot==1)
  
df.2013.2015_undamaged$noSeedsPerDamFruit <- as.numeric(as.character("NA"))
df.2013.2015_undamaged$noSeedsPerUndFruit <- as.numeric(as.character(df.2013.2015_undamaged$noSeedsPerUndFruit))
df.2013.2015_undamaged$permanentPlot <- as.numeric(df.2013.2015_undamaged$permanentPlot)

df.2013.2015_damaged<-cbind(select(df_2013.2015_damaged,FILENAMEVAR),
                              select(df_2013.2015_damaged,-FILENAMEVAR))
yrs <- substring(df.2013.2015_damaged[,1],7,10)
df.2013.2015_damaged$FILENAMEVAR <- as.numeric(yrs)

df.2013.2015_damaged <- df.2013.2015_damaged %>%
  dplyr::rename(year = FILENAMEVAR) %>%
  dplyr::rename(noSeedsPerDamFruit = Seed.no.per.damaged.fruit) %>%
  dplyr::rename(site = Site) %>%
  dplyr::filter(site != "<NA>")
  
df.2013.2015_damaged$noSeedsPerUndFruit <- as.numeric(as.character("NA"))
df.2013.2015_damaged$noSeedsPerDamFruit <- as.numeric(as.character(df.2013.2015_damaged$noSeedsPerDamFruit))
df.2013.2015_damaged$permanentPlot <- as.numeric(1)

df.2013.2015 <- df.2013.2015_undamaged %>%
  dplyr::full_join(df.2013.2015_damaged)

df.2013.2015_damaged$noSeedsPerFruit <- as.numeric(as.character("NA"))

# clean up
df.2011.2012<-df.2011.2012[,1:3]
df_full <- df.2006.2010 %>%
  full_join(df.2011.2012)
df_full$noSeedsPerUndFruit <- as.numeric(as.character("NA"))
df_full$noSeedsPerDamFruit <- as.numeric(as.character("NA"))

df_full <- df_full %>%
dplyr::full_join(df.2013.2015)

###########################################################################
# Write data frames as .csv files
###########################################################################
setwd("/Users/Gregor/Dropbox/projects/clarkiaScripts/reshapeData/")
# fruits per plants - ALL plots
write.csv(df_full, file = "seeds.csv", row.names = FALSE)