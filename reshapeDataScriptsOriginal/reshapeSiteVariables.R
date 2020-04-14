###########################################################################
# library(d packages
###########################################################################
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
#library(XLConnect)
library(xlsx)

###########################################################################
# Read in data
###########################################################################
setwd("/Users/Gregor/Dropbox/Clarkia-LTREB/20_demography_sites")
surv_fec<-read.xlsx(file="Surv_Fec_env_06-15.xls",sheetIndex = 1, header = TRUE)
#surv_fec <- read.csv("Surv_Fec_env_06-15.csv", header=TRUE)
surv_fec <- surv_fec %>%
  dplyr::rename(year=year_spring)
#biotic<-read.xlsx(file="20_site_env_climate_2005-2015.xls",sheetIndex = 1, header = TRUE,startRow=2)

setwd("/Users/Gregor/Dropbox/Clarkia-LTREB/41_sites")
abioticVars<-read.xlsx(file="41_site_env_climate_2005-2015_complete.xlsx",
                       sheetIndex = 1, header = TRUE,startRow=2,colIndex=1:16)
abioticVarsMissing<-read.csv(file="abioticData41Sites.csv",header=TRUE)

###########################################################################
# Site-Level Variables (names, abiotic variables, climate)
###########################################################################
siteLevel <- surv_fec %>%
  dplyr::select(site:P_spring)
  
siteYearList <- siteLevel %>%
  dplyr::select(c(site,year,transect,position)) %>%
  unique

siteAbiotic <- siteLevel %>%
  dplyr::select(c(site,easting,northing,elevation,
                  area,dominant.surface.rock.type,
                  raw.aspect,linear_azimuth,average_slope,
                  spring.equinox.rad)) %>%
  unique %>%
  dplyr::slice(1:20)

siteAbiotic$otherClarkia<-biotic[,4]   

siteClimate <- siteLevel %>%
  dplyr::select(c(site,year,T_winter,T_spring,P_winter,P_spring)) %>%
  unique() 
siteClimate<-siteClimate[!duplicated(siteClimate[1:2]),]

allSitesAbiotic <- abioticVars %>%
  dplyr::select(SITE:spring.equinox.rad)

names(allSitesAbiotic) <- c("site","demography","otherClarkia",
                            "otherClarkiaSp","admixture","easting",
                            "northing","elevation","area","rock","aspect",
                            "azimuth","slope","winterSolsticeRad","springEquinoxRad")

# replace missing data
siteAbiotic$average_slope[9] <- abioticVarsMissing$average_slope[30]
siteAbiotic$raw.aspect[17] <- abioticVarsMissing$raw.aspect[38]

allSitesAbiotic$slope[30] <- abioticVarsMissing$average_slope[30]
allSitesAbiotic$aspect[38] <- abioticVarsMissing$raw.aspect[38]

# missingData <- abioticVars %>%
#   dplyr::select(c(site,raw.aspect,linear_azimuth,average_slope))
# 
# write.csv(missingData,file="missingData.csv")

###########################################################################
# Write data frames as .csv files
###########################################################################
setwd("/Users/Gregor/Dropbox/projects/clarkiaScripts/reshapeData/")
# site-level abiotic variables
write.csv(siteAbiotic, file = "siteAbiotic.csv",row.names=FALSE)
# site*year climate variables
write.csv(siteClimate, file = "siteClimate.csv", row.names=FALSE)
# all site-level abiotic variables
write.csv(allSitesAbiotic, file= "allSitesAbiotic.csv",row.names=FALSE)
