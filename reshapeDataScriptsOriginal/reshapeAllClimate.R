rm(list=ls(all=TRUE)) # clear R environment
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import climate
# -------------------------------------------------------------------
# -------------------------------------------------------------------
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
# Read in data
###########################################################################
setwd("~/Dropbox/Clarkia-LTREB/41_sites")
siteEnv<-read.xlsx(file="41_site_env_climate_2005-2016_complete.xlsx",sheetIndex = 1, header = TRUE,startRow=1)

siteEnv <- siteEnv %>%
  dplyr::select(c(SITE,intense.demography.,Twinter05.06:Pspring16)) %>%
  dplyr::rename(site=SITE,demography=intense.demography.)

envMolten<-melt(siteEnv, id.vars = c("site", "demography"))

year<-rep(rep(2006:2016,each=41*3)[1:1312],2)
season<-c(rep(rep(c("T_winter","T_spring","T_summer"),each=41),
              11)[1:1312],
          rep(rep(c("P_winter","P_spring","P_summer"),each=41),
              11)[1:1312])

reboundEnv <- cbind(envMolten,year,season)

order <- as.vector(c("site","demography",
                     "T_winter","T_spring","T_summer",
                     "P_winter","P_spring","P_summer"))

allClimate <- reboundEnv %>%
  dplyr::select(-variable) %>%
  tidyr::spread(season,value) 

allClimate<-allClimate[,c(1,2,3,9,7,8,6,4,5)]
allClimate<-allClimate %>%
  dplyr::mutate(T_springSummer=((T_spring*150)+(T_summer*122))/272) %>%
  dplyr::mutate(P_springSummer=(P_spring+P_summer))

###########################################################################
# Write data frames as .csv files
###########################################################################
setwd("~/Dropbox/projects/clarkiaScripts/reshapeData")
write.csv(allClimate,file="allClimate.csv",row.names=FALSE)
