# Libraries ---------------------------------------------------------------
library(tidyverse)
library(knitr)
library(ggplot2)
library(rjags)
library(tidybayes)
library(MCMCvis)

# Process seed bag data ---------------------------------------------------

# import data
seedBagsData <- readRDS("~/Dropbox/dataLibrary/postProcessingData/seedBagsData.rds")

# clean and organize data
datasetSize=nrow(seedBagsData)

# remove rows without January data
seedBagsData<-subset(seedBagsData,!is.na(seedBagsData$totalJan))
datasetSize <- cbind(datasetSize,nrow(seedBagsData))

# remove rows with more than 100 seeds initially
seedBagsData <- seedBagsData %>%
  dplyr::filter(totalJan<=100)
datasetSize <- cbind(datasetSize,nrow(seedBagsData))

# remove rows without October data
seedBagsData<-subset(seedBagsData,!is.na(seedBagsData$intactOct))
datasetSize <- cbind(datasetSize,nrow(seedBagsData))

# remove rows where there are more seeds in October than January
seedBagsData<-subset(seedBagsData,!(seedBagsData$intactJan<seedBagsData$intactOct))
datasetSize <- cbind(datasetSize,nrow(seedBagsData))

# Create a variable for the number of seeds at the start of the trial
seedBagsData$seedStart<-as.double(100)

# relabel seed bag number manually
seedBagsData[seedBagsData$site=="DLW" & seedBagsData$bagNo==43 & seedBagsData$age==1,]$bagNo = 42

saveRDS(seedBagsData,"~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/seedBagMaster.RDS")

# Import viability trial data
viabilityRawData <- readRDS("~/Dropbox/dataLibrary/postProcessingData/viabilityRawData.rds")

# clean and organize data
viabilityRawData <- viabilityRawData %>% 
  dplyr::select(-c(germPerc,germNot,viabPerc,viabPerc2,condTest))

viabilityRawData$bag <-as.integer(as.numeric(viabilityRawData$bagNo))

# recode rows with NA for the start of the viability trials to 0
viabilityRawData[is.na(viabilityRawData$viabStart),]$viabStart = 0

# remove data with incosistent sum of seeds across germination and viability trials
viabilityRawData %>% dplyr::filter(germStart - germCount - viabStart<0) 

viabilityRawData<-viabilityRawData %>% 
  dplyr::filter(germStart - germCount - viabStart >= 0)

# recode bags to account for duplicates (discovered by 'hand')
duplicates<-viabilityRawData %>%
  dplyr::select(site,round,age,bagNo) %>%
  unique()%>%
  dplyr::group_by(site,bagNo) %>%
  dplyr::summarise(n=n()) %>% 
  dplyr::filter(n>1) %>%
  tidyr::unite(col='id', c(site,bagNo), sep="-", remove=FALSE)
dim(duplicates)

duplicateBags<-viabilityRawData %>%
  tidyr::unite(col='id', c(site,bagNo), sep="-", remove=FALSE) %>%
  dplyr::filter(id %in% duplicates$id)

relabelBags <- function(site,round,age,bagNo,bagNoNew){
  x=viabilityRawData
  x[x$site==site & x$round==round & x$age==age & x$bagNo==bagNo,]$bagNo = bagNoNew
  return(x)
}

viabilityRawData=relabelBags("BR",1,3,9,71)
viabilityRawData=relabelBags("BR",1,2,23,72)
viabilityRawData=relabelBags("CF",1,3,64,71)
viabilityRawData=relabelBags("CP3",1,3,17,71)
viabilityRawData=relabelBags("GCN",2,2,50,71)
viabilityRawData=relabelBags("KYE",1,2,8,71)
viabilityRawData=relabelBags("LCE",1,3,3,71)
viabilityRawData=relabelBags("LO",1,3,13,78)
viabilityRawData=relabelBags("MC",1,3,10,71)
viabilityRawData=relabelBags("OKRE",1,3,5,71)
viabilityRawData=relabelBags("SM",1,3,3,71)

duplicates<-viabilityRawData %>%
  dplyr::select(site,round,age,bagNo) %>%
  unique()%>%
  dplyr::group_by(site,bagNo) %>%
  dplyr::summarise(n=n()) %>% 
  dplyr::filter(n>1) %>%
  tidyr::unite(col='id', c(site,bagNo), sep="-", remove=FALSE)
dim(duplicates)

# also remove bags that are in the viability dataset but not the seed bag dataset
vComp=viabilityRawData %>% 
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="-", remove=FALSE)
sComp=seedBagsData %>% dplyr::select(site,round,age,bagNo) %>%
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="-", remove=FALSE)

viabilityRawData<-viabilityRawData %>% 
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="-", remove=FALSE) %>%
  dplyr::filter((vComp$id %in% sComp$id)) %>% 
  dplyr::select(-id)

# these are no longer the "raw data" but the ones where bags have been relabeled
saveRDS(viabilityRawData,"~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/viabilityFinal.RDS")

# index for JAGS
# filterData<-function(x, ageNumber) {
#   x %>%
#     dplyr::filter(age==ageNumber)
# }

# Subset data by age
# seedBagsData1<-filterData(seedBagsData,ageNumber=1)
# viabilityRawData1<-filterData(viabilityRawData,ageNumber=1)

# seedBagsData2<-filterData(seedBagsData,ageNumber=2)
# viabilityRawData2<-filterData(viabilityRawData,ageNumber=2)
# 
# seedBagsData3<-filterData(seedBagsData,ageNumber=3)
# viabilityRawData3<-filterData(viabilityRawData,ageNumber=3)

# age 1 data
seedBagsData<-seedBagsData %>%
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="-", remove=FALSE) %>%
  tidyr::unite(col='siteBag', c(site,bagNo), sep="-", remove=FALSE) %>%
  dplyr::mutate(siteBag = as.factor(siteBag)) 

viabilityRawData<-viabilityRawData %>%
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="-", remove=FALSE) %>%
  tidyr::unite(col='siteBag', c(site,bagNo), sep="-", remove=FALSE) %>%
  dplyr::mutate(siteBag = as.factor(siteBag)) 

# link these via a reference table
referenceTable<-data.frame(id=union(seedBagsData$id, viabilityRawData$id)) %>%
  dplyr::mutate(idNo = 1:length(id)) 

seedBagsData<-seedBagsData %>%
  dplyr::left_join(referenceTable,by="id")

viabilityRawData<-viabilityRawData %>%
  dplyr::left_join(referenceTable,by="id")

seedBagsData = seedBagsData %>%
  dplyr::mutate(year = as.factor(yearStart),
                age = as.factor(age)) %>%
  dplyr::select(site,year,age,totalJan,seedStart,seedlingJan,intactOct) %>%
  dplyr::rename(siteBags = site,
                yearBags = year,
                ageBags = age)

viabilityRawData = viabilityRawData %>%
  dplyr::mutate(year = as.factor(round),
                age = as.factor(age)) %>%
  dplyr::select(id, site, year, age, germStart, germCount, viabStart, viabStain, idNo) %>%
  dplyr::rename(siteViab = site,
                yearViab = year,
                ageViab = age,
                bag = idNo) %>%
  dplyr::mutate(germStart = ifelse(is.na(germCount), NA, germStart),
                viabStart = ifelse(is.na(viabStain), NA, viabStart)) %>% 
  dplyr::group_by(id, siteViab,yearViab, ageViab, bag) %>% 
  # sum observations in each bag; this is ignoring some variation
  dplyr::summarise(germStart = sum(germStart,na.rm=TRUE),
                   germCount = sum(germCount,na.rm=TRUE),
                   viabStart = sum(viabStart,na.rm=TRUE),
                   viabStain = sum(viabStain,na.rm=TRUE))

# # Age 2 data
# seedBagsData2<-seedBagsData2 %>%
#   tidyr::unite(col='id', c(site,bagNo,round,age), sep="-", remove=FALSE) %>%
#   tidyr::unite(col='siteBag', c(site,bagNo), sep="-", remove=FALSE) %>%
#   dplyr::mutate(siteBag = as.factor(siteBag)) 
# 
# viabilityRawData2<-viabilityRawData2 %>%
#   tidyr::unite(col='id', c(site,bagNo,round,age), sep="-", remove=FALSE) %>%
#   tidyr::unite(col='siteBag', c(site,bagNo), sep="-", remove=FALSE) %>%
#   dplyr::mutate(siteBag = as.factor(siteBag)) 
# 
# referenceTable<-data.frame(id=union(seedBagsData2$id, viabilityRawData2$id)) %>%
#   dplyr::mutate(idNo = 1:length(id)) 
# 
# seedBagsData2<-seedBagsData2 %>%
#   dplyr::left_join(referenceTable,by="id")
# 
# viabilityRawData2<-viabilityRawData2 %>%
#   dplyr::left_join(referenceTable,by="id")
# 
# seedBagsData2 = seedBagsData2 %>%
#   dplyr::mutate(year = as.factor(yearStart),
#                 age = as.factor(age)) %>%
#   dplyr::select(site,year,age,totalJan,seedStart,seedlingJan,intactOct) %>%
#   dplyr::rename(siteBags = site,
#                 yearBags = year,
#                 ageBags = age)
# 
# viabilityRawData2 = viabilityRawData2 %>%
#   dplyr::mutate(year = as.factor(round),
#                 age = as.factor(age)) %>%
#   dplyr::select(site, year, age, germStart, germCount, viabStart, viabStain, idNo) %>%
#   dplyr::rename(siteViab = site,
#                 yearViab = year,
#                 ageViab = age,
#                 bag = idNo) %>%
#   #dplyr::mutate(bag = as.factor(bag)) %>%
#   dplyr::group_by(siteViab,yearViab, ageViab, bag) %>%
#   dplyr::mutate(germStart = ifelse(is.na(germCount), NA, germStart),
#                 viabStart = ifelse(is.na(viabStain), NA, viabStart)) %>%
#   # sum observations in each bag; this is ignoring some variation
#   dplyr::summarise(germStart = sum(germStart),
#                    germCount = sum(germCount),
#                    viabStart = sum(viabStart),
#                    viabStain = sum(viabStain))
# 
# names(seedBagsData2)=paste(names(seedBagsData2),"2",sep="")
# names(viabilityRawData2)=paste(names(viabilityRawData2),"2",sep="")

# pass to JAGS
viabilityRawData<-viabilityRawData[complete.cases(viabilityRawData),] %>% dplyr::mutate(bag = as.factor(bag))
#viabilityRawData2<-viabilityRawData2[complete.cases(viabilityRawData2),] %>% dplyr::mutate(bag2 = as.factor(bag2))

viabilityRawDataRefTable<-viabilityRawData %>%
  dplyr::left_join(referenceTable,by="id")

saveRDS(viabilityRawDataRefTable,"~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/viabilityBagRefTable.RDS")

viabilityRawData <- viabilityRawData %>%
  ungroup %>% 
  dplyr::select(-bag) %>%
  dplyr::rename(bag=id)

# filter to one site for testing
#seedBagsData <- seedBagsData %>% dplyr::filter(siteBags=="S22")
#viabilityRawData <- viabilityRawData %>% dplyr::filter(siteViab=="S22") 
#viabilityRawData <- viabilityRawData %>% dplyr::filter(siteViab2=="S22")

saveRDS(viabilityRawData,"~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/viabilityRawData.RDS")
saveRDS(seedBagsData,"~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/seedBagsRawData.RDS")

data <- tidybayes::compose_data(seedBagsData,viabilityRawData)

data$n = dim(seedBagsData)[1]
data$n_bag  = length(data$bag)
data$n_total = 6
data$roundIndex = c(1,1,1,2,2,3)

data$round= ifelse(data$yearViab==1&data$ageViab==1,1,
         ifelse(data$yearViab==2&data$ageViab==1,2,
                ifelse(data$yearViab==3&data$ageViab==1,3,
                       ifelse(data$yearViab==1&data$ageViab==2,4,
                              ifelse(data$yearViab==2&data$ageViab==2,5,6)))))

data$roundBags = ifelse(data$yearBags==1&data$ageBags==1,1,
                   ifelse(data$yearBags==2&data$ageBags==1,2,
                          ifelse(data$yearBags==3&data$ageBags==1,3,
                                 ifelse(data$yearBags==1&data$ageBags==2,4,
                                        ifelse(data$yearBags==2&data$ageBags==2,5,6)))))

saveRDS(data,"~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/viabilityData.RDS")

# 
# # Viability ----------------------------------------------------------------
# 
# viabilityRawData1 %>%
#   dplyr::mutate(pv = viabStain/viabStart,
#                 ng = germStart-germCount,
#                 pvtot = (germStart-ng*(pv-1))/germStart)
# 
# viabilityRawData1 %>%
#   dplyr::group_by(yearViab) %>%
#   dplyr::mutate(mu.tj = mean(totalJan)/100,
#                 mu.sd = mean(seedlingJan)/100,
#                 mu.io = mean(intactOct)/100) %>%
#   dplyr::mutate(p.tj = totalJan/seedStart,
#                 p.sd = seedlingJan/seedStart,
#                 p.io = intactOct/seedStart)
# 
# 
# # Plotting ----------------------------------------------------------------
# 
# par(mfrow=c(1,3))
# plot(seedBagsData1$yearBags,seedBagsData1$totalJan,xlab="Year",ylab="Total January")
# plot(seedBagsData1$yearBags,seedBagsData1$seedlingJan,xlab="Year",ylab="Seedlings January")
# plot(seedBagsData1$yearBags,seedBagsData1$intactOct,xlab="Year",ylab="Intact October")
# 
# age1Summary <- seedBagsData1 %>%
#   dplyr::group_by(yearBags) %>%
#   dplyr::mutate(mu.tj = mean(totalJan)/100,
#                 mu.sd = mean(seedlingJan)/100,
#                 mu.io = mean(intactOct)/100) %>%
#   dplyr::mutate(p.tj = totalJan/seedStart,
#                 p.sd = seedlingJan/seedStart,
#                 p.io = intactOct/seedStart)
# 
# par(mfrow=c(1,3))
# plot(age1Summary$mu.tj,age1Summary$p.tj,
#      xlim=c(0,1),ylim=c(0,1),pch=16,cex=.75,
#      xlab="Mean prop January total",
#      ylab="Observed prop January total")
# abline(a=0,b=1)
# points(unique(age1Summary$mu.tj),unique(age1Summary$mu.tj),pch=1,col='red')
# 
# plot(age1Summary$mu.sd,age1Summary$p.sd,
#      xlim=c(0,.5),ylim=c(0,.5),pch=16,cex=.75,
#      xlab="Mean prop January seedlings",
#      ylab="Observed prop January seedlings")
# abline(a=0,b=1)
# points(unique(age1Summary$mu.sd),unique(age1Summary$mu.sd),pch=1,col='red')
# 
# plot(age1Summary$mu.io,age1Summary$p.io,
#      xlim=c(0,1),ylim=c(0,1),pch=16,cex=.75,
#      xlab="Mean prop October intact",
#      ylab="Observed prop October intact")
# abline(a=0,b=1)
# points(unique(age1Summary$mu.io),unique(age1Summary$mu.io),pch=1,col='red')
# 
# 
# seedBagsData1 %>%
#   dplyr::group_by(yearBags) %>%
#   dplyr::summarise(mu.tj = mean(totalJan)/100,
#                    mu.sd = mean(seedlingJan)/100,
#                    mu.io = mean(intactOct)/100) 
# 
