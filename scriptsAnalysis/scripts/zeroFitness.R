# -------------------------------------------------------------------
# Checking individual models
# to look for missing data values
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)
# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(rjags) # jags interface
library(MCMCvis)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(HDInterval)
library(bayesplot)
library(stringr)

set.seed(10)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Seedling survival to fruiting
# -------------------------------------------------------------------
# -------------------------------------------------------------------
censusSeedlingsFruitingPlants <- readRDS("~/Dropbox/dataLibrary/postProcessingData/censusSeedlingsFruitingPlants.RDS")

censusSeedlingsFruitingPlants <- censusSeedlingsFruitingPlants %>% 
  # more fruiting plants than seedlings
  # recode s.t. number of seedlings equals number of fruiting plants
  dplyr::mutate(seedlingNumber=ifelse(fruitplNumber>seedlingNumber,fruitplNumber,seedlingNumber)) %>% 
  # NA for seedlings
  # filter these out; these are true missing data
  # recode s.t. number of seedlings equals number of fruiting plants
  #dplyr::mutate(seedlingNumber=ifelse(is.na(seedlingNumber),fruitplNumber,seedlingNumber)) %>%
  dplyr::filter(!is.na(seedlingNumber)) %>% 
  # NA for fruiting plants
  # still filter these out; missing response
  dplyr::filter(!is.na(fruitplNumber)) 

sigma.df = censusSeedlingsFruitingPlants %>%
  dplyr::mutate( p = fruitplNumber/seedlingNumber )# %>% 
 # dplyr::filter( p <= 1 | is.na(p))

sigma.noObs = sigma.df %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n.obs = n(),
                   n.na = sum(is.na(p)),
                   n.zero = sum(p==0|is.na(p)))

sigmaSummary = sigma.noObs %>%
  dplyr::mutate(trueNA = ifelse(n.obs==n.na,1,0)) %>%
  dplyr::mutate(obsZero = ifelse(n.obs==n.zero,1,0))

siteNames=unique(sigmaSummary$site)

plot(NA,NA,xlim=c(2006,2020),ylim=c(0,20),
     axes=FALSE,frame=FALSE,xaxt='n',yaxt='n',
     xlab='',ylab='')
y.pt=20:1
for(i in 20:1){
  tmp=sigmaSummary[sigmaSummary$site==siteNames[i],]
  points(tmp$year,y=rep(y.pt[i],length(tmp$year)),
         pch=ifelse(tmp$trueNA==0,19,4),
         col=ifelse(tmp$obsZero==1,"gray","black"))
}
axis(1,  2006:2020, col.ticks = 1,las=2)
axis(2, (1:20),
     labels = rev(siteNames), las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)

zeroFitness = sigmaSummary %>% dplyr::filter(obsZero==1) %>%
  dplyr::mutate(site.year = paste0(site,year))
zeroFitness.df=sigma.df %>%
  dplyr::mutate(site.year = paste0(site,year)) %>%
  dplyr::filter(site.year %in% zeroFitness$site.year) %>% 
  dplyr::filter(p==0) 

f = function(x){
  max(seq(0,1,length.out=1000)[dbinom(0,x,seq(0,1,length.out=1000))>.5])
}

combos=unique(zeroFitness.df$site.year)
combos.list = list()
for(i in 1:length(combos)){
  tmp = zeroFitness.df[zeroFitness.df$site.year==combos[i],]
  combos.list[[i]] = sapply(tmp$seedlingNumber,f)
}
prob.sigma=lapply(combos.list, prod )



# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Fruits per plant
# -------------------------------------------------------------------
# -------------------------------------------------------------------
countFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countFruitsPerPlantAllPlots.RDS")

countFruitsPerPlantAllPlots <- countFruitsPerPlantAllPlots %>%
  dplyr::rename(y_tfe = countFruitNumberPerPlant) %>%
  dplyr::select(site,year,y_tfe)

countFruitsPerPlantAllPlots$year <- as.character(countFruitsPerPlantAllPlots$year)

countUndamagedDamagedFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countUndamagedDamagedFruitsPerPlantAllPlots.RDS")

countUndamagedDamagedFruitsPerPlantAllPlots <- countUndamagedDamagedFruitsPerPlantAllPlots %>%
  dplyr::rename(y_und = countUndamagedFruitNumberPerPlant) %>%
  dplyr::rename(y_dam = countDamagedFruitNumberPerPlant) %>%
  dplyr::rename(site2 = site) %>%
  dplyr::rename(year2 = year) %>%
  dplyr::select(site2,year2,y_und,y_dam) 

countUndamagedDamagedFruitsPerPlantAllPlots$year2 <- as.character(countUndamagedDamagedFruitsPerPlantAllPlots$year2)

fec.df = countFruitsPerPlantAllPlots

fec.noObs = fec.df %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n.obs = n(),
                   n.na = sum(is.na(y_tfe)),
                   n.zero = sum(y_tfe==0|is.na(y_tfe)))

fecSummary1 = fec.noObs %>%
  dplyr::mutate(trueNA = ifelse(n.obs==n.na,1,0)) %>%
  dplyr::mutate(obsZero = ifelse(n.obs==n.zero,1,0)) %>%
  dplyr::select(site,year,n.obs,trueNA,obsZero)

fec.df = countUndamagedDamagedFruitsPerPlantAllPlots %>%
  dplyr::rename(site=site2,year=year2)

fec.noObs = fec.df %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n.obs = n(),
                   n.na1 = sum(is.na(y_und)),n.na2=sum(is.na(y_dam)),
                   n.zero1 = sum(y_und==0|is.na(y_und)),n.zero2 = sum(y_dam==0|is.na(y_dam)) )

fecSummary2 = fec.noObs %>%
  dplyr::mutate(trueNA = ifelse(n.obs==((n.na1+n.na2)/2),1,0)) %>% 
  dplyr::mutate(obsZero = ifelse(n.obs==((n.zero1+n.zero2)/2),1,0)) %>%
  dplyr::select(site,year,n.obs,trueNA,obsZero)

fecSummary <- fecSummary1 %>%
  dplyr::bind_rows(fecSummary2)

siteNames=unique(fecSummary$site)

plot(NA,NA,xlim=c(2006,2020),ylim=c(0,20),
     axes=FALSE,frame=FALSE,xaxt='n',yaxt='n',
     xlab='',ylab='')
y.pt=20:1
for(i in 20:1){
  tmp=fecSummary[fecSummary$site==siteNames[i],]
  points(tmp$year,y=rep(y.pt[i],length(tmp$year)),
         pch=ifelse(tmp$trueNA==0,19,4),
         col=ifelse(tmp$obsZero==1,"gray","black"))
}
axis(1,  2006:2020, col.ticks = 1,las=2)
axis(2, (1:20),
     labels = rev(siteNames), las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Seeds per fruit
# -------------------------------------------------------------------
# -------------------------------------------------------------------

countSeedPerFruit <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countSeedPerFruit.RDS")

countSeedPerUndamagedFruit <- countSeedPerFruit %>%
  dplyr::filter(demography==1) %>%
  dplyr::filter(damaged==0) %>%
  dplyr::select(site,year,sdno)

countSeedPerDamagedFruit <- countSeedPerFruit %>%
  dplyr::filter(demography==1) %>%
  dplyr::filter(damaged==1) %>%
  dplyr::rename(sdno_dam = sdno) %>%
  dplyr::select(site,year,sdno_dam)

seeds.df = countSeedPerUndamagedFruit

seeds.noObs = seeds.df %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n.obs = n(),
                   n.na = sum(is.na(sdno)),
                   n.zero = sum(sdno==0|is.na(sdno)))

seedsSummary1 = seeds.noObs %>%
  dplyr::mutate(trueNA = ifelse(n.obs==n.na,1,0)) %>%
  dplyr::mutate(obsZero = ifelse(n.obs==n.zero,1,0)) %>%
  dplyr::select(site,year,n.obs,trueNA,obsZero)

countSeedPerDamagedFruit <- countSeedPerFruit %>%
  dplyr::filter(demography==1) %>%
  dplyr::filter(damaged==1) %>%
  dplyr::rename(sdno_dam = sdno) %>%
  dplyr::select(site,year,sdno_dam)

seeds.df = countSeedPerDamagedFruit

seeds.noObs = seeds.df %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n.obs = n(),
                   n.na = sum(is.na(sdno_dam)),
                   n.zero = sum(sdno_dam==0|is.na(sdno_dam)))

seedsSummary2 = seeds.noObs %>%
  dplyr::mutate(trueNA = ifelse(n.obs==n.na,1,0)) %>%
  dplyr::mutate(obsZero = ifelse(n.obs==n.zero,1,0)) %>%
  dplyr::select(site,year,n.obs,trueNA,obsZero)

# join undamaged and damaged datasets
seedsSummary <- seedsSummary2 %>%
  dplyr::left_join(seedsSummary1,by=c('site','year')) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(trueNA=ifelse((trueNA.x+trueNA.y)==2,1,0),
               obsZero=ifelse((obsZero.x+obsZero.y)==2,1,0),
               n.obs = n.obs.x+n.obs.y) %>% 
  dplyr::select(site,year,n.obs,trueNA,obsZero) %>%
  dplyr::bind_rows(seedsSummary1 %>% dplyr::filter(year<2013))

siteNames=unique(seedsSummary1$site)

plot(NA,NA,xlim=c(2006,2020),ylim=c(0,20),
     axes=FALSE,frame=FALSE,xaxt='n',yaxt='n',
     xlab='',ylab='')
y.pt=20:1
for(i in 20:1){
  tmp=seedsSummary[seedsSummary$site==siteNames[i],]
  points(tmp$year,y=rep(y.pt[i],length(tmp$year)),
         pch=ifelse(tmp$trueNA==0,19,4),
         col=ifelse(tmp$obsZero==1,"gray","black"))
}
axis(1,  2006:2020, col.ticks = 1,las=2)
axis(2, (1:20),
     labels = rev(siteNames), las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)

# combine plots

pdf("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/analysis/zero-fitness.pdf",width=8,height=6)
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(2006,2020),ylim=c(0,22),
     axes=FALSE,frame=FALSE,xaxt='n',yaxt='n',
     xlab='',ylab='')
polygon(x=c(2006,2007,2007,2006)+.5,
        y=c(0,0,21,21),
        col='gray97',border='gray97')
polygon(x=c(2008,2009,2009,2008)+.5,
        y=c(0,0,21,21),
        col='gray97',border='gray97')
polygon(x=c(2010,2011,2011,2010)+.5,
        y=c(0,0,21,21),
        col='gray97',border='gray97')
polygon(x=c(2012,2013,2013,2012)+.5,
        y=c(0,0,21,21),
        col='gray97',border='gray97')
polygon(x=c(2014,2015,2015,2014)+.5,
        y=c(0,0,21,21),
        col='gray97',border='gray97')
polygon(x=c(2016,2017,2017,2016)+.5,
        y=c(0,0,21,21),
        col='gray97',border='gray97')
polygon(x=c(2018,2019,2019,2018)+.5,
        y=c(0,0,21,21),
        col='gray97',border='gray97')

y.pt=20:1
for(i in 20:1){
  tmp=sigmaSummary[sigmaSummary$site==siteNames[i],]
  points(tmp$year-.25,y=rep(y.pt[i],length(tmp$year)),
         pch=ifelse(tmp$trueNA==0,19,4),
         col=ifelse(tmp$obsZero==1,"orange","black"))
  
  tmp=fecSummary[fecSummary$site==siteNames[i],]
  points(tmp$year,y=rep(y.pt[i],length(tmp$year)),
         pch=ifelse(tmp$trueNA==0,19,4),
         col=ifelse(tmp$obsZero==1,"orange","black"))
  
  tmp=seedsSummary[seedsSummary$site==siteNames[i],]
  points(tmp$year+.25,y=rep(y.pt[i],length(tmp$year)),
         pch=ifelse(tmp$trueNA==0,19,4),
         col=ifelse(tmp$obsZero==1,"orange","black"))
}
axis(1,  2006:2020, col.ticks = 1,las=2)
axis(2, (1:20),
     labels = rev(siteNames), las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)

legend("topleft",
       legend = c("Observations of fitness component",
                  "0 seedlings survive in plots",
                  "No observations"),
       pch=c(19,19,4),
       col=c("black","orange","orange"),
       cex=.75,
       horiz=TRUE,
       bty='n',
      x.intersp=0.5,
      text.width=c(2,5,4.5))

text(x=2009,y=21,cex=.5,
     "Seedling survival")
text(x=2011,y=21,cex=.5,
     "Fruits per plant")
text(x=2013,y=21,cex=.5,
     "Seeds per fruit")
segments(x0=2009.85,x1=2011-.25,y0=21,y1=20)
segments(x0=2011,y0=20.75,y1=20)
segments(x0=2011.25,x1=2012.3,y0=20,y1=21)

dev.off()

missing.list = list()
for(i in 1:20){
  tmp.sigma=sigmaSummary[sigmaSummary$site==siteNames[i],][1:13,]
  tmp.fec=fecSummary[fecSummary$site==siteNames[i],]
  tmp.seeds=seedsSummary[seedsSummary$site==siteNames[i],]
  tmp.seeds= tmp.seeds[ order( tmp.seeds$year),]
  na.obs=tmp.sigma$trueNA==1&tmp.fec$trueNA==1&tmp.seeds$trueNA==1
  na.obs2=tmp.sigma$obsZero==1&tmp.fec$trueNA==1&tmp.seeds$trueNA==1
  tmp=tmp.sigma[na.obs|na.obs2,1:2]
  missing.list[[i]] = tmp
}
lowFitnessYears=do.call(rbind,missing.list)
saveRDS(lowFitnessYears,"/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/output/lowFitnessYears.RDS")


missing.list = list()
for(i in 1:20){
  tmp.sigma=sigmaSummary[sigmaSummary$site==siteNames[i],][1:13,]
  tmp.fec=fecSummary[fecSummary$site==siteNames[i],]
  tmp.seeds=seedsSummary[seedsSummary$site==siteNames[i],]
  tmp.seeds= tmp.seeds[ order( tmp.seeds$year),]
  na.obs=tmp.sigma$trueNA==1
  na.obs2=tmp.sigma$obsZero==1
  tmp=tmp.sigma[na.obs|na.obs2,1:2]
  missing.list[[i]] = tmp
}
lowFitnessYearsPlots=do.call(rbind,missing.list)
saveRDS(lowFitnessYearsPlots,"/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/output/lowFitnessYearsPlots.RDS")
