#################################################################################
################################################################################
################################################################################
# Code for figures of total fruit equivalents per plant
# 
# Scripts by Gregor Siegmund
# fausto.siegmund@gmail.com
# last updated 03-07-2021
#################################################################################
#################################################################################
#################################################################################
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)

library(MCMCvis)
library(tidybayes)
library(tidyverse)
library(magrittr)
library(bayesplot)
library(rethinking)

directory = "/Users/Gregor/Dropbox/dataLibrary/posteriors/"
simFiles <- paste0(directory,list.files(directory))

mcmcSamples <- readRDS(simFiles[[grep("fitnessSamples.rds",simFiles)]])

dirPars = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/parameterSummary/"
dir.create(file.path(dirPars), showWarnings = FALSE)

################################################################################
# Data directory
#################################################################################

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/fitness-weakly-informative/"
dataFiles <- paste0(directory,list.files(directory))

data <- readRDS(dataFiles[[grep("data.rds",dataFiles)]])


################################################################################
# Create composite
#################################################################################

countFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countFruitsPerPlantAllPlots.RDS")

siteNames = unique(countFruitsPerPlantAllPlots$site)

# recover site indices and year indices
siteIndex <- data.frame(siteIndex=unique(countFruitsPerPlantAllPlots$site),site=unique(data$site))
yearIndex <- data.frame(yearIndex=unique(countFruitsPerPlantAllPlots$year),
                        year=unique(data$year)) 

summary.fun = function(x){
  quant=quantile(x,c(.025,.975))
  mu = mean(x)
  med = median(x)
  hpdi.interval = HPDI(x, prob = .95)
  return(c(mu,quant,med,hpdi.interval))
}

parm.mu_tfe=MCMCchains(mcmcSamples, params=c("mu_tfe"))
parm.mu_tfe.sum = apply(parm.mu_tfe,2,summary.fun)

df.list = list()
for(i in 1:20){
  obj = parm.mu_tfe.sum
  index=grep(paste0("\\[",i,","),colnames(obj))
  tmp = signif(obj[,index],3)
  tmp.df=data.frame(site=(siteIndex[,1])[i],year=yearIndex[,1],(t(tmp)))
  names(tmp.df) = c("site","year","mu","ci.lo95","ci.hi95","med","hpdi.lo95","hpdi.hi95")
  rownames(tmp.df) = NULL
  df.list[[i]] = tmp.df
}

summary.mu_tfe.df=do.call(rbind,df.list)


# tfeDF<-mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(mu_py_tfe[site,year]) %>%
#   dplyr::group_by(site,year) %>%
#   dplyr::summarise(mu = mean(mu_py_tfe),
#                    ci.lo95 = quantile(mu_py_tfe,probs=0.025), 
#                    ci.hi95 = quantile(mu_py_tfe,probs=0.975),
#                    med = median(mu_py_tfe), 
#                    hpdi.lo95 = HPDI(mu_py_tfe,.95)[1], 
#                    hpdi.hi95 = HPDI(mu_py_tfe,.95)[2]
#   ) %>%
#   dplyr::ungroup() %>%
#   dplyr::left_join(siteIndex,by="site") %>%
#   dplyr::left_join(yearIndex,by="year") %>%
#   dplyr::select(-c(site,year)) %>%
#   dplyr::rename(site = siteIndex) %>%
#   dplyr::rename(year = yearIndex)


directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/"
dataFiles <- paste0(directory,list.files(directory))

countSeedPerDamagedFruit <- readRDS(dataFiles[[grep("countSeedPerDamagedFruit",dataFiles)]])

siteNames = unique(countSeedPerDamagedFruit$site4) 

# recover site indices and year indices
siteIndex <- data.frame(siteIndex=unique(countSeedPerDamagedFruit$site4),site=unique(data$site4))
yearIndex <- data.frame(yearIndex=unique(countSeedPerDamagedFruit$year4),
                        year=unique(data$year4)) 

# tfeCompDF<- mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(mu_py_tfe_comp[site,year]) %>%
#   dplyr::group_by(site,year) %>%
#   dplyr::summarise(mu = mean(mu_py_tfe_comp),
#                    ci.lo95 = quantile(mu_py_tfe_comp,probs=0.025), 
#                    ci.hi95 = quantile(mu_py_tfe_comp,probs=0.975),
#                    med = median(mu_py_tfe_comp), 
#                    hpdi.lo95 = HPDI(mu_py_tfe_comp,.95)[1], 
#                    hpdi.hi95 = HPDI(mu_py_tfe_comp,.95)[2]
#   ) %>%
#   dplyr::ungroup() %>%
#   dplyr::left_join(siteIndex,by="site") %>%
#   dplyr::left_join(yearIndex,by="year") %>%
#   dplyr::select(-c(site,year)) %>%
#   dplyr::rename(site = siteIndex) %>%
#   dplyr::rename(year = yearIndex) 


parm.mu_und=MCMCchains(mcmcSamples, params=c("mu_und"))
parm.mu_und.sum = apply(parm.mu_und,2,summary.fun)

df.list = list()
for(i in 1:20){
  obj = parm.mu_und.sum
  index=grep(paste0("\\[",i,","),colnames(obj))
  tmp = signif(obj[,index],3)
  tmp.df=data.frame(site=(siteIndex[,1])[i],year=yearIndex[,1],(t(tmp)))
  names(tmp.df) = c("site","year","mu","ci.lo95","ci.hi95","med","hpdi.lo95","hpdi.hi95")
  rownames(tmp.df) = NULL
  df.list[[i]] = tmp.df
}

summary.parm.mu_und.df=do.call(rbind,df.list)


parm.mu_dam=MCMCchains(mcmcSamples, params=c("mu_dam"))
parm.mu_dam.sum = apply(parm.mu_dam,2,summary.fun)

df.list = list()
for(i in 1:20){
  obj = parm.mu_dam.sum
  index=grep(paste0("\\[",i,","),colnames(obj))
  tmp = signif(obj[,index],3)
  tmp.df=data.frame(site=(siteIndex[,1])[i],year=yearIndex[,1],(t(tmp)))
  names(tmp.df) = c("site","year","mu","ci.lo95","ci.hi95","med","hpdi.lo95","hpdi.hi95")
  rownames(tmp.df) = NULL
  df.list[[i]] = tmp.df
}

summary.parm.mu_dam.df=do.call(rbind,df.list)

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/"
dataFiles <- paste0(directory,list.files(directory))

countSeedPerUndamagedFruit <- readRDS(dataFiles[[grep("countSeedPerUndamagedFruit",dataFiles)]])

siteNames = unique(countSeedPerDamagedFruit$site3) 

# recover site indices and year indices
siteIndex <- data.frame(siteIndex=unique(countSeedPerUndamagedFruit$site3),site=unique(data$site3))
yearIndex <- data.frame(yearIndex=unique(countSeedPerUndamagedFruit$year3),
                        year=unique(data$year3)) 

parm.mu_seeds=MCMCchains(mcmcSamples, params=c("mu_seeds"))
parm.mu_seeds.sum = apply(parm.mu_seeds,2,summary.fun)

df.list = list()
for(i in 1:20){
  obj = parm.mu_seeds.sum
  index=grep(paste0("\\[",i,","),colnames(obj))
  tmp = signif(obj[,index],3)
  tmp.df=data.frame(site=(siteIndex[,1])[i],year=yearIndex[,1],(t(tmp)))
  names(tmp.df) = c("site","year","mu","ci.lo95","ci.hi95","med","hpdi.lo95","hpdi.hi95")
  rownames(tmp.df) = NULL
  df.list[[i]] = tmp.df
}

summary.parm.mu_seeds.df=do.call(rbind,df.list)


directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/"
dataFiles <- paste0(directory,list.files(directory))

countSeedPerDamagedFruit <- readRDS(dataFiles[[grep("countSeedPerDamagedFruit",dataFiles)]])

siteNames = unique(countSeedPerDamagedFruit$site4) 

# recover site indices and year indices
siteIndex <- data.frame(siteIndex=unique(countSeedPerDamagedFruit$site4),site=unique(data$site4))
yearIndex <- data.frame(yearIndex=unique(countSeedPerDamagedFruit$year4),
                        year=unique(data$year4)) 

parm.mu_dam_seeds=MCMCchains(mcmcSamples, params=c("mu_dam_seeds"))
parm.mu_dam_seeds.sum = apply(parm.mu_dam_seeds,2,summary.fun)

df.list = list()
for(i in 1:20){
  obj = parm.mu_dam_seeds.sum
  index=grep(paste0("\\[",i,","),colnames(obj))
  tmp = signif(obj[,index],3)
  tmp.df=data.frame(site=(siteIndex[,1])[i],year=yearIndex[,1],(t(tmp)))
  names(tmp.df) = c("site","year","mu","ci.lo95","ci.hi95","med","hpdi.lo95","hpdi.hi95")
  rownames(tmp.df) = NULL
  df.list[[i]] = tmp.df
}

summary.parm.mu_dam_seeds.df=do.call(rbind,df.list)

# calculate fraction of seeds in a damaged verus undamaged fruit
unique(summary.parm.mu_dam_seeds.df$year)

# recover site indices and year indices
countSeedPerDamagedFruit <- readRDS(dataFiles[[grep("countSeedPerDamagedFruit",dataFiles)]])
siteIndexDam <- data.frame(siteIndex=unique(countSeedPerDamagedFruit$site4),site=unique(data$site4))
yearIndexDam <- data.frame(yearIndex=unique(countSeedPerDamagedFruit$year4),
                        year=unique(data$year4)) 

countSeedPerUndamagedFruit <- readRDS(dataFiles[[grep("countSeedPerUndamagedFruit",dataFiles)]])
siteIndexUnd <- data.frame(siteIndex=unique(countSeedPerUndamagedFruit$site3),site=unique(data$site3))
yearIndexUnd <- data.frame(yearIndex=unique(countSeedPerUndamagedFruit$year3),
                        year=unique(data$year3)) 

ratioIndex=yearIndexUnd %>%
  dplyr::left_join(yearIndexDam,by='yearIndex')

parm.mu_dam_seeds=MCMCchains(mcmcSamples, params=c("mu_dam_seeds"))
mu_seeds.lastyears=parm.mu_seeds[,141:260]

ratio.seeds=parm.mu_dam_seeds/mu_seeds.lastyears

parm.mu_comp=parm.mu_und + parm.mu_dam*ratio.seeds


directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/"
dataFiles <- paste0(directory,list.files(directory))

countSeedPerDamagedFruit <- readRDS(dataFiles[[grep("countSeedPerDamagedFruit",dataFiles)]])

# recover site indices and year indices
siteIndex <- data.frame(siteIndex=unique(countSeedPerDamagedFruit$site4),site=unique(data$site4))
yearIndex <- data.frame(yearIndex=unique(countSeedPerDamagedFruit$year4),
                        year=unique(data$year4)) 

parm.mu_comp.sum = apply(parm.mu_comp,2,summary.fun)

df.list = list()
for(i in 1:20){
  obj = parm.mu_comp.sum
  index=grep(paste0("\\[",i,","),colnames(obj))
  tmp = signif(obj[,index],3)
  tmp.df=data.frame(site=(siteIndex[,1])[i],year=yearIndex[,1],(t(tmp)))
  names(tmp.df) = c("site","year","mu","ci.lo95","ci.hi95","med","hpdi.lo95","hpdi.hi95")
  rownames(tmp.df) = NULL
  df.list[[i]] = tmp.df
}

summary.parm.mu_comp.df=do.call(rbind,df.list)

################################################################################
# Summarize
#################################################################################


# manually filter out S22 2013
# tfeCompDF<- tfeCompDF %>% 
#   dplyr::filter(site!="S22"|year!=2013) %>%
#   dplyr::filter(site!="OSR"|year!=2013) %>%
#   dplyr::filter(site!="LCE"|year!=2013) %>%
#   dplyr::filter(site!="GCN"|year!=2013)

reference<-summary.parm.mu_comp.df %>% dplyr::select(site,year) %>% unique()
reference$year <- as.integer(reference$year)
summary.parm.mu_comp.df$year<-as.integer(summary.parm.mu_comp.df$year)

# undamaged
# countFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/countUndamagedDamagedFruitsPerPlantAllPlots.rds")
# siteIndex <- data.frame(siteIndex=unique(countFruitsPerPlantAllPlots$site),site=unique(data$site))
# yearIndex <- data.frame(yearIndex=unique(countFruitsPerPlantAllPlots$year),
#                         year=unique(data$year)) 
# 
# ufDF<-mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(mu_py_und[site,year]) %>%
#   dplyr::group_by(site,year) %>%
#   dplyr::summarise(mu = mean(mu_py_und),
#                    ci.lo95 = quantile(mu_py_und,probs=0.025), 
#                    ci.hi95 = quantile(mu_py_und,probs=0.975),
#                    med = median(mu_py_und), 
#                    hpdi.lo95 = HPDI(mu_py_und,.95)[1], 
#                    hpdi.hi95 = HPDI(mu_py_und,.95)[2]
#   ) %>%
#   dplyr::ungroup() %>%
#   dplyr::left_join(siteIndex,by="site") %>%
#   dplyr::left_join(yearIndex,by="year") %>%
#   dplyr::select(-c(site,year)) %>%
#   dplyr::rename(site = siteIndex) %>%
#   dplyr::rename(year = yearIndex)

summary.parm.mu_und.df$year<-as.integer(summary.parm.mu_und.df$year)



referenceUF<-summary.parm.mu_und.df %>% dplyr::select(site,year) %>% unique()
# 
# ufOnly<-setdiff(referenceUF,reference) %>%
#   dplyr::mutate(site.year=paste0(site,year))
# 
# ufDF_missing<-ufDF %>%
#   dplyr::mutate(site.year=paste0(site,year)) %>%
#   dplyr::filter(site.year %in%ufOnly$site.year) %>%
#   dplyr::select(-site.year)

# bind all datasets
tfeFullDF<-bind_rows(summary.mu_tfe.df,summary.parm.mu_comp.df)


# fill in with NAs
referenceDF <- expand.grid(site=unique(tfeFullDF$site),year=unique(tfeFullDF$year))

tfeFullDF<-tfeFullDF %>%
  dplyr::full_join(referenceDF,by=c("site","year"))

write.csv(tfeFullDF,file=paste0(dirPars,"fruitEquivalentsPerPlantSummary.csv"))

################################################################################
# Make summary plots
#################################################################################

summary.fun = function(x){
  med = median(x)
  hpdi.interval = HPDI(x, prob = .89)
  hpdi.interval2 = HPDI(x, prob = .5)
  return(c(med,hpdi.interval,hpdi.interval2))
}

countFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countFruitsPerPlantAllPlots.RDS")
# recover site indices and year indices
siteIndex <- data.frame(siteIndex=unique(countFruitsPerPlantAllPlots$site),site=unique(data$site))
yearIndex <- data.frame(yearIndex=unique(countFruitsPerPlantAllPlots$year),
                        year=unique(data$year)) 

parm.mu_tfe.sum = apply(parm.mu_tfe,2,summary.fun)

df.list = list()
for(i in 1:20){
  obj = parm.mu_tfe.sum
  index=grep(paste0("\\[",i,","),colnames(obj))
  tmp = signif(obj[,index],3)
  tmp.df=data.frame(site=(siteIndex[,1])[i],year=yearIndex[,1],(t(tmp)))
  names(tmp.df) = c("site","year","med","ci.lo","ci.hi","ci.lo2","ci.hi2")
  rownames(tmp.df) = NULL
  df.list[[i]] = tmp.df
}

summary.mu_tfe.df=do.call(rbind,df.list)

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/"
dataFiles <- paste0(directory,list.files(directory))

countSeedPerDamagedFruit <- readRDS(dataFiles[[grep("countSeedPerDamagedFruit",dataFiles)]])

# recover site indices and year indices
siteIndex <- data.frame(siteIndex=unique(countSeedPerDamagedFruit$site4),site=unique(data$site4))
yearIndex <- data.frame(yearIndex=unique(countSeedPerDamagedFruit$year4),
                        year=unique(data$year4)) 

parm.mu_comp.sum = apply(parm.mu_comp,2,summary.fun)

df.list = list()
for(i in 1:20){
  obj = parm.mu_comp.sum
  index=grep(paste0("\\[",i,","),colnames(obj))
  tmp = signif(obj[,index],3)
  tmp.df=data.frame(site=(siteIndex[,1])[i],year=as.numeric(yearIndex[,1]),(t(tmp)))
  names(tmp.df) = c("site","year","med","ci.lo","ci.hi","ci.lo2","ci.hi2")
  rownames(tmp.df) = NULL
  df.list[[i]] = tmp.df
}

summary.parm.mu_comp.df=do.call(rbind,df.list)


# bind all datasets
tfeFullDF<-bind_rows(summary.mu_tfe.df,summary.parm.mu_comp.df)

# fill in with NAs
referenceDF <- expand.grid(site=unique(tfeFullDF$site),year=unique(tfeFullDF$year))

tfeFullDF<-tfeFullDF %>%
  dplyr::full_join(referenceDF,by=c("site","year"))


# time plot
ribbonTFEfull <-   tfeFullDF %>%
  ggplot(aes(x=year,y=med,group=site,ymin=ci.lo,ymax=ci.hi)) +
  geom_ribbon(alpha=0.5,colour="grey50") +
  geom_line() +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  xlab("Year") +
  ylab("Total fruit equivalents per plant") +
  facet_wrap(~site,scales='free_y') +
  labs(title="Total fruit equivalents per plant, composite (2006-2018)") 


# time plot
tfeFullDFsorted <-  tfeFullDF %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(desc(med),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  mutate(id = row_number())

interannualTFEfull <-   ggplot(tfeFullDFsorted,aes(x = id , y = med)) + 
  geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +  
  geom_point(aes(color=year),size=1) +
  coord_flip() +
  facet_grid(site ~ ., scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ylab("Fruits per plant (total fruit equivalents)")  +
  labs(color="Year") 


position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,dominant.surface.rock.type) %>%
  dplyr::mutate(easting=easting/1000)

tfeFullDFsorted <-  tfeFullDFsorted %>%
  dplyr::left_join(position,by="site")

spatialTFEfull <-   ggplot(tfeFullDFsorted,aes(x = easting , y = med)) + 
  geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +  
  geom_point(size=1) +
  theme_bw() +
  facet_wrap( ~ year, scales="free_y",nrow=2) +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  xlab("Easting") +
  ylab("Fruits per plant (total fruit equivalents)")  +
  labs(color="Year")

dirFigures = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/"

ggsave(filename=paste0(dirFigures,"interannualTFEfull.pdf"),
       plot=interannualTFEfull,width=6,height=12)

ggsave(filename=paste0(dirFigures,"spatialTFEfull.pdf"),
       plot=spatialTFEfull,width=12,height=8)


################################################################################
# Make summary plots
#################################################################################

summary_py = tfeFullDF
#sigma_p = summary.pop.df

pdf("~/Dropbox/clarkiaSeedBanks/products/figures/TFE-population.pdf",width=8,height=6)
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1    )

time.sample = 1:13+2005
for(i in 1:20){
  
  tmp=summary_py[summary_py$site==siteNames[i],]
 # tmp.pop=sigma_p[sigma_p$site==siteNames[i],]
  
  tmp.vec=c()
  # need 2 for-loops because datasets indexed separately
  for(j in 1:7){
    tmp.vec[j]=sum(data$y_tfe[data$site==i&data$year==j])
  }
  
  # second loop indexed to start at 8
  for(j in 1:6){
    tmp.vec[7+j]=sum(data$y_und[data$site2==i&data$year2==j],
                   data$y_dam[data$site2==i&data$year2==j])
    
  }
  
  tmp.vec=ifelse(is.na(tmp.vec),0,tmp.vec)


  
  plot(NA,NA,
       ylim=c(0,max(tmp[,3:7])),pch=16,xlim=c(2006,2019),
       ylab='',xlab='',xaxt='n',yaxt='n')
  
  index=tmp.vec==0
  if (sum(index)>0) {for(j in 1:length(time.sample[index])){
    ts=time.sample[index]
    polygon(x=c(ts[j]-.5,ts[j]-.5,
                ts[j]+.5, ts[j]+.5),
            y=c(-.1,max(tmp[,3:7])*2,max(tmp[,3:7])*2,-.1),col='gray95',border='gray95')
  }} else {NA}
  
  
  # polygon(x=c(2005,2020,2020,2005),
  #         y=c(tmp.pop$ci.lo,tmp.pop$ci.lo,tmp.pop$ci.hi,tmp.pop$ci.hi),
  #         col='gray95',border='gray95')
  
  #abline(h=tmp.pop$med,col='gray')
  segments(time.sample,y0=tmp$ci.lo,y1=tmp$ci.hi)
  
  points(time.sample,
         tmp$med,pch=21,cex=1,
         col="black",bg='white')
  
  
  text(x=2005.5,y=.9*max(tmp[,3:7]),siteNames[i],pos=4)
  ifelse(i%in%c(16:20),axis(1L),NA)
  axis(2, seq(0,max(tmp[,3:7]),by=5), tick=FALSE,
       labels = seq(0,max(tmp[,3:7]),by=5), las = 1, 
       cex.axis = 1,line=0,mgp=c(3,.25,0))
  #ifelse(i%in%c(1,6,11,16),,NA)
  ifelse(i%in%c(5), legend(x = 15, y = 1,
                           col = c('gray','orange'),
                           lty = c(1,1),
                           legend = c("Persistence only","Persistence & viability"),
                           cex=.55,
                           box.lty=0), NA)
}
mtext("Year", side = 1, outer = TRUE, line = 2.2)
mtext("Total fruit equivalents per plant", side = 2, outer = TRUE, line = 2.2)
#mtext("Population*year-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)
dev.off()

# 
# sigma0=MCMCchains(mcmcSamples,params="sigma0")
# 
# x.sum=apply(sigma0,2,quantile,probs=c(0.025,.25,.5,.75,.975))
# index=order(x.sum[3,],decreasing=TRUE)
# x.sum=x.sum[,index]
# par(mfrow=c(1,1))
# plot(NA,NA,type='n',xlim=c(0,3),ylim=c(0,20),
#      axes=FALSE,frame=FALSE,
#      xlab="",ylab="")
# y.pt = 1:20
# for(i in 1:20){
#   tmp<-x.sum[,i]
#   segments(x0=tmp[1],x1=tmp[5],y0=y.pt[i])
#   segments(x0=tmp[2],x1=tmp[4],y0=y.pt[i],lwd=3)
#   points(x=tmp[3],y=y.pt[i],pch=21,bg='white')
#   
# }
# axis(1,  seq(0,3,by=.5), col.ticks = 1)
# axis(2, (1:20),
#      labels = siteNames[index], las = 1, 
#      col = NA, col.ticks = 1, cex.axis = 1)
# mtext("sigma0",
#       side=1,line=2.5,adj=.5,col='black',cex=1)
# 
# df.sum=censusSeedlingsFruitingPlants %>%
#   dplyr::group_by(site) %>%
#   dplyr::summarise(var(seedlingNumber,na.rm=TRUE))
# 
# sigma=MCMCchains(mcmcSamples,params="sigma")
# x.sum=apply(sigma,2,quantile,probs=c(0.025,.25,.5,.75,.975))
# 
# par(mfrow = c(4,5),
#     oma = c(5,4,0,0) + 0.1,
#     mar = c(0,0,1,1) + 0.1)
# 
# time.sample = 1:14+2005
# for(i in 1:20){
#   
#   index=grep(paste0("\\[",i,","),colnames(x.sum))
#   tmp=x.sum[,index]
#   
#   plot(NA,NA,
#        ylim=c(0,4),pch=16,xlim=c(2006,2019),
#        ylab='',xlab='',xaxt='n',yaxt='n')
#   
#   segments(time.sample,y0=tmp[1,],y1=tmp[5,])
#   points(time.sample,
#          tmp[3,],pch=21,cex=1,
#          col="black",bg='white')
#   
#   
#   text(x=2005.5,y=3.75,siteNames[i],pos=4)
#   ifelse(i%in%c(16:20),axis(1L),NA)
#   ifelse(i%in%c(1,6,11,16),axis(2L),NA)
#   ifelse(i%in%c(5), legend(x = 15, y = 1,
#                            col = c('gray','orange'),
#                            lty = c(1,1),
#                            legend = c("Persistence only","Persistence & viability"),
#                            cex=.55,
#                            box.lty=0), NA)
# }
# mtext("Year", side = 1, outer = TRUE, line = 2.2)
# mtext("sigma", side = 2, outer = TRUE, line = 2.2)
# #mtext("Population*year-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)
# 
# sigma0=MCMCchains(mcmcSamples,params="sigma0")
# 
