#################################################################################
################################################################################
################################################################################
# Code for figures to compare the following modeling approaches for the seedling survivorship data:
# Filtering out data or recoding it
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

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/parallel/"
simFiles <- paste0(directory,list.files(directory))

mcmcSamplesFilter <- readRDS(simFiles[[grep("seedlingSurvivalSamples-filter",simFiles)]])
mcmcSamplesRecode <- readRDS(simFiles[[grep("seedlingSurvivalSamples-recode",simFiles)]])

################################################################################
# Data directory
#################################################################################

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/parallel/"
dataFiles <- paste0(directory,list.files(directory))

dataFilter <- readRDS(dataFiles[[grep("Data-filter",dataFiles)]])
dataRecode <- readRDS(dataFiles[[grep("Data-recode",dataFiles)]])

# -------------------------------------------------------------------
# Get site names and position
# -------------------------------------------------------------------
directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedlingSurvival/"
dataFiles <- paste0(directory,list.files(directory))

censusSeedlingsFruitingPlants <- readRDS(dataFiles[[grep("censusSeedlings",dataFiles)]])

siteNames = unique(censusSeedlingsFruitingPlants$site)
years = as.numeric(unique(censusSeedlingsFruitingPlants$year))
yearNames = years[order(years)]

position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,dominant.surface.rock.type) %>%
  dplyr::mutate(easting=easting/1000)

################################################################################
# Calculate summary objects
#################################################################################


summary.fun = function(x){
  med = median(x)
  hpdi.interval = HPDI(x, prob = .89)
  hpdi.interval2 = HPDI(x, prob = .5)
  return(c(med,hpdi.interval,hpdi.interval2))
}

f=function(x){
  
  parm.mu0=MCMCchains(x,params='mu0')
  parm.prob0=apply(parm.mu0,2,boot::inv.logit)
  parm.prob0.sum = apply(parm.prob0,2,summary.fun)
  
  df.list = list()
  for(i in 1:20){
    obj = parm.prob0.sum
    index=grep(paste0("\\[",i,"\\]"),colnames(obj))
    tmp = signif(parm.prob0.sum[,index],3)
    tmp.df=data.frame(site=siteNames[i],(t(tmp)))
    names(tmp.df) = c("site","med","ci.lo","ci.hi","ci.lo2","ci.hi2")
    rownames(tmp.df) = NULL
    df.list[[i]] = tmp.df
  }
  
  summary.pop.df=do.call(rbind,df.list)
  return(summary.pop.df)
}

parm.mu0.filter = f(mcmcSamplesFilter)
parm.mu0.recode = f(mcmcSamplesRecode)


f2 = function(x){
  parm.mu=MCMCchains(x,params='mu')
  
  parm.prob=apply(parm.mu,2,boot::inv.logit)
  parm.prob.sum = apply(parm.prob,2,summary.fun)
  
  df.list = list()
  for(i in 1:20){
    obj = parm.prob.sum
    index=grep(paste0("\\[",i,","),colnames(obj))
    tmp = signif(parm.prob.sum[,index],3)
    tmp.df=data.frame(site=siteNames[i],year=yearNames,(t(tmp)))
    names(tmp.df) = c("site","year","med","ci.lo","ci.hi","ci.lo2","ci.hi2")
    rownames(tmp.df) = NULL
    df.list[[i]] = tmp.df
  }
  
  summary.df=do.call(rbind,df.list)
  return(summary.df)
}

parm.mu.filter <- f2(mcmcSamplesFilter)
parm.mu.recode <- f2(mcmcSamplesRecode)

# rename
sigma_p = summary.pop.df %>% dplyr::select(site,med)
sigma_py = summary.df 

################################################################################
# Make summary plots
#################################################################################


pdf("~/Dropbox/clarkiaSeedBanks/products/figures/modelChecks/compareSurvivalModels.pdf",width=6,height=6)
par(mfrow = c(1,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

time.sample = 1:14+2005
for(i in 1:20){
  
  # filtering out data with errors
  sigma_py = parm.mu.filter
  sigma_p = parm.mu0.filter
  data = dataFilter
  
  tmp=sigma_py[sigma_py$site==siteNames[i],]
  tmp.pop=sigma_p[sigma_p$site==siteNames[i],]
  
  tmp.vec=c()
  for(j in 1:14){
    tmp.vec[j]=sum(data$seedlingNumber[data$site==i&data$year==j])
  }
  
  plot(NA,NA,
       ylim=c(0,1),pch=16,xlim=c(2006,2019),
       ylab='',xlab='',xaxt='n',yaxt='n')
  
  index=tmp.vec==0
  if (sum(index)>0) {for(j in 1:length(time.sample[index])){
    ts=time.sample[index]
    polygon(x=c(ts[j]-.5,ts[j]-.5,
                ts[j]+.5, ts[j]+.5),
            y=c(-.1,1.1,1.1,-.1),col='gray95',border='gray95')
  }} else {NA}
  
  
  polygon(x=c(2005,2020,2020,2005),
          y=c(tmp.pop$ci.lo,tmp.pop$ci.lo,tmp.pop$ci.hi,tmp.pop$ci.hi),
          col='gray95',border='gray95')
  
  # recoding data with errors
  sigma_py.recode = parm.mu.recode
  sigma_p.recode = parm.mu0.recode
  data.recode = dataRecode
  
  tmp.recode=sigma_py.recode[sigma_py.recode$site==siteNames[i],]
  tmp.pop.recode=sigma_p.recode[sigma_p.recode$site==siteNames[i],]

  polygon(x=c(2005,2020,2020,2005),
          y=c(tmp.pop.recode$ci.lo,tmp.pop.recode$ci.lo,
              tmp.pop.recode$ci.hi,tmp.pop.recode$ci.hi),
          col=rgb(1,.647,0,.25),border=rgb(1,.647,0,.25))
  
  abline(h=tmp.pop.recode$med,col=rgb(1,.647,0,.25))
  
  abline(h=tmp.pop$med,col='gray')
  segments(time.sample-.15,y0=tmp$ci.lo,y1=tmp$ci.hi)
  points(time.sample-.15,
         tmp$med,pch=21,cex=1,
         col="black",bg='white')
  
  
  segments(time.sample+.15,y0=tmp.recode$ci.lo,y1=tmp.recode$ci.hi,col='orange')
  points(time.sample+.15,
         tmp.recode$med,pch=21,cex=1,
         col="orange",bg='white')
  
  # label
  text(x=2005.5,y=1,siteNames[i],pos=4)
  axis(1L);axis(2L)
  legend(x = 2005.5, y = 1,
         col = c('gray','orange'),
         lty = c(1,1),
         legend = c("Filter out undercounts","Recode undercounts"),
         cex=.55,
         box.lty=0)
#  ifelse(i%in%c(16:20),axis(1L),NA)
#  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
  # ifelse(i%in%c(5), , NA)
  mtext("Year", side = 1, outer = TRUE, line = 2.2)
  mtext("Probability of seedling survival to fruiting", side = 2, outer = TRUE, line = 2.2)
  
}
#mtext("Population*year-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)
dev.off()

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
