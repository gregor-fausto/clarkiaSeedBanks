# -------------------------------------------------------------------
# Analysis of fitness models
# Outputs reproductive success estimates
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE,max.print=100000)
# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(rjags) # jags interface
library(MCMCvis)
library(tidyverse)
library(reshape2)
library(HDInterval)
library(bayesplot)
library(coda)
library(rethinking)

# -------------------------------------------------------------------
# Seedling survival to fruiting
# -------------------------------------------------------------------
directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/parallel/"
sampleFiles <- paste0(directory,list.files(directory))

sigmaMCMCsamples <- readRDS(sampleFiles[grep("seedlingSurvivalSamples-recode-2.rds",sampleFiles)])
sigmaMCMCsamples <- MCMCchains(sigmaMCMCsamples,params="mu")
sampleInd = sample(1:45000,30000,replace=FALSE)
sigmaMCMCsamples <- sigmaMCMCsamples[sampleInd,]
# -------------------------------------------------------------------
# Fruits per plant
# -------------------------------------------------------------------
directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/"
sampleFiles <- paste0(directory,list.files(directory))

observedTfeMCMCsamples <- readRDS(sampleFiles[grep("observedTFE.RDS",sampleFiles)])
compositeTfeMCMCsamples <- readRDS(sampleFiles[grep("compositeTFE.RDS",sampleFiles)])

# -------------------------------------------------------------------
# Seeds per fruit
# -------------------------------------------------------------------
directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/"
sampleFiles <- paste0(directory,list.files(directory))

observedSeedsMCMCsamples <- readRDS(sampleFiles[grep("observedSeeds.RDS",sampleFiles)])

################################################################################
# Create composite
#################################################################################

# drop last year, change if working with full dataset
sigmaMCMCsamples = sigmaMCMCsamples[,c(1:260)]
tfeMCMCsamples = cbind(observedTfeMCMCsamples,compositeTfeMCMCsamples)

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/"
saveRDS(boot::inv.logit(sigmaMCMCsamples),file="~/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/sigma-analysis.RDS")
saveRDS(tfeMCMCsamples,file="~/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/tfe-analysis.RDS")
saveRDS(observedSeedsMCMCsamples,file="~/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/phi-analysis.RDS")

# check to make sure data frames are the same size
# this is not automatic at the moment
dim(sigmaMCMCsamples)==dim(tfeMCMCsamples)
dim(sigmaMCMCsamples)==dim(observedSeedsMCMCsamples)
dim(tfeMCMCsamples)==dim(observedSeedsMCMCsamples)

rsMCMCsamples = boot::inv.logit(sigmaMCMCsamples)*tfeMCMCsamples*observedSeedsMCMCsamples
saveRDS(rsMCMCsamples,file="~/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/rsPosterior.RDS")

seedsPerPlant = tfeMCMCsamples*observedSeedsMCMCsamples
saveRDS(seedsPerPlant,file="~/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/seedsPerPlantPosterior.RDS")

# -------------------------------------------------------------------
# Summarize
# -------------------------------------------------------------------

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/"
dataFiles <- paste0(directory,list.files(directory))

countSeedPerUndamagedFruit <- readRDS(dataFiles[[grep("countSeedPerUndamagedFruit",dataFiles)]])
siteIndex <- data.frame(siteIndex=unique(countSeedPerUndamagedFruit$site3),site=1:20)
yearIndex <- data.frame(yearIndex=unique(countSeedPerUndamagedFruit$year3),
                        year=1:13) 
siteNames=unique(siteIndex$siteIndex)

# calculate posterior mode; 
# directly from https://github.com/mikemeredith/BEST/blob/master/R/plotPost.R
posterior.mode = function(x){
  x.max=max(x)
  x.min=min(x)
  dres <- density( x ,from = x.min, to = x.max)
  modeParam <- dres$x[which.max(dres$y)]
  return(modeParam)
}

summary.fun = function(x){
  quant=quantile(x,c(.025,.975))
  mu = mean(x)
  med = median(x)
  mode = posterior.mode(x)
  hpdi.interval = HPDI(x, prob = .95)
  return(c(mu,med,mode,quant,hpdi.interval))
}

rsSummary.tmp = apply(rsMCMCsamples,2,summary.fun)

par(mfrow=c(1,3))
plot(x=rsSummary.tmp[1,],y=rsSummary.tmp[2,],xlab="Mean",ylab="Median")
plot(x=rsSummary.tmp[1,],y=rsSummary.tmp[3,],xlab="Mean",ylab="Mode")
plot(x=rsSummary.tmp[2,],y=rsSummary.tmp[3,],xlab="Median",ylab="Mode")

# -------------------------------------------------------------------
# Examine relationship between climate and per-capita reproductive success
# Sensitivity of RS to precipitation
# -------------------------------------------------------------------
climate=readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/climateData.RDS")
df <- climate %>% dplyr::filter(intenseDemography==1)
siteNames <- unique(df$site)
position <- unique(df$easting)
siteNames <- siteNames[order(position,decreasing=FALSE)]


climate=df %>%
  dplyr::filter(year>2005&year<2019)
climate.winter = climate[climate$season=="winter",]
climate.spring = climate[climate$season=="spring",]

#dev.off()
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,0,0) + 0.1)
ind=(1:20)[order(position)]
b = c()
for(i in 1:20){
  index=grep(paste0("\\[",ind[i],","),colnames(rsSummary.tmp))
  x=climate.spring[climate.spring$site==siteNames[ind[i]],]$p
  x=(x-mean(x))#/(2*sd(x))
  y=rsSummary.tmp[3,index]
  y=log(y)
  mod=lm(y~x)
  
  plot(x,y,xaxt='n',yaxt='n',pch=19,
       xlim=c(-120,120),ylim=c(-7,7))
  abline(v=0,lty='dotted')
  #abline(a=coef(mod)[1],b=coef(mod)[2])
  f=function(a,b,x){
    x.min=min(x);x.max=max(x);
    segments(y0=a+x.min*b,y1=a+x.max*b,x0=x.min,x1=x.max)
  }
  f(a=coef(mod)[1],b=coef(mod)[2],x=x)
  b[i]=coef(mod)[2]
  #text(.9*max(log(x+1)),.9*max(log(y+.5)),signif(cor(x,y),2),col='red')
  text(-125,-6,paste0("R2=",signif(summary(mod)$r.squared,2)),col='red',pos=4)
  text(-125,-4,paste0("b=",signif(coef(mod)[2],2)),col='red',pos=4)
  text(-125,-2,paste0(siteNames[i]),col='red',pos=4)
    ifelse(i %in% 16:20, axis(1L),NA)
    ifelse(i %in% c(1,6,11,16), axis(2,las=2),NA)
}

dev.off()
plot(position[order(position,decreasing=FALSE)],b)
text(position[order(position,decreasing=FALSE)],b,
     siteNames)

# standardizing 
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,0,0) + 0.1)
ind=(1:20)[order(position)]
b = c()
g = c()
for(i in 1:20){
  index=grep(paste0("\\[",ind[i],","),colnames(rsSummary.tmp))
  x=climate.spring[climate.spring$site==siteNames[ind[i]],]$p
  x=(x-mean(x))/(2*sd(x))
  y=rsSummary.tmp[3,index]
  g[i]=gsd.am(y)
  y=log(y)
  mod=lm(y~x)
  
  plot(x,y,xaxt='n',yaxt='n',pch=19,
       xlim=c(-2,2),ylim=c(-7,7))
  abline(v=0,lty='dotted')
  #abline(a=coef(mod)[1],b=coef(mod)[2])
  f=function(a,b,x){
    x.min=min(x);x.max=max(x);
    segments(y0=a+x.min*b,y1=a+x.max*b,x0=x.min,x1=x.max)
  }
  f(a=coef(mod)[1],b=coef(mod)[2],x=x)
  b[i]=coef(mod)[2]

  #text(.9*max(log(x+1)),.9*max(log(y+.5)),signif(cor(x,y),2),col='red')
  text(-2,-6,paste0("R2=",signif(summary(mod)$r.squared,2)),col='red',pos=4)
  text(-2,-4,paste0("b=",signif(coef(mod)[2],2)),col='red',pos=4)
  text(-2,-2,paste0(siteNames[i]),col='red',pos=4)
  ifelse(i %in% 16:20, axis(1L),NA)
  ifelse(i %in% c(1,6,11,16), axis(2,las=2),NA)
}

dev.off()
plot(position[order(position,decreasing=FALSE)],b,type='n')
text(position[order(position,decreasing=FALSE)],b,
     siteNames)

gsd.am <- function(x){
  x=x
  n = length(x[!is.na(x)])
  mu = exp(mean(log(x),na.rm=TRUE))
  y <- exp(sqrt(sum((log(x/mu))^2,na.rm=TRUE)/(n-1)))
  return(y)
}


plot(b,g,type='n',ylim=c(0,20))
text(b,g,
     siteNames)

cor(b,g)

par(mfrow = c(1,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
plot(NA,NA,pch=19,xlim=c(-250,250),ylim=c(-100,200))
for(i in 1:20){
  index=grep(paste0("\\[",i,","),colnames(rsSummary.tmp))
  x.tmp=climate.spring[climate.spring$site==siteNames[i],]$p
  x.tmp = x.tmp-mean(x.tmp)
  y.tmp=rsSummary.tmp[3,index]
  mod=lm(y.tmp~x.tmp)
  #abline(a=coef(mod)[1],b=coef(mod)[2])
  f=function(a,b,x){
    x.min=min(x);x.max=max(x);
    segments(y0=a+x.min*b,y1=a+x.max*b,x0=x.min,x1=x.max)
  }
  f(a=coef(mod)[1],b=coef(mod)[2],x=x.tmp)
 # text(.9*max(x),.9*max(y),signif(cor(x,y),2),col='red')
  #text(.9*max(x),.8*max(y),siteNames[ind[i]],col='red')
}

plot(NA,NA,pch=19,xlim=c(340,375),ylim=c(-1,1))
for(i in 1:20){
  index=grep(paste0("\\[",i,","),colnames(rsSummary.tmp))
  x.tmp=climate.spring[climate.spring$site==siteNames[i],]$p
  x.tmp = x.tmp-mean(x.tmp)
  y.tmp=rsSummary.tmp[3,index]
  mod=lm(y.tmp~x.tmp)
  #abline(a=coef(mod)[1],b=coef(mod)[2])
  f=function(a,b,x){
    x.min=min(x);x.max=max(x);
    segments(y0=a+x.min*b,y1=a+x.max*b,x0=x.min,x1=x.max)
  }
  easting=climate.spring[climate.spring$site==siteNames[i],]$easting/1000
  points(easting[1],coef(mod)[2])
  #f(a=coef(mod)[1],b=coef(mod)[2],x=x.tmp)
  # text(.9*max(x),.9*max(y),signif(cor(x,y),2),col='red')
  #text(.9*max(x),.8*max(y),siteNames[ind[i]],col='red')
}



# -------------------------------------------------------------------
# Two ways to calculate reproductive success
# 1. take posterior mode
# -------------------------------------------------------------------
sigma.mode=apply(boot::inv.logit(sigmaMCMCsamples),2,posterior.mode)
fec.mode=apply(tfeMCMCsamples,2,posterior.mode)
phi.mode=apply(observedSeedsMCMCsamples,2,posterior.mode)

dev.off()
rs1=sigma.mode*fec.mode*phi.mode

rsMCMCsamples = boot::inv.logit(sigmaMCMCsamples)*tfeMCMCsamples*observedSeedsMCMCsamples
rs2=apply(rsMCMCsamples,2,posterior.mode)

plot(rs1,rs2)

#pdf(file="~/Dropbox/clarkiaSeedBanks/products/figures/analysis/fecundity-distributions.pdf",height=16,width=8)
# par(mfrow = c(4,5),
#     oma = c(5,4,0,0) + 0.1,
#     mar = c(0,0,1,1) + 0.1)
# time.sample = 2006:2018
# for(i in 1:20){
#   index=grep(paste0("\\[",i,","),colnames(rsMCMCsamples))
#   tmp = (boot::inv.logit(sigmaMCMCsamples))[,index]
#   # tmp.f = (tfeMCMCsamples)[,index]
#   # tmp.phi = (observedSeedsMCMCsamples)[,index]
#   # tmp  = rsMCMCsamples[,index]
# 
#   xmax=c()
#   for(k in 1:13){
#     dhist=hist(tmp[,k],plot=FALSE,breaks=50)
#     xmax[k]=1.1*max(dhist$mids)
#   }
#   xmax=max(xmax)
#   plot(0,0, type="n", xlim=c(0,xmax), ylab="",xaxt='n',yaxt='n', xlab="Value",ylim=c(0,13))
#   
#   y.pt=rev(1:13)
#   for(k in 13:1){
#     
#     dhist=hist(tmp[,k],plot=FALSE,breaks=25)
#     x.max=max(tmp[,k])
#     x.min=min(tmp[,k])
#     dres <- density( tmp[,k])
#     xrange=range(dres$x);yrange=range(dres$y)
#     xrange=range(dres$x)
#     yrange=range(dres$y)
#    # polygon(dres$x, dres$y/(1.1*max(dres$y))+k-1, col="light blue")
#     segments(x0=dhist$mids,y0=y.pt[k]-1,y1=dhist$density/(1.25*max(dhist$density))+y.pt[k]-1,col='gray75')
#     points(dhist$mids,dhist$density/(1.25*max(dhist$density))+y.pt[k]-1,pch=16,cex=.5,col='gray50')
#   }
#   ifelse(i %in% 16:20, axis(1L),NA)
#   ifelse(i %in% c(1,6,11,16), axis(2,0:12,y.pt,las=2),NA)
#   
# }
# #dev.off()
# 
# pdf(file="~/Dropbox/clarkiaSeedBanks/products/figures/analysis/fecundity-distributions.pdf",height=16,width=8)
# par(mfrow = c(13,4),
#     oma = c(5,4,0,0) + 0.1,
#     mar = c(0,0,1,1) + 0.1)
# time.sample = 2006:2018
# for(i in 1:20){
#   index=grep(paste0("\\[",i,","),colnames(rsMCMCsamples))
#   tmp.sigma = (boot::inv.logit(sigmaMCMCsamples))[,index]
#   tmp.f = (tfeMCMCsamples)[,index]
#   tmp.phi = (observedSeedsMCMCsamples)[,index]
#   tmp  = rsMCMCsamples[,index]
#   for(k in 1:13){
#     hist(tmp.sigma[,k],breaks=100,main="",xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(0,1))
#     mtext(time.sample[k], side = 2,  line = 1, cex = .5, las =2)
#     hist(tmp.f[,k],breaks=100,main="",xaxt='n',yaxt='n',xlab='',ylab='')
#     hist(tmp.phi[,k],breaks=100,main="",xaxt='n',yaxt='n',xlab='',ylab='')
#     hist(tmp[,k],breaks=100,main="",xaxt='n',yaxt='n',xlab='',ylab='')
#   }
# }
# dev.off()
# 
# 
# pdf(file="~/Dropbox/clarkiaSeedBanks/products/figures/analysis/rs-distributions.pdf",height=10,width=10)
# par(mfrow = c(4,5),
#     oma = c(5,4,0,0) + 0.1,
#     mar = c(0,0,1,1) + 0.1)
# for(k in 1:13){
#   index=grep(paste0(",",k,"\\]"),colnames(rsMCMCsamples))
#   tmp  = rsMCMCsamples[,index]
#   for(i in 1:20){
#     hist(tmp[,i],breaks=100,main="",xaxt='n',yaxt='n',xlab='',ylab='')
#     abline(v=rsSummary.tmp[1,index[i]],col='red',lwd=2)
#     abline(v=rsSummary.tmp[2,index[i]],col='orange',lwd=2)
#     abline(v=rsSummary.tmp[3,index[i]],col='pink',lwd=2)
#     hist.df=hist(tmp[,i],breaks=100,plot=FALSE)
#     legend(x=.6*max(hist.df$breaks),
#            y=max(hist.df$counts),
#            legend = c("mean","median","mode"),
#            col = c("red","orange","pink"),cex=.75,
#            lty = 1,border=0)
#   }
# }
# dev.off()

par(mfrow=c(1,1))
plot(x=rep(1,260),rsSummary.tmp[1,],xlim=c(0,4),ylim=c(0,400))
points(x=rep(2,260),rsSummary.tmp[2,])
points(x=rep(3,260),rsSummary.tmp[3,])


df.list = list()
for(i in 1:20){
  obj = rsSummary.tmp
  index=grep(paste0("\\[",i,","),colnames(obj))
  tmp = signif(obj[,index],3)
  tmp.df=data.frame(site=(siteIndex[,1])[i],year=yearIndex[,1],(t(tmp)))
  names(tmp.df) = c("site","year","mu","med","mode","ci.lo95","ci.hi95","hpdi.lo95","hpdi.hi95")
  rownames(tmp.df) = NULL
  df.list[[i]] = tmp.df
}

rsSummary=do.call(rbind,df.list)


#pdf("~/Dropbox/clarkiaSeedBanks/products/figures/seedsPerFruit-population.pdf",width=8,height=6)
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1    )

time.sample = 1:13+2005
for(i in 1:20){
  
  tmp=rsSummary[rsSummary$site==siteNames[i],]

  # tmp.vec=c()
  # for(j in 1:length(time.sample)){
  #   tmp.vec[j]=sum(data$sdno[data$site3==i&data$year3==j],na.rm=TRUE)
  # }
  # 
  # tmp.vec=ifelse(is.na(tmp.vec),0,tmp.vec)
  # 
  
  
  plot(NA,NA,
       ylim=c(0,max(tmp[,3:7])),pch=16,xlim=c(2006,2019),
       ylab='',xlab='',xaxt='n',yaxt='n')
  
  # index=tmp.vec==0
  # if (sum(index)>0) {for(j in 1:length(time.sample[index])){
  #   ts=time.sample[index]
  #   polygon(x=c(ts[j]-.5,ts[j]-.5,
  #               ts[j]+.5, ts[j]+.5),
  #           y=c(-.1,max(tmp[,3:7])*2,max(tmp[,3:7])*2,-.1),col='gray95',border='gray95')
  # }} else {NA}
  # 
  # 
  # polygon(x=c(2005,2020,2020,2005),
  #         y=c(tmp.pop$ci.lo,tmp.pop$ci.lo,tmp.pop$ci.hi,tmp.pop$ci.hi),
  #         col='gray95',border='gray95')
  # 
  # abline(h=tmp.pop$med,col='gray')
  segments(time.sample,y0=tmp$ci.lo,y1=tmp$ci.hi)
  
  points(time.sample,tmp$mode,pch=21,cex=1,
         col="black",bg='pink')
  points(time.sample,tmp$mu,pch=21,cex=1,
         col="black",bg='red')
  points(time.sample,tmp$med,pch=21,cex=1,
         col="black",bg='orange')
  
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
mtext("Seeds per fruit", side = 2, outer = TRUE, line = 2.2)
#mtext("Population*year-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)
dev.off()

# could add the summary here
# histograms by year per population
par(mfrow = c(4,4),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1    )

time.sample = 1:13+2005
for(j in 1:13){

  index = grep(paste0("\\[",20,","),colnames(rsMCMCsamples))
  tmp=rsMCMCsamples[,index[j]]

  hist(tmp,freq=FALSE,breaks=100)

}






# 
# # fill in with NAs
# referenceDF <- expand.grid(site=unique(tfeFullDF$site),year=unique(tfeFullDF$year))
# 
# tfeFullDF<-tfeFullDF %>%
#   dplyr::full_join(referenceDF,by=c("site","year"))
# 
# 
# 
# data <- readRDS("~/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/data.rds")
# 
# 
# countFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countFruitsPerPlantAllPlots.RDS")
# siteIndex <- data.frame(siteIndex=unique(countFruitsPerPlantAllPlots$site),site=unique(data$site))
# yearIndex <- data.frame(yearIndex=unique(countFruitsPerPlantAllPlots$year),
#                         year=unique(data$year)) 
# 
# tfeDF<-fitness_mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(mu_py_tfe[site,year]) %>%
#   dplyr::ungroup() %>%
#   dplyr::left_join(siteIndex,by="site") %>%
#   dplyr::left_join(yearIndex,by="year") %>%
#   #dplyr::select(-c(site,year)) %>%
#   # dplyr::rename(site = siteIndex) %>%
#   # dplyr::rename(year = yearIndex) %>%
#   dplyr::mutate(mu_py_fruits=mu_py_tfe)
# 
# countSeedPerDamagedFruit <- readRDS("~/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/countSeedPerDamagedFruit.rds")
# siteIndex <- data.frame(siteIndex=unique(countSeedPerDamagedFruit$site4),site=unique(data$site4))
# yearIndex <- data.frame(yearIndex=unique(countSeedPerDamagedFruit$year4),
#                         year=unique(data$year4)) 
# 
# tfeCompDF<- fitness_mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(mu_py_tfe_comp[site,year])  %>%
#   dplyr::ungroup() %>%
#   dplyr::left_join(siteIndex,by="site") %>%
#   dplyr::left_join(yearIndex,by="year") %>%
#   #dplyr::select(-c(site,year)) %>%
#   # dplyr::rename(site = siteIndex) %>%
#   # dplyr::rename(year = yearIndex) %>%
#   dplyr::mutate(mu_py_fruits=mu_py_tfe_comp)
#   
# countUndamagedDamagedFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/countUndamagedDamagedFruitsPerPlantAllPlots.rds")
# siteIndex <- data.frame(siteIndex=unique(countUndamagedDamagedFruitsPerPlantAllPlots$site2),site=unique(data$site2))
# yearIndex <- data.frame(yearIndex=unique(countUndamagedDamagedFruitsPerPlantAllPlots$year2),
#                         year=unique(data$year2)) 
# 
# ufDF<-fitness_mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(mu_py_und[site,year])   %>%
#   dplyr::ungroup() %>%
#   dplyr::left_join(siteIndex,by="site") %>%
#   dplyr::left_join(yearIndex,by="year") %>%
#   #dplyr::select(-c(site,year)) %>%
#   # dplyr::rename(site = siteIndex) %>%
#   # dplyr::rename(year = yearIndex)  %>%
#   dplyr::mutate(mu_py_fruits=mu_py_und)
# 
# ufRef<-ufDF %>%
#   dplyr::select(siteIndex,yearIndex) %>%
#   unique() %>% 
#   dplyr::mutate(site.year=paste0(siteIndex,".",yearIndex))
# tfeCompRef<-tfeCompDF %>%
#   dplyr::select(siteIndex,yearIndex) %>%
#   unique() %>% 
#   dplyr::mutate(site.year=paste0(siteIndex,".",yearIndex))
# tfeRef<-tfeDF %>%
#   dplyr::select(siteIndex,yearIndex) %>%
#   unique() %>% 
#   dplyr::mutate(site.year=paste0(siteIndex,".",yearIndex))
# 
# # find index with appropriate samples
# findIndex<-function(i=siteNumber,j=yearNumber){
#   site<-siteIndex$siteIndex[i]
#   year<-(2006:2018)[j]
#    paste0(site,".",year)
# }
# 
# # find corresponding matrix
# findMat<-function(i=siteNumber,j=yearNumber){
# if (findIndex(i,j) %in% tfeRef$site.year) {
#   return(tfeDF)
# } else if (findIndex(i,j) %in% tfeCompRef$site.year) { 
#   return(tfeCompDF)
# } else {
#   return(ufDF)
# }
# }
# 
# # get samples corresponding to particular site and years
# f<-function(chains=zc,site=i,year=j){
#   tmp <- chains %>% dplyr::filter(site==i&year==j)
#   return(tmp)
# }
# 
# rsEstimates<-list()
# #iter<-max(dim(mcmcSigma)[1],dim(mcmcFec)[1],dim(mcmcPhi)[1])
# iter=1000
# rs<-matrix(NA,nrow=iter,ncol=13)
# 
# # change if number of years changes
# for(i in 1:20){
#   for(j in 1:13){
#     sigma=f(chains=mcmcSigma,site=i,year=j)$p_1
#     fec=f(chains=findMat(i,j),site=i,year=j)$mu_py_fruits
#     phi=f(chains=mcmcPhi,site=i,year=j)$mu_py_seeds
#     
#     # ignore missing data for now?
#     sigma2 = if(length(sigma)>0) {sample(sigma,iter)} else {rep(NA,iter)}
#     fec2= if(length(fec)>0) {sample(fec,iter)} else {rep(NA,iter)}
#     phi2= if(length(phi)>0) {sample(phi,iter)} else {rep(NA,iter)}
#     
#     rs[,j] <- sigma2*fec2*phi2
#   }
#   rsEstimates[[i]] <- rs
# }
# 
# saveRDS(rsEstimates,file="~/Dropbox/clarkiaSeedBanks/products/dataFiles/rsVarFullPosterior.RDS")

# # -------------------------------------------------------------------
# # Import data
# # -------------------------------------------------------------------
# sigmaSummary<-read.csv("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/parameterSummary/sigmaSummary.csv")[,-1]
# fecSummary<-read.csv("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/parameterSummary/fruitEquivalentsPerPlantSummary.csv")[,-1]
# phiSummary<-read.csv("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/parameterSummary/seedsSummary.csv")[,-1]
# 
# # -------------------------------------------------------------------
# # Reproductive success: medians
# # -------------------------------------------------------------------
# nsite <- 20
# nyear <- 6
# 
# sigmaSummary <- sigmaSummary %>%
#   dplyr::select(site,year,med) %>%
#   dplyr::rename(med_sigma = med)
# 
# fecSummary <- fecSummary %>%
#   dplyr::select(site,year,med) %>%
#   dplyr::rename(med_fec = med)
# 
# phiSummary <- phiSummary %>%
#   dplyr::select(site,year,med) %>%
#   dplyr::rename(med_phi = med)
# 
# abovegroundMedians <- sigmaSummary %>% 
#   dplyr::full_join(fecSummary,by=c("site","year")) %>%
#   dplyr::full_join(phiSummary,by=c("site","year")) %>%
#   dplyr::filter(year<2019)
# 
# rsMedianEstimates <- abovegroundMedians %>% 
#   dplyr::mutate(rs = med_sigma*med_fec*med_phi) %>%
#   dplyr::filter(year<2019)
# 
# # save estimates of RS with sites/years where there is missing data excluded
# saveRDS(rsMedianEstimates,file="~/Dropbox/clarkiaSeedBanks/products/dataFiles/rsMedianEstimates.RDS")
