rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)

library(MCMCvis)
library(tidybayes)
library(tidyverse)
library(magrittr)

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/fruitsPerPlant/"
simFiles <- paste0(directory,list.files(directory))
dirFigures = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/"

samples.rjags <- readRDS(simFiles[[3]])

data <- readRDS(simFiles[[2]])
countFruitsPerPlantTransects <- readRDS(simFiles[[1]])

# # -------------------------------------------------------------------
# # -------------------------------------------------------------------
# # Visualize
# # -------------------------------------------------------------------
# # -------------------------------------------------------------------
summary<-apply(exp(MCMCchains(samples.rjags,params="gamma")),2,quantile, probs = c(.05,.5,.95))

summaryList<-list(summary)

#pdf(paste0(dirFigures,"shrinkage-fruitsPerPlant.pdf"), width=8, height=6)
est=countFruitsPerPlantTransects %>%
  dplyr::group_by(year) %>% 
  dplyr::summarise(mu=mean(countFruitsPerPlant))

grandMean = countFruitsPerPlantTransects %>%
  dplyr::summarise(mu = mean(countFruitsPerPlant))

par(mfrow=c(1,1))
for(i in 1:1){
  plot(x = est$mu,
       y = t(summaryList[[i]])[,2],
       xlab= "Mean count of fruits per plant",
       ylab="MCMC estimate and 95% CI for probability",
       pch=16,cex=.5,
       xlim=c(0,7),ylim=c(0,7))
  segments(x0=est$mu,
           y0= t(summaryList[[i]])[,1],
           y1= t(summaryList[[i]])[,3])
  abline(a=0,b=1)
  #abline(h=grandMean$mu)
}

summary2<-apply(exp(MCMCchains(samples.rjags,params="mu0")),2,quantile, probs = c(.05,.5,.95))

grandMeanAll = exp(rnorm(10000,MCMCchains(samples.rjags,params="mu0"),MCMCchains(samples.rjags,params="sigma0")))
#summary2<-quantile(grandMeanAll,probs = c(.05,.5,.95))

points(x=grandMean$mu,y=t(summary2)[,2],pch=16)
segments(x0=grandMean$mu,
         y0= t(summary2)[,1],
         y1= t(summary2)[,3])
#dev.off()




# fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBurial/")
# dir.create(file.path(fileDirectory), showWarnings = FALSE)
# 
# saveRDS(samples.rjags,file=paste0(fileDirectory,"samples.rjags.rds"))

data$freq <- data$seedlingJan/data$totalJan

summary.LogitCentered<-apply(MCMCchains(samples.rjags,params="theta_2"),2,quantile, probs = c(.05,.5,.95))

summaryList<-list(summary.LogitCentered)

#pdf(paste0(dirFigures,"shrinkage.pdf"), width=8, height=6)
est=seedBagsData %>%
  dplyr::group_by(siteBags) %>% 
  dplyr::summarise(totalJan=sum(totalJan),seedlingJan=sum(seedlingJan)) %>%
  dplyr::mutate(prop = seedlingJan/totalJan)
par(mfrow=c(1,1))
for(i in 1:1){
  plot(x = data$freq,
       y = t(summaryList[[i]])[,2],
       xlab= "Proportion of germinants remaining in bag",
       ylab="MCMC estimate and 95% CI for probability",
       pch=16,cex=.5,
       xlim=c(0,1),ylim=c(0,1),
       col=data$siteBags)
  segments(x0=data$freq,
           y0= t(summaryList[[i]])[,1],
           y1= t(summaryList[[i]])[,3],
           col=data$siteBags)
  abline(a=0,b=1)
  abline(h=est$prop,col=as.factor(est$siteBags))
}
#dev.off()

# germination trials
data$freq <- data$germCount/data$germStart

summary.LogitCentered<-apply(MCMCchains(samples.rjags,params="theta_g"),2,quantile, probs = c(.05,.5,.95))

summaryList<-list(summary.LogitCentered)

#pdf(paste0(dirFigures,"shrinkage.pdf"), width=8, height=6)
est=viabilityRawData %>%
  dplyr::group_by(siteViab) %>% 
  dplyr::summarise(germCount=sum(germCount),germStart=sum(germStart)) %>%
  dplyr::mutate(prop = germCount/germStart)
par(mfrow=c(1,1))
for(i in 1:1){
  plot(x = data$freq,
       y = t(summaryList[[i]])[,2],
       xlab= "Proportion of germinants in test",
       ylab="MCMC estimate and 95% CI for probability",
       pch=16,cex=.5,
       xlim=c(0,1),ylim=c(0,1),
       col=data$siteViab)
  segments(x0=data$freq,
           y0= t(summaryList[[i]])[,1],
           y1= t(summaryList[[i]])[,3],
           col=data$siteViab)
  abline(a=0,b=1)
  abline(h=est$prop,col=as.factor(est$siteViab))
}

# viability trials
data$freq <- data$viabStain/data$viabStart

summary.LogitCentered<-apply(MCMCchains(samples.rjags,params="theta_v"),2,quantile, probs = c(.05,.5,.95))

summaryList<-list(summary.LogitCentered)

#pdf(paste0(dirFigures,"shrinkage.pdf"), width=8, height=6)
est=viabilityRawData %>%
  dplyr::group_by(siteViab) %>% 
  dplyr::summarise(viabStain=sum(viabStain,na.rm=TRUE),viabStart=sum(viabStart)) %>%
  dplyr::mutate(prop = viabStain/viabStart)
par(mfrow=c(1,1))
for(i in 1:1){
  plot(x = data$freq,
       y = t(summaryList[[i]])[,2],
       xlab= "Proportion of germinants in test",
       ylab="MCMC estimate and 95% CI for probability",
       pch=16,cex=.5,
       xlim=c(0,1),ylim=c(0,1),
       col=data$siteViab)
  segments(x0=data$freq,
           y0= t(summaryList[[i]])[,1],
           y1= t(summaryList[[i]])[,3],
           col=data$siteViab)
  abline(a=0,b=1)
  abline(h=est$prop,col=as.factor(est$siteViab))
}
