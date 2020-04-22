rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)

library(MCMCvis)
library(tidybayes)
library(tidyverse)
library(magrittr)

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBurial/"
simFiles <- paste0(directory,list.files(directory))
dirFigures = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/")

samples.rjags <- readRDS(simFiles[[3]])

data <- readRDS(simFiles[[1]])
seedBagsData <- readRDS(simFiles[[2]])
viabilityRawData <- readRDS(simFiles[[4]])

# # -------------------------------------------------------------------
# # -------------------------------------------------------------------
# # Visualize
# # -------------------------------------------------------------------
# # -------------------------------------------------------------------

plot(density((boot::inv.logit(MCMCchains(samples.rjags,params="s1")[,1]))),
     ylim=c(0,10),type='n',xlim=c(0,1))
for(i in 1:2){
  lines(density(MCMCchains(samples.rjags,params="s1")[,i]),col=i)
  lines(density(MCMCchains(samples.rjags,params="s2")[,i]),col=i,lty="dashed")
}


# fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBurial/")
# dir.create(file.path(fileDirectory), showWarnings = FALSE)
# 
# saveRDS(samples.rjags,file=paste0(fileDirectory,"samples.rjags.rds"))

data$freq <- data$totalJan/data$seedStart

summary.LogitCentered<-apply(MCMCchains(samples.rjags,params="theta_1"),2,quantile, probs = c(.05,.5,.95))

summaryList<-list(summary.LogitCentered)

#pdf(paste0(dirFigures,"shrinkage.pdf"), width=8, height=6)
est=seedBagsData %>%
  dplyr::group_by(siteBags) %>% 
  dplyr::summarise(totalJan=sum(totalJan),seedStart=sum(seedStart)) %>%
  dplyr::mutate(prop = totalJan/seedStart)
par(mfrow=c(1,1))
for(i in 1:1){
  plot(x = data$freq,
       y = t(summaryList[[i]])[,2],
       xlab= "Proportion of seeds remaining in bag",
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

plot(density((boot::inv.logit(MCMCchains(samples.rjags,params="s1")[,1]))),
     ylim=c(0,10),type='n',xlim=c(0,1))
for(i in 1:2){
  lines(density(MCMCchains(samples.rjags,params="s1")[,i]),col=i)
  lines(density(MCMCchains(samples.rjags,params="s2")[,i]),col=i,lty="dashed")
}


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
