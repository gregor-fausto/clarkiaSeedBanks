rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)

library(MCMCvis)
library(tidybayes)
library(tidyverse)
library(magrittr)
library(bayesplot)
library(rethinking)

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedlingSurvival/"
modelFittingFiles <- paste0(directory,list.files(directory))

mcmcSamples <- readRDS(modelFittingFiles[[grep("seedSurvivalSamplesChecksPooling.rds",modelFittingFiles)]])
data <- readRDS(modelFittingFiles[[grep("data.rds",modelFittingFiles)]])


censusSeedlingsFruitingPlants <- readRDS("~/Dropbox/dataLibrary/postProcessingData/censusSeedlingsFruitingPlants.RDS")

siteNames = unique(censusSeedlingsFruitingPlants$site)

MCMCsummary(mcmcSamples,params="e.mu")
e.mu=MCMCchains(mcmcSamples,params="e.mu")

par(mfrow=c(4,4))
for(i in 1:14){
  hist(e.mu[,i],breaks=100)
  abline(v=0,lwd=2,col='red')
}

lambda.mu <- 1-var(apply(e.mu,2,mean))/mean(apply(e.mu,1,var))

mu=MCMCchains(mcmcSamples,params="mu")

p.pop=c()
eps.mat = matrix(NA,nrow=dim(mu)[1],ncol=14)
par(mfrow=c(4,4))
for(i in 1:14){
  index=data$site==1&data$year==i
  p.pop[i]=sum(data$fruitplNumber[index],na.rm=TRUE)/sum(data$seedlingNumber[index],na.rm=TRUE)
  mu.pop=boot::inv.logit(mu[,i])
  eps.mat[,i]=mu.pop
 # hist(eps.mat[,i],breaks=100)
 # abline(v=0,col='red',lwd=2)
}
plot(p.pop,apply(eps.mat,2,mean),xlim=c(0,1),ylim=c(0,1))
abline(a=0,b=1)
