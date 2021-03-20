#################################################################################
################################################################################
################################################################################
# Code for figures to compare the following modeling approaches for the seed bag data
# 
# Scripts by Gregor Siegmund
# fausto.siegmund@gmail.com
# last updated 04-22-2020
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
library(gridExtra)


directory = "/Users/Gregor/Dropbox/dataLibrary/posteriors/"
simFiles <- paste0(directory,list.files(directory))

index<-grep("belowgroundSamplesAllYears.RDS",simFiles)

mcmcSamples <- readRDS(simFiles[[index]])

################################################################################
# Data directory
#################################################################################

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBurial/"
dataFiles <- paste0(directory,list.files(directory))

data <- readRDS(dataFiles[[1]])
seedBagExperiment <- readRDS(dataFiles[[2]])

################################################################################
# P(G|S)
#################################################################################

df <- seedBagExperiment %>%
  dplyr::filter(siteBags=="S22") %>%
  dplyr::rename(site=siteBags,year=yearBags)

ggplot() +
  geom_point(data=df,aes(x=year,y=seedlingJan)) +
  geom_point(data=df %>%
               dplyr::group_by(year) %>%
               dplyr::summarise(p.mu = mean(seedlingJan)),
             aes(x=year,y=p.mu),size=2,color='red')

################################################################################
# Test seedling counts for one site
#################################################################################

color_scheme_set("brightblue")

zc <- MCMCchains(mcmcSamples,params="g1.0")
zc.s1 <-  MCMCchains(mcmcSamples,params="s1.0")

siteNames <- unique(seedBagExperiment$siteBags)

extractChains<-function(x){
  spl <- strsplit(sub("]", "[", colnames(x)), "[", fixed = TRUE)
  ind <- strsplit(sapply(spl, "[", 2), ",", fixed = TRUE)
  siteVec <- sapply(ind, "[", 1)
  yearVec <- sapply(ind, "[", 2)
  return(list(siteVec,yearVec))
}
siteVec=extractChains(zc)[[1]];yearVec=extractChains(zc)[[2]]
zc.df<-zc[,siteVec==(1:20)[siteNames=="S22"]]

siteVec=extractChains(zc.s1)[[1]];yearVec=extractChains(zc.s1)[[2]]
zc.s1.df<-zc.s1[,siteVec==(1:20)[siteNames=="S22"]]


par(mfrow=c(3,3))
out=df %>% dplyr::filter(year==2006) %>% dplyr::mutate(p=seedlingJan/seedStart)
hist(zc.df[,1],breaks=100);points(x=out$p,y=rep(0,dim(out)[1]),col='red')

out=df %>% dplyr::filter(year==2007) %>% dplyr::mutate(p=seedlingJan/seedStart)
hist(zc.df[,2],breaks=100);points(x=out$p,y=rep(0,dim(out)[1]),col='red')

out=df %>% dplyr::filter(year==2008) %>% dplyr::mutate(p=seedlingJan/seedStart)
hist(zc.df[,3],breaks=100);points(x=out$p,y=rep(0,dim(out)[1]),col='red')

# s1
hist(zc.s1.df[,1],breaks=100);
hist(zc.s1.df[,2],breaks=100);
hist(zc.s1.df[,3],breaks=100);
vitalRates$g1*vitalRates$s1

#g1s1
out=df %>% dplyr::filter(year==2006) %>% dplyr::mutate(p=seedlingJan/seedStart)
hist(zc.df[,1]*zc.s1.df[,1],breaks=100);points(x=out$p,y=rep(0,dim(out)[1]),col='red')

out=df %>% dplyr::filter(year==2007) %>% dplyr::mutate(p=seedlingJan/seedStart)
hist(zc.df[,2]*zc.s1.df[,2],breaks=100);points(x=out$p,y=rep(0,dim(out)[1]),col='red')

out=df %>% dplyr::filter(year==2008) %>% dplyr::mutate(p=seedlingJan/seedStart)
hist(zc.df[,3]*zc.s1.df[,3],breaks=100);points(x=out$p,y=rep(0,dim(out)[1]),col='red')


median(zc.df[,1]*zc.s1.df[,1])
median(zc.df[,2]*zc.s1.df[,2])
median(zc.df[,3]*zc.s1.df[,3])


iter<-dim(zc.df)[1]

colnames(zc.df) = c('2006','2007','2008')
zc.df = zc.df %>%
  data.frame() %>%
  tidyr::pivot_longer(cols=1:3) %>%
  dplyr::mutate(year = ifelse(name=="X2006",2006,
                              ifelse(name=="X2007",2007,2008)))


#boxplot(zc.df[sample(iter,1000), ][,1],zc.df[sample(iter,1000), ][,2],zc.df[sample(iter,1000), ][,3])

ggplot() +
  geom_boxplot(data=zc.df[sample(iter,1000),], aes(x=as.factor(year),y=value),size=0.25,alpha=0.5) +
  geom_point(data=df,aes(x=year,y=seedlingJan/100)) +
  geom_point(data=df %>%
               dplyr::group_by(year) %>%
               dplyr::summarise(p.mu = mean(seedlingJan)/100),
             aes(x=year,y=p.mu),size=2,color='red') +
  theme_bw()

med<-zc.df %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(p.med = median(value))

df.med <- df %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(p.mu = mean(seedlingJan)/100)

df=df %>%
  dplyr::rename(year=year) %>%
  dplyr::mutate(year=as.numeric(as.character(year)))

df.join<-df %>% dplyr::left_join(med,by='year')

ggplot() +
  geom_jitter(data=df.join,aes(x=p.med,y=seedlingJan/seedStart),
              width=0.005,alpha=0.5,height=0) +
  geom_point(data=df.join %>%
               dplyr::group_by(year,p.med) %>%
               dplyr::summarise(p.mu = mean(seedlingJan)/100),
             aes(x=p.med,y=p.mu),size=2,color='red') +
  geom_abline(intercept=0,slope=1) +
  xlim(c(0,0.6)) + ylim(c(0,0.6)) +
  theme_bw()

par(mfrow=c(1,1))
plot(df.join$p.med,df.join$seedlingJan/df.join$seedStart,xlim=c(0,.5),ylim=c(0,.5));abline(a=0,b=1)

plot(med$p.med,df.med$p.mu,xlim=c(0,.5),ylim=c(0,.5));abline(a=0,b=1)


################################################################################
# Load Am Nat estimates
#################################################################################
vitalRates <- readRDS(file="~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/modelCheckingSeeds/eckhartEstimates.RDS")

names(vitalRates)
vr.df<-do.call(rbind, vitalRates)
vr.df<-t(vr.df)
colnames(vr.df) = names(vitalRates)
vr.df = data.frame(vr.df[,c("g1","s1","s2")])

zc.g1 <- MCMCchains(mcmcSamples,params="g1.0")
zc.s1 <-  MCMCchains(mcmcSamples,params="s1.0")
zc.s2 <-  MCMCchains(mcmcSamples,params="s2.0")

extractChains<-function(x){
  spl <- strsplit(sub("]", "[", colnames(x)), "[", fixed = TRUE)
  ind <- strsplit(sapply(spl, "[", 2), ",", fixed = TRUE)
  siteVec <- sapply(ind, "[", 1)
  yearVec <- sapply(ind, "[", 2)
  return(list(siteVec,yearVec))
}
siteVec=extractChains(zc)[[1]];yearVec=extractChains(zc)[[2]]
zc.g1.df<-zc.g1[,siteVec==(1:20)[siteNames=="S22"]]
zc.s1.df<-zc.s1[,siteVec==(1:20)[siteNames=="S22"]]
zc.s2.df<-zc.s2[,siteVec==(1:20)[siteNames=="S22"]]


## Total January == s1
par(mfrow=c(1,3))
out=df %>% dplyr::filter(year==2006) %>% dplyr::mutate(p=totalJan/seedStart)
hist(zc.s1.df[,1],breaks=100);points(x=out$p,y=rep(0,dim(out)[1]),col='red')
abline(v=vr.df$s1[1],col='red',lty='dotted');abline(v=median(zc.s1.df[,1]),col='orange')

out=df %>% dplyr::filter(year==2007) %>% dplyr::mutate(p=totalJan/seedStart)
hist(zc.s1.df[,2],breaks=100);points(x=out$p,y=rep(0,dim(out)[1]),col='red')
abline(v=vr.df$s1[2],col='red',lty='dotted');abline(v=median(zc.s1.df[,2]),col='orange')

out=df %>% dplyr::filter(year==2008) %>% dplyr::mutate(p=totalJan/seedStart)
hist(zc.s1.df[,3],breaks=100);points(x=out$p,y=rep(0,dim(out)[1]),col='red')
abline(v=vr.df$s1[3],col='red',lty='dotted');abline(v=median(zc.s1.df[,3]),col='orange')

## Seedlings January == s1*g1
par(mfrow=c(1,3))
out=df %>% dplyr::filter(year==2006) %>% dplyr::mutate(p=seedlingJan/seedStart)
hist(zc.g1.df[,1]*zc.s1.df[,1],breaks=100);points(x=out$p,y=rep(0,dim(out)[1]),col='red')
abline(v=vr.df$g1[1]*vr.df$s1[1],col='red',lty='dotted');abline(v=median(zc.g1.df[,1]*zc.s1.df[,1]),col='orange')

out=df %>% dplyr::filter(year==2007) %>% dplyr::mutate(p=seedlingJan/seedStart)
hist(zc.g1.df[,2]*zc.s1.df[,2],breaks=100);points(x=out$p,y=rep(0,dim(out)[1]),col='red')
abline(v=vr.df$g1[2]*vr.df$s1[2],col='red',lty='dotted');abline(v=median(zc.g1.df[,2]*zc.s1.df[,2]),col='orange')

out=df %>% dplyr::filter(year==2008) %>% dplyr::mutate(p=seedlingJan/seedStart)
hist(zc.g1.df[,3]*zc.s1.df[,3],breaks=100);points(x=out$p,y=rep(0,dim(out)[1]),col='red')
abline(v=vr.df$g1[3]*vr.df$s1[3],col='red',lty='dotted');abline(v=median(zc.g1.df[,3]*zc.s1.df[,3]),col='orange')

## October intact == s1*(1-g1)*s2
par(mfrow=c(1,3))
out=df %>% dplyr::filter(year==2006) %>% dplyr::mutate(p=intactOct/seedStart)
hist(zc.s2.df[,1]*(1-zc.g1.df[,1])*zc.s1.df[,1],breaks=100);points(x=out$p,y=rep(0,dim(out)[1]),col='red')
abline(v=vr.df$s2[1]*(1-vr.df$g1[1])*vr.df$s1[1],col='red',lty='dotted');abline(v=median(zc.s2.df[,1]*(1-zc.g1.df[,1])*zc.s1.df[,1]),col='orange')

out=df %>% dplyr::filter(year==2007) %>% dplyr::mutate(p=intactOct/seedStart)
hist(zc.s2.df[,2]*(1-zc.g1.df[,2])*zc.s1.df[,2],breaks=100);points(x=out$p,y=rep(0,dim(out)[1]),col='red')
abline(v=vr.df$s2[2]*(1-vr.df$g1[2])*vr.df$s1[2],col='red',lty='dotted');abline(v=median(zc.s2.df[,2]*(1-zc.g1.df[,2])*zc.s1.df[,2]),col='orange')

out=df %>% dplyr::filter(year==2008) %>% dplyr::mutate(p=intactOct/seedStart)
hist(zc.s2.df[,3]*(1-zc.g1.df[,3])*zc.s1.df[,3],breaks=100);points(x=out$p,y=rep(0,dim(out)[1]),col='red')
abline(v=vr.df$s2[3]*(1-vr.df$g1[3])*vr.df$s1[3],col='red',lty='dotted');abline(v=median(zc.s2.df[,3]*(1-zc.g1.df[,3])*zc.s1.df[,3]),col='orange')


## compare all estimates
par(mfrow=c(1,1))
plot(vr.df$g1,apply(zc.g1.df,2,median),xlim=c(0,1),ylim=c(0,1),type='n')
abline(a=0,b=1,lty='dotted')
points(vr.df$g1,apply(zc.g1.df,2,median),col='red')
points(vr.df$s1,apply(zc.s1.df,2,median),col='orange')
points(vr.df$s2,apply(zc.s2.df,2,median),col='black')

