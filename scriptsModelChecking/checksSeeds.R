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
  dplyr::filter(siteBags=="GCN") %>%
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
zc.df<-zc[,siteVec==(1:20)[siteNames=="GCN"]]

siteVec=extractChains(zc.s1)[[1]];yearVec=extractChains(zc.s1)[[2]]
zc.s1.df<-zc.s1[,siteVec==(1:20)[siteNames=="GCN"]]


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

#g1s1
out=df %>% dplyr::filter(year==2006) %>% dplyr::mutate(p=seedlingJan/seedStart)
hist(zc.df[,1]*zc.s1.df[,1],breaks=100);points(x=out$p,y=rep(0,dim(out)[1]),col='red')

out=df %>% dplyr::filter(year==2007) %>% dplyr::mutate(p=seedlingJan/seedStart)
hist(zc.df[,2]*zc.s1.df[,2],breaks=100);points(x=out$p,y=rep(0,dim(out)[1]),col='red')

out=df %>% dplyr::filter(year==2008) %>% dplyr::mutate(p=seedlingJan/seedStart)
hist(zc.df[,3]*zc.s1.df[,3],breaks=100);points(x=out$p,y=rep(0,dim(out)[1]),col='red')


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
  geom_point(data=df,aes(x=yearBags,y=seedlingJan/100)) +
  geom_point(data=df %>%
               dplyr::group_by(yearBags) %>%
               dplyr::summarise(p.mu = mean(seedlingJan)/100),
             aes(x=yearBags,y=p.mu),size=2,color='red') +
  theme_bw()

med<-zc.df %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(p.med = median(value))

df.med <- df %>%
  dplyr::group_by(yearBags) %>%
  dplyr::summarise(p.mu = mean(seedlingJan)/100)

df=df %>%
  dplyr::rename(year=yearBags) %>%
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
  


plot(df.join$p.med,df.join$seedlingJan/df.join$seedStart,xlim=c(0,.5),ylim=c(0,.5));abline(a=0,b=1)

plot(med$p.med,df.med$p.mu,xlim=c(0,.5),ylim=c(0,.5));abline(a=0,b=1)


################################################################################
# Seedling counts for multiple sites
# now s1g1
#################################################################################

df <- seedBagExperiment

siteNames <- unique(df$siteBags)

ggplot() +
  geom_point(data=df,aes(x=siteBags,y=seedlingJan)) +
  geom_point(data=df %>%
               dplyr::group_by(siteBags,yearBags) %>%
               dplyr::summarise(p.mu = mean(seedlingJan)),
             aes(x=siteBags,y=p.mu),size=2,color='red') +
  facet_wrap(~yearBags)

color_scheme_set("brightblue")

zc <- MCMCchains(mcmcSamples,params="g1.0")
zc.s1 <-  MCMCchains(mcmcSamples,params="s1.0")

extractChains<-function(x){
  spl <- strsplit(sub("]", "[", colnames(x)), "[", fixed = TRUE)
  ind <- strsplit(sapply(spl, "[", 2), ",", fixed = TRUE)
  siteVec <- sapply(ind, "[", 1)
  yearVec <- sapply(ind, "[", 2)
  return(list(siteVec,yearVec))
}

siteVec=extractChains(zc)[[1]];yearVec=extractChains(zc)[[2]]
# zc.df<-zc[,siteVec==(1:20)[siteNames=="GCN"]]
# zc.s1.df<-zc.s1[,siteVec==(1:20)[siteNames=="GCN"]]

zc.df <- list()
zc.s1.df <- list()
for(i in 1:20){
  tmp = zc[,siteVec==i]
  tmp.s1 = zc.s1[,siteVec==i]
  tmp = tmp*tmp.s1
  colnames(tmp) = c('2006','2007','2008')
  tmp = tmp %>%
    data.frame() %>%
    tidyr::pivot_longer(cols=1:3) %>%
    dplyr::mutate(year = ifelse(name=="X2006",2006,
                                ifelse(name=="X2007",2007,2008))) %>%
    dplyr::mutate(site = siteNames[i]) %>%
    dplyr::select(site,year,value) 
  zc.df[[i]]<-tmp
}

zc.df<-data.table::rbindlist(zc.df)


med<-zc.df %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(p.med = median(value))

df.med <- df %>%
  dplyr::group_by(siteBags,yearBags) %>%
  dplyr::summarise(p.mu = mean(seedlingJan)/100)

df=df %>%
  dplyr::rename(year=yearBags,site=siteBags) %>%
  dplyr::mutate(year=as.numeric(as.character(year)))

df.join<-df %>% dplyr::left_join(med,by=c('site','year'))

ggplot() +
  geom_point(data=df.join,aes(x=p.med,y=seedlingJan/seedStart),
              alpha=0.5) +
  geom_point(data=df.join %>%
               dplyr::group_by(site,year,p.med) %>%
               dplyr::summarise(p.mu = mean(seedlingJan)/100),
             aes(x=p.med,y=p.mu),size=2,color='red') +
  geom_abline(intercept=0,slope=1) +
#  xlim(c(0,1)) + ylim(c(0,1)) +
  theme_bw() +
  facet_wrap(~year,scales='free')

ggplot() +
  geom_point(data=df.join,aes(x=p.med,y=seedlingJan/seedStart),
             alpha=0.5) +
  geom_point(data=df.join %>%
               dplyr::group_by(site,year,p.med) %>%
               dplyr::summarise(p.mu = mean(seedlingJan)/100),
             aes(x=p.med,y=p.mu,color=as.factor(year)),size=2) +
  geom_abline(intercept=0,slope=1) +
  theme_bw() +
  facet_wrap(~site,scales='free')


zc.pool <- MCMCchains(mcmcSamples,params="g1")
zc.pool.s1 <- MCMCchains(mcmcSamples,params="s1")
zc.pool=zc.pool*zc.pool.s1
colnames(zc.pool) = siteNames

zc.pool = zc.pool %>% 
  data.frame() %>%
  tidyr::pivot_longer(1:20) %>%
  dplyr::rename(site=name)

ggplot(data=zc.pool) +
  geom_boxplot(aes(x=site,y=value))

zc.pool.sum = zc.pool %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(site.med = median(value))

ggplot() +
  geom_hline(data=zc.pool.sum,aes(yintercept=site.med),linetype='dotted') +
  geom_point(data=df.join,aes(x=p.med,y=seedlingJan/seedStart),
             alpha=0.5) +
  geom_point(data=df.join %>%
               dplyr::group_by(site,year,p.med) %>%
               dplyr::summarise(p.mu = mean(seedlingJan)/100),
             aes(x=p.med,y=p.mu,color=as.factor(year)),size=2) +
  geom_abline(intercept=0,slope=1) +
  theme_bw() +
  facet_wrap(~site,scales='free')

plot(zc.pool$site,zc.pool$site.med)
## What explains the high residuals in S22 2008 and GCN 2006?

residuals = df.join %>%
  dplyr::group_by(site,year,p.med) %>%
  dplyr::summarise(p.hat = mean(seedlingJan)/100) %>%
  dplyr::mutate(resid = p.med-p.hat)

## Check if these sites*years had few seed bags

residuals=df.join %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::left_join(residuals,by=c('site','year'))

plot(residuals$resid,residuals$count)

## check if these sites*years had high variance

residuals=df.join %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(var = var(seedlingJan)) %>%
  dplyr::left_join(residuals,by=c('site','year'))

plot(residuals$resid,residuals$var)

## check if the site*years just had few seedlings?

residuals=df.join %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(tot = sum(seedlingJan)) %>%
  dplyr::left_join(residuals,by=c('site','year'))

plot(residuals$resid,residuals$tot)

# high resid because high number of seeds
# check the site level estimate - seems like S22 2008 is same as prior?

# higher residuals for higher median estimates

plot(residuals$resid,residuals$p.med);abline(a=0,b=1)
