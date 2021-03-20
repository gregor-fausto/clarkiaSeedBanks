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
  dplyr::rename(site=siteBags,year=yearBags)

ggplot() +
  geom_point(data=df,aes(x=year,y=seedlingJan)) +
  geom_point(data=df %>%
               dplyr::group_by(year) %>%
               dplyr::summarise(p.mu = mean(seedlingJan)),
             aes(x=year,y=p.mu),size=2,color='red')

################################################################################
# Comparison for multiple sites
#################################################################################

siteNames <- unique(df$site)

ggplot() +
  geom_point(data=df,aes(x=site,y=seedlingJan)) +
  geom_point(data=df %>%
               dplyr::group_by(site,site) %>%
               dplyr::summarise(p.mu = mean(seedlingJan)),
             aes(x=site,y=p.mu),size=2,color='red') +
  facet_wrap(~site)

color_scheme_set("brightblue")

zc.g1 <- MCMCchains(mcmcSamples,params="g1.0")
zc.s1 <-  MCMCchains(mcmcSamples,params="s1.0")
zc.s2 <-  MCMCchains(mcmcSamples,params="s2.0")
zc.nu <- MCMCchains(mcmcSamples,params="nu0_1")

extractChains<-function(x){
  spl <- strsplit(sub("]", "[", colnames(x)), "[", fixed = TRUE)
  ind <- strsplit(sapply(spl, "[", 2), ",", fixed = TRUE)
  siteVec <- sapply(ind, "[", 1)
  yearVec <- sapply(ind, "[", 2)
  return(list(siteVec,yearVec))
}

siteVec=extractChains(zc.g1)[[1]];yearVec=extractChains(zc.g1)[[2]]

zc.g1.df <- list()
zc.s1.df <- list()
zc.s2.df <- list()
zc.nu.df <- list()

for(i in 1:20){
  tmp.g1 = zc.g1[,siteVec==i]
  tmp.s1 = zc.s1[,siteVec==i]
  tmp.s2 = zc.s2[,siteVec==i]
  tmp.nu = zc.nu[,siteVec==i]
  colnames(tmp.g1) = colnames(tmp.s1) = colnames(tmp.s2) = colnames(tmp.nu) = c('2006','2007','2008')
  tmp.g1 = tmp.g1 %>%
    data.frame() %>%
    tidyr::pivot_longer(cols=1:3) %>%
    dplyr::mutate(year = ifelse(name=="X2006",2006,
                                ifelse(name=="X2007",2007,2008))) %>%
    dplyr::mutate(site = siteNames[i]) %>%
    dplyr::select(site,year,value) 
  tmp.s1 = tmp.s1 %>%
    data.frame() %>%
    tidyr::pivot_longer(cols=1:3) %>%
    dplyr::mutate(year = ifelse(name=="X2006",2006,
                                ifelse(name=="X2007",2007,2008))) %>%
    dplyr::mutate(site = siteNames[i]) %>%
    dplyr::select(site,year,value) 
  tmp.s2 = tmp.s2 %>%
    data.frame() %>%
    tidyr::pivot_longer(cols=1:3) %>%
    dplyr::mutate(year = ifelse(name=="X2006",2006,
                                ifelse(name=="X2007",2007,2008))) %>%
    dplyr::mutate(site = siteNames[i]) %>%
    dplyr::select(site,year,value) 
  tmp.nu = tmp.nu %>%
    data.frame() %>%
    tidyr::pivot_longer(cols=1:3) %>%
    dplyr::mutate(year = ifelse(name=="X2006",2006,
                                ifelse(name=="X2007",2007,2008))) %>%
    dplyr::mutate(site = siteNames[i]) %>%
    dplyr::select(site,year,value) 
  zc.g1.df[[i]]<-tmp.g1
  zc.s1.df[[i]]<-tmp.s1
  zc.s2.df[[i]]<-tmp.s2
  zc.nu.df[[i]]<-tmp.nu
  
}

zc.g1.df<-data.table::rbindlist(zc.g1.df)
zc.s1.df<-data.table::rbindlist(zc.s1.df)
zc.s2.df<-data.table::rbindlist(zc.s2.df)
zc.nu.df<-data.table::rbindlist(zc.nu.df)

# 
# med<-zc.df %>%
#   dplyr::group_by(site,year) %>%
#   dplyr::summarise(p.med = median(value))
# 
# df.med <- df %>%
#   dplyr::group_by(site,site) %>%
#   dplyr::summarise(p.mu = mean(seedlingJan)/100)
# 
# df=df %>%
#   dplyr::rename(year=site,site=site) %>%
#   dplyr::mutate(year=as.numeric(as.character(year)))
# 
# df.join<-df %>% dplyr::left_join(med,by=c('site','year'))
# 
# ggplot() +
#   geom_point(data=df.join,aes(x=p.med,y=seedlingJan/seedStart),
#               alpha=0.5) +
#   geom_point(data=df.join %>%
#                dplyr::group_by(site,year,p.med) %>%
#                dplyr::summarise(p.mu = mean(seedlingJan)/100),
#              aes(x=p.med,y=p.mu),size=2,color='red') +
#   geom_abline(intercept=0,slope=1) +
# #  xlim(c(0,1)) + ylim(c(0,1)) +
#   theme_bw() +
#   facet_wrap(~year,scales='free')
# 
# ggplot() +
#   geom_point(data=df.join,aes(x=p.med,y=seedlingJan/seedStart),
#              alpha=0.5) +
#   geom_point(data=df.join %>%
#                dplyr::group_by(site,year,p.med) %>%
#                dplyr::summarise(p.mu = mean(seedlingJan)/100),
#              aes(x=p.med,y=p.mu,color=as.factor(year)),size=2) +
#   geom_abline(intercept=0,slope=1) +
#   theme_bw() +
#   facet_wrap(~site,scales='free')
# 
# 
# zc.pool <- MCMCchains(mcmcSamples,params="g1")
# zc.pool.s1 <- MCMCchains(mcmcSamples,params="s1")
# zc.pool=zc.pool*zc.pool.s1
# colnames(zc.pool) = siteNames
# 
# zc.pool = zc.pool %>% 
#   data.frame() %>%
#   tidyr::pivot_longer(1:20) %>%
#   dplyr::rename(site=name)
# 
# ggplot(data=zc.pool) +
#   geom_boxplot(aes(x=site,y=value))
# 
# zc.pool.sum = zc.pool %>%
#   dplyr::group_by(site) %>%
#   dplyr::summarise(site.med = median(value))
# 
# ggplot() +
#   geom_hline(data=zc.pool.sum,aes(yintercept=site.med),linetype='dotted') +
#   geom_point(data=df.join,aes(x=p.med,y=seedlingJan/seedStart),
#              alpha=0.5) +
#   geom_point(data=df.join %>%
#                dplyr::group_by(site,year,p.med) %>%
#                dplyr::summarise(p.mu = mean(seedlingJan)/100),
#              aes(x=p.med,y=p.mu,color=as.factor(year)),size=2) +
#   geom_abline(intercept=0,slope=1) +
#   theme_bw() +
#   facet_wrap(~site,scales='free')
# 
# plot(zc.pool$site,zc.pool$site.med)
# ## What explains the high residuals in S22 2008 and GCN 2006?
# 
# residuals = df.join %>%
#   dplyr::group_by(site,year,p.med) %>%
#   dplyr::summarise(p.hat = mean(seedlingJan)/100) %>%
#   dplyr::mutate(resid = p.med-p.hat)
# 
# ## Check if these sites*years had few seed bags
# 
# residuals=df.join %>%
#   dplyr::group_by(site,year) %>%
#   dplyr::summarise(count = n()) %>%
#   dplyr::left_join(residuals,by=c('site','year'))
# 
# plot(residuals$resid,residuals$count)
# 
# ## check if these sites*years had high variance
# 
# residuals=df.join %>%
#   dplyr::group_by(site,year) %>%
#   dplyr::summarise(var = var(seedlingJan)) %>%
#   dplyr::left_join(residuals,by=c('site','year'))
# 
# plot(residuals$resid,residuals$var)
# 
# ## check if the site*years just had few seedlings?
# 
# residuals=df.join %>%
#   dplyr::group_by(site,year) %>%
#   dplyr::summarise(tot = sum(seedlingJan)) %>%
#   dplyr::left_join(residuals,by=c('site','year'))
# 
# plot(residuals$resid,residuals$tot)
# 
# # high resid because high number of seeds
# # check the site level estimate - seems like S22 2008 is same as prior?
# 
# # higher residuals for higher median estimates
# 
# plot(residuals$resid,residuals$p.med);abline(a=0,b=1)

################################################################################
# Load Am Nat estimates
#################################################################################
vitalRates <- readRDS(file="~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/eckhartEstimates.RDS")

# nu comparison
amnat.nu = vitalRates$PV1
amnat.nu=data.frame(amnat.nu) %>%
  bind_cols(site=siteNames) %>%
  dplyr::rename(`2006`=X1,`2007`=X2,`2008`=X3) %>%
  tidyr::pivot_longer(cols=1:3) %>%
  dplyr::mutate(year=as.numeric(name)) %>%
  dplyr::rename(amnat.est=value) %>%
  dplyr::select(-name)

zc.med<-zc.nu.df %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(p.med = median(value))

nu.join = zc.med %>%
  dplyr::left_join(amnat.nu,by=c('site','year'))

ggplot(data=nu.join) +
  geom_abline(intercept=0,slope=1) +
  geom_point(aes(x=amnat.est,p.med,color=as.factor(year)))+
  facet_wrap(~site)

ggplot(data=nu.join) +
  geom_abline(intercept=0,slope=1) +
  geom_point(aes(x=amnat.est,p.med,color=as.factor(year)))

nu.resid = nu.join %>%
  dplyr::mutate(residual.nu=abs(p.med-amnat.est)) %>%
  dplyr::select(site,year,residual.nu) 

nu.join %>% dplyr::left_join(nu.resid,by=c("site","year")) %>% View

viabFinal <- readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/viabilityFinal.RDS")

df<-viabFinal %>%
  dplyr::group_by(site,round) %>%
  dplyr::summarise(germStart = sum(germStart),viabStart=sum(viabStart)) %>%
  dplyr::mutate(year=round+2005) %>%
  dplyr::left_join(nu.resid,by=c("site","year"))

plot(df$germStart,df$residual.nu)
plot(df$viabStart+df$germStart,df$residual.nu)

df<-viabFinal %>%
  dplyr::group_by(site,round) %>%
  dplyr::summarise(n=n_distinct(bagNo)) %>%
  dplyr::mutate(year=round+2005) %>%
  dplyr::left_join(nu.resid,by=c("site","year"))

plot(df$n,df$residual.nu)

# g1 comparison
amnat.g1 = vitalRates$g1
amnat.g1=data.frame(amnat.g1) %>%
  bind_cols(site=siteNames) %>%
  dplyr::rename(`2006`=X1,`2007`=X2,`2008`=X3) %>%
  tidyr::pivot_longer(cols=1:3) %>%
  dplyr::mutate(year=as.numeric(name)) %>%
  dplyr::rename(amnat.est=value) %>%
  dplyr::select(-name)

zc.med<-zc.g1.df %>%
     dplyr::group_by(site,year) %>%
     dplyr::summarise(p.med = median(value))

g1.join = zc.med %>%
  dplyr::left_join(amnat.g1,by=c('site','year'))

ggplot(data=g1.join) +
  geom_abline(intercept=0,slope=1) +
  geom_point(aes(x=amnat.est,p.med,color=as.factor(year)))

ggplot(data=g1.join %>%
  dplyr::left_join(nu.resid ,by=c('site','year')))+
  geom_point(aes(x=residual.nu,y=p.med-amnat.est)) 

# s1 comparison
amnat.s1 = vitalRates$s1
amnat.s1=data.frame(amnat.s1) %>%
  bind_cols(site=siteNames) %>%
  dplyr::rename(`2006`=X1,`2007`=X2,`2008`=X3) %>%
  tidyr::pivot_longer(cols=1:3) %>%
  dplyr::mutate(year=as.numeric(name)) %>%
  dplyr::rename(amnat.est=value) %>%
  dplyr::select(-name)

zc.med<-zc.s1.df %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(p.med = median(value))

s1.join = zc.med %>%
  dplyr::left_join(amnat.s1,by=c('site','year'))

ggplot(data=s1.join) +
  geom_abline(intercept=0,slope=1) +
  geom_point(aes(x=amnat.est,p.med,color=as.factor(year)))

ggplot(data=s1.join %>%
         dplyr::left_join(nu.resid ,by=c('site','year')))+
  geom_point(aes(x=residual.nu,y=p.med-amnat.est))

# s2 comparison
amnat.s2 = vitalRates$s2
amnat.s2=data.frame(amnat.s2) %>%
  bind_cols(site=siteNames) %>%
  dplyr::rename(`2006`=X1,`2007`=X2,`2008`=X3) %>%
  tidyr::pivot_longer(cols=1:3) %>%
  dplyr::mutate(year=as.numeric(name)) %>%
  dplyr::rename(amnat.est=value) %>%
  dplyr::select(-name)

zc.med<-zc.s2.df %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(p.med = median(value))

s2.join = zc.med %>%
  dplyr::left_join(amnat.s2,by=c('site','year'))

ggplot(data=s2.join) +
  geom_abline(intercept=0,slope=1) +
  geom_point(aes(x=amnat.est,p.med,color=as.factor(year)))

ggplot(data=s2.join %>%
         dplyr::left_join(nu.resid ,by=c('site','year')))+
  geom_point(aes(x=residual.nu,y=p.med-amnat.est))


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

