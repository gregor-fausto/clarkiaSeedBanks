rm(list=ls(all=TRUE)) # clear R environment

# load packages

library(MCMCvis)
library(tidybayes)
library(tidyverse)
library(magrittr)
library(bayesplot)
library(rethinking)


# Summary -----------------------------------------------------------------

# Load MCMC samples -------------------------------------------------------

mcmcSamples<-readRDS(file=paste0("/Users/Gregor/Dropbox/dataLibrary/posteriors/modelChecking/seedSamplesAllAges.RDS"))

#directory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/"
#simFiles <- paste0(directory,list.files(directory))

data <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/viabilityData.RDS")

position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,northing,elevation) %>%
  dplyr::mutate(easting=easting/1000,northing=northing/1000)

# load seed bag data passed that gets passed to JAGS
seedBagsData<-readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/seedBagsRawData.RDS")

# Summarize per-bag intact seeds -------------------------------------------------------

# summarize mcmc samples for per-bag estimates of intact seeds
summary <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(theta_1[bag]) %>%
  dplyr::group_by(bag) %>%
  dplyr::summarise(p.med = median(theta_1)) 

# organize and join data
seedBagsEstimates <- seedBagsData %>%
  #dplyr::rename(id=bag ) %>%
  tibble::rowid_to_column( "bag") %>%
  dplyr::left_join(summary,by='bag') %>%
  dplyr::rename(site=siteBags,round=yearBags,age=ageBags) %>%
  dplyr::mutate(round=as.double(round),
                age=as.double(age)) 

# plot observed vs. predicted estimates of per-bag viability (predicted is from fit)
ggplot(data=seedBagsEstimates) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(y=totalJan/seedStart,x=p.med)) +
  xlab("Predicted per-bag intact to January") + ylab("Observed per-bag intact to January") +
  theme_bw()

# PLOT 
ggplot(data=seedBagsEstimates) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(y=totalJan/seedStart,x=p.med)) +
  xlab("Predicted per-bag intact to January") + ylab("Observed per-bag intact to January") +
  theme_bw() +
  facet_wrap(~age) 

# Summarize per-bag seedlings -------------------------------------------------------

# summarize mcmc samples for per-bag estimates of intact seeds
summary <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(theta_2[bag]) %>%
  dplyr::group_by(bag) %>%
  dplyr::summarise(p.med = median(theta_2)) 

# organize and join data
seedBagsEstimates <- seedBagsData %>%
  tibble::rowid_to_column( "bag") %>%
  dplyr::left_join(summary,by='bag') %>%
  dplyr::rename(site=siteBags,round=yearBags,age=ageBags) %>%
  dplyr::mutate(round=as.double(round),
                age=as.double(age)) 

# plot observed vs. predicted estimates of per-bag viability (predicted is from fit)
ggplot(data=seedBagsEstimates) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(y=seedlingJan/totalJan,x=p.med)) +
  xlab("Predicted per-bag emergence rate") + ylab("Observed per-bag emergence rate") +
  theme_bw()

# PLOT 
ggplot(data=seedBagsEstimates) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(y=seedlingJan/totalJan,x=p.med,size=totalJan),alpha=0.25) +
  xlab("Predicted per-bag emergence rate") + ylab("Observed per-bag emergence rate") +
  theme_bw() +
  scale_size_continuous(name="n Trials") +
  facet_wrap(~age) 

# Summarize per-bag January to October -------------------------------------------------------

# summarize mcmc samples for per-bag estimates of intact seeds
summary <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(theta_3[bag]) %>%
  dplyr::group_by(bag) %>%
  dplyr::summarise(p.med = median(theta_3)) 

# organize and join data
seedBagsEstimates <- seedBagsData %>%
  tibble::rowid_to_column( "bag") %>%
  dplyr::left_join(summary,by='bag') %>%
  dplyr::rename(site=siteBags,round=yearBags,age=ageBags) %>%
  dplyr::mutate(round=as.double(round),
                age=as.double(age)) 

# plot observed vs. predicted estimates of per-bag viability (predicted is from fit)
ggplot(data=seedBagsEstimates) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(y=intactOct/(totalJan-seedlingJan),x=p.med)) +
  xlab("Predicted per-bag survival rate") + ylab("Observed per-bag survival rate") +
  theme_bw()

# PLOT 
ggplot(data=seedBagsEstimates) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(y=intactOct/(totalJan-seedlingJan),x=p.med,size=(totalJan-seedlingJan)),alpha=0.25) +
  xlab("Predicted per-bag survival rate") + ylab("Observed per-bag survival rate") +
  theme_bw() +
  scale_size_continuous(name="n Trials") +
  facet_wrap(~age) 

# Summarize per-bag intact to October seeds -------------------------------------------------------

# summarize mcmc samples for per-bag estimates of intact seeds
summary <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(theta_4[bag]) %>%
  dplyr::group_by(bag) %>%
  dplyr::summarise(p.med = median(theta_4)) 

# organize and join data
seedBagsEstimates <- seedBagsData %>%
  tibble::rowid_to_column( "bag") %>%
  dplyr::left_join(summary,by='bag') %>%
  dplyr::rename(site=siteBags,round=yearBags,age=ageBags) %>%
  dplyr::mutate(round=as.double(round),
                age=as.double(age)) 

# plot observed vs. predicted estimates of per-bag viability (predicted is from fit)
ggplot(data=seedBagsEstimates) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(y=intactOct/seedStart,x=p.med)) +
  xlab("Predicted per-bag survival rate") + ylab("Observed per-bag survival rate") +
  theme_bw()

# PLOT 
ggplot(data=seedBagsEstimates) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(y=intactOct/(seedStart),x=p.med,size=(seedStart)),alpha=0.25) +
  xlab("Predicted per-bag survival rate") + ylab("Observed per-bag survival rate") +
  theme_bw() +
  scale_size_continuous(name="n Trials") +
  facet_wrap(~age) 

# Summarize intact seeds -------------------------------------------------------

# summarize mcmc samples for per-bag estimates of viability (round 2 of trials)
summary <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p0_1[siteBags,matIndex]) %>%
  dplyr::group_by(siteBags,matIndex) %>%
  dplyr::summarise(p.med = median(p0_1),
                   ci.lo = HPDI(p0_1,0.95)[1], 
                   ci.hi = HPDI(p0_1,0.95)[2],) 

# load viability data passed that gets passed to JAGS
seedBagsData<-readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/seedBagsRawData.RDS")

refTable <- data.frame(site=rep(unique(seedBagsData$siteBags),each=6),
                       siteBags=rep(1:20,each=6),
                       matIndex=rep(1:6,20),
                       round=rep(c(1,2,3,1,2,1),20),age=rep(c(1,1,1,2,2,3),20))

seedSummary <- summary %>%
  dplyr::left_join(refTable,by=c('siteBags','matIndex'))

ggplot(seedSummary) +
  geom_pointrange(aes(x=age,y=p.med,
                      ymin=ci.lo,ymax=ci.hi,
                      color=as.factor(round))) +
  ylim(c(0,1)) + theme_bw() + 
  facet_wrap(~site) +
  xlab("Age (years)") + ylab("Median predicted survival to January (95% HPDI)") +
  scale_color_discrete(name="Round")


# summarize mcmc samples for per-bag estimates of viability (round 2 of trials)
summary <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p_1[siteBags,ageBags]) %>%
  dplyr::group_by(siteBags,ageBags) %>%
  dplyr::summarise(p.med = median(p_1),
                   ci.lo = HPDI(p_1,0.95)[1], 
                   ci.hi = HPDI(p_1,0.95)[2],) 

refTable <- data.frame(site=rep(unique(seedBagsData$siteBags),each=6),
                       siteBags=rep(1:20,each=6),
                       matIndex=rep(1:6,20),
                       round=rep(c(1,2,3,1,2,1),20),
                       ageBags=rep(c(1,1,1,2,2,3),20))

seedSummary <- summary %>%
  dplyr::left_join(refTable,by=c('siteBags','ageBags'))

ggplot(seedSummary) +
  geom_pointrange(aes(x=ageBags,y=p.med,
                      ymin=ci.lo,ymax=ci.hi)) +
  ylim(c(0,1)) + theme_bw() + 
  facet_wrap(~site) +
  xlab("Age (years)") + ylab("Median predicted survival to January (95% HPDI)") +
  scale_color_discrete(name="Round")

seedSummary <- summary %>%
  dplyr::left_join(refTable,by=c('siteBags','ageBags')) %>%
  dplyr::left_join(position,by='site')

ggplot(seedSummary) +
  geom_pointrange(aes(x=easting,y=p.med,ymin=ci.lo,ymax=ci.hi,color=as.factor(ageBags)), alpha=0.5) +
  ylim(c(0,1)) + theme_bw() +   facet_wrap(~ageBags) +
  xlab("Easting") + ylab("Median predicted survival to January (95% HPDI)") +
  scale_color_discrete(name="Age")

# Summarize seedlings -------------------------------------------------------

# summarize mcmc samples for per-bag estimates of viability (round 2 of trials)
summary <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p0_2[siteBags,matIndex]) %>%
  dplyr::group_by(siteBags,matIndex) %>%
  dplyr::summarise(p.med = median(p0_2),
                   ci.lo = HPDI(p0_2,0.95)[1], 
                   ci.hi = HPDI(p0_2,0.95)[2],) 

# load viability data passed that gets passed to JAGS
seedBagsData<-readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/seedBagsRawData.RDS")

refTable <- data.frame(site=rep(unique(seedBagsData$siteBags),each=6),
                       siteBags=rep(1:20,each=6),
                       matIndex=rep(1:6,20),
                       round=rep(c(1,2,3,1,2,1),20),age=rep(c(1,1,1,2,2,3),20))

seedSummary <- summary %>%
  dplyr::left_join(refTable,by=c('siteBags','matIndex'))

ggplot(seedSummary) +
  geom_pointrange(aes(x=age,y=p.med,
                      ymin=ci.lo,ymax=ci.hi,
                      color=as.factor(round))) +
  ylim(c(0,1)) + theme_bw() + 
  facet_wrap(~site) +
  xlab("Age (years)") + ylab("Median predicted emergence (95% HPDI)") +
  scale_color_discrete(name="Round")

# summarize mcmc samples for per-bag estimates of viability (round 2 of trials)
summary <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p_2[siteBags,ageBags]) %>%
  dplyr::group_by(siteBags,ageBags) %>%
  dplyr::summarise(p.med = median(p_2),
                   ci.lo = HPDI(p_2,0.95)[1], 
                   ci.hi = HPDI(p_2,0.95)[2],) 

refTable <- data.frame(site=rep(unique(seedBagsData$siteBags),each=6),
                       siteBags=rep(1:20,each=6),
                       matIndex=rep(1:6,20),
                       round=rep(c(1,2,3,1,2,1),20),
                       ageBags=rep(c(1,1,1,2,2,3),20))

seedSummary <- summary %>%
  dplyr::left_join(refTable,by=c('siteBags','ageBags'))

ggplot(seedSummary) +
  geom_pointrange(aes(x=ageBags,y=p.med,
                      ymin=ci.lo,ymax=ci.hi)) +
  ylim(c(0,1)) + theme_bw() + 
  facet_wrap(~site) +
  xlab("Age (years)") + ylab("Median predicted emergence (95% HPDI)") +
  scale_color_discrete(name="Round")


seedSummary <- summary %>%
  dplyr::left_join(refTable,by=c('siteBags','ageBags')) %>%
  dplyr::left_join(position,by='site')

ggplot(seedSummary) +
  geom_pointrange(aes(x=easting,y=p.med,ymin=ci.lo,ymax=ci.hi,color=as.factor(ageBags)), alpha=0.5) +
  ylim(c(0,1)) + theme_bw() +   facet_wrap(~ageBags) +
  xlab("Easting") + ylab("Median predicted emergence (95% HPDI)") +
  scale_color_discrete(name="Age")

# Summarize survival to October -------------------------------------------------------

# summarize mcmc samples for per-bag estimates of viability (round 2 of trials)
summary <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p0_3[siteBags,matIndex]) %>%
  dplyr::group_by(siteBags,matIndex) %>%
  dplyr::summarise(p.med = median(p0_3),
                   ci.lo = HPDI(p0_3,0.95)[1], 
                   ci.hi = HPDI(p0_3,0.95)[2],) 

# load viability data passed that gets passed to JAGS
seedBagsData<-readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/seedBagsRawData.RDS")

sum<-seedBagsData %>% 
  dplyr::mutate(theta = intactOct/seedStart) %>%
  dplyr::group_by(siteBags,yearBags,ageBags) %>%
  dplyr::summarise(mu =mean(theta))

ggplot(sum) +
  geom_point(aes(x=ageBags,y=mu,
                      color=as.factor(yearBags))) +
  ylim(c(0,1)) + theme_bw() + 
  facet_wrap(~siteBags) +
  xlab("Age (years)") + ylab("Median predicted survival from January to October (95% HPDI)") +
  scale_color_discrete(name="Round")

refTable <- data.frame(site=rep(unique(seedBagsData$siteBags),each=6),
                       siteBags=rep(1:20,each=6),
                       matIndex=rep(1:6,20),
                       round=rep(c(1,2,3,1,2,1),20),age=rep(c(1,1,1,2,2,3),20))

seedSummary <- summary %>%
  dplyr::left_join(refTable,by=c('siteBags','matIndex'))

ggplot(seedSummary) +
  geom_pointrange(aes(x=age,y=p.med,
                      ymin=ci.lo,ymax=ci.hi,
                      color=as.factor(round))) +
  ylim(c(0,1)) + theme_bw() + 
  facet_wrap(~site) +
  xlab("Age (years)") + ylab("Median predicted survival from January to October (95% HPDI)") +
  scale_color_discrete(name="Round")

# summarize mcmc samples for per-bag estimates of viability (round 2 of trials)
summary <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p_3[siteBags,ageBags]) %>%
  dplyr::group_by(siteBags,ageBags) %>%
  dplyr::summarise(p.med = median(p_3),
                   ci.lo = HPDI(p_3,0.95)[1], 
                   ci.hi = HPDI(p_3,0.95)[2],) 

refTable <- data.frame(site=rep(unique(seedBagsData$siteBags),each=6),
                       siteBags=rep(1:20,each=6),
                       matIndex=rep(1:6,20),
                       round=rep(c(1,2,3,1,2,1),20),
                       ageBags=rep(c(1,1,1,2,2,3),20))

seedSummary <- summary %>%
  dplyr::left_join(refTable,by=c('siteBags','ageBags'))

ggplot(seedSummary) +
  geom_pointrange(aes(x=ageBags,y=p.med,
                      ymin=ci.lo,ymax=ci.hi)) +
  ylim(c(0,1)) + theme_bw() + 
  facet_wrap(~site) +
  xlab("Age (years)") + ylab("Median predicted survival from January to October (95% HPDI)") +
  scale_color_discrete(name="Round")


seedSummary <- summary %>%
  dplyr::left_join(refTable,by=c('siteBags','ageBags')) %>%
  dplyr::left_join(position,by='site')

ggplot(seedSummary) +
  geom_pointrange(aes(x=easting,y=p.med,ymin=ci.lo,ymax=ci.hi,color=as.factor(ageBags)), alpha=0.5) +
  ylim(c(0,1)) + theme_bw() +   facet_wrap(~ageBags) +
  xlab("Easting") + ylab("Median predicted survival from January to October (95% HPDI)") +
  scale_color_discrete(name="Age")


# Summarize survival to October -------------------------------------------------------

# summarize mcmc samples for per-bag estimates of viability (round 2 of trials)
summary <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p0_4[siteBags,matIndex]) %>%
  dplyr::group_by(siteBags,matIndex) %>%
  dplyr::summarise(p.med = median(p0_4),
                   ci.lo = HPDI(p0_4,0.95)[1], 
                   ci.hi = HPDI(p0_4,0.95)[2],) 

# load viability data passed that gets passed to JAGS
seedBagsData<-readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/seedBagsRawData.RDS")

refTable <- data.frame(site=rep(unique(seedBagsData$siteBags),each=6),
                       siteBags=rep(1:20,each=6),
                       matIndex=rep(1:6,20),
                       round=rep(c(1,2,3,1,2,1),20),age=rep(c(1,1,1,2,2,3),20))

seedSummary <- summary %>%
  dplyr::left_join(refTable,by=c('siteBags','matIndex'))

ggplot(seedSummary) +
  geom_pointrange(aes(x=age,y=p.med,
                      ymin=ci.lo,ymax=ci.hi,
                      color=as.factor(round))) +
  ylim(c(0,1)) + theme_bw() + 
  facet_wrap(~site) +
  xlab("Age (years)") + ylab("Median predicted survival from start to October (95% HPDI)") +
  scale_color_discrete(name="Round")

# summarize mcmc samples for per-bag estimates of viability (round 2 of trials)
summary <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p_4[siteBags,ageBags]) %>%
  dplyr::group_by(siteBags,ageBags) %>%
  dplyr::summarise(p.med = median(p_4),
                   ci.lo = HPDI(p_4,0.95)[1], 
                   ci.hi = HPDI(p_4,0.95)[2],) 

refTable <- data.frame(site=rep(unique(seedBagsData$siteBags),each=6),
                       siteBags=rep(1:20,each=6),
                       matIndex=rep(1:6,20),
                       round=rep(c(1,2,3,1,2,1),20),
                       ageBags=rep(c(1,1,1,2,2,3),20))

seedSummary <- summary %>%
  dplyr::left_join(refTable,by=c('siteBags','ageBags'))

ggplot(seedSummary) +
  geom_pointrange(aes(x=ageBags,y=p.med,
                      ymin=ci.lo,ymax=ci.hi)) +
  ylim(c(0,1)) + theme_bw() + 
  facet_wrap(~site) +
  xlab("Age (years)") + ylab("Median predicted survival from start to October (95% HPDI)") +
  scale_color_discrete(name="Round")


seedSummary <- summary %>%
  dplyr::left_join(refTable,by=c('siteBags','ageBags')) %>%
  dplyr::left_join(position,by='site')

ggplot(seedSummary) +
  geom_pointrange(aes(x=easting,y=p.med,ymin=ci.lo,ymax=ci.hi,color=as.factor(ageBags)), alpha=0.5) +
  ylim(c(0,1)) + theme_bw() +   facet_wrap(~ageBags) +
  xlab("Easting") + ylab("Median predicted survival from start to October (95% HPDI)") +
  scale_color_discrete(name="Age")

# Composite comparison ---------------------------------------------------------------
theta_1 = MCMCchains(mcmcSamples,params="theta_1")
theta_2 = MCMCchains(mcmcSamples,params="theta_2")
theta_3 = MCMCchains(mcmcSamples,params="theta_3")
theta_4 = MCMCchains(mcmcSamples,params="theta_4")

theta_comp = theta_1*(1-theta_2)*theta_3
par(mfrow=c(1,1))
plot(apply(theta_4,2,median),apply(theta_comp,2,median))

p0_4 = MCMCchains(mcmcSamples,params="p0_4")

p0_1 = MCMCchains(mcmcSamples,params="p0_1")

s3 = apply(p0_1[,61:80],2,median)/apply(p0_4[,1:20],2,median)

par(mfrow=c(1,1))
plot(s3)

# Composite ---------------------------------------------------------------

p_1 = MCMCchains(mcmcSamples,params="p_1")
p_2 = MCMCchains(mcmcSamples,params="p_2")
p_3 = MCMCchains(mcmcSamples,params="p_3")
p_4 = MCMCchains(mcmcSamples,params="p_4")

composite0.1 = p_1[,1:20]
composite1 = p_1[,1:20]*(1-p_2[,1:20])*p_3[,1:20]
p_4.1 = p_4[,1:20]
composite0.2 = p_1[,21:40]
composite2 = p_1[,21:40]*(1-p_2[,21:40])*p_3[,21:40]
composite0.3 = p_1[,41:60]
composite3 = p_1[,41:60]*(1-p_2[,41:60])*p_3[,41:60]


summaryComposite<- function(x,ageValue=1){
  
  tmp<-data.frame(siteBags=1:20,
                  age = rep(ageValue,20),
                            p.med=apply(x,2,median),
                            ci.lo=apply(x,2,HPDI,prob=.95)[1,],
                            ci.hi=apply(x,2,HPDI,prob=.95)[2,])
  return(tmp)
}

# composite<-summaryComposite(composite0.1,ageValue=1/3) %>%
#   bind_rows(summaryComposite(composite1,ageValue=1) ) %>%
#   bind_rows(summaryComposite(composite0.2,ageValue=(1+1/3)) ) %>%
#   bind_rows(summaryComposite(composite2,ageValue=2) ) %>%
#   bind_rows(summaryComposite(composite0.3,ageValue=(2+1/3)) ) %>%
#   bind_rows(summaryComposite(composite3,ageValue=3) )


composite<-summaryComposite(composite1,ageValue=1)  %>%
  bind_rows(summaryComposite(composite2,ageValue=2) ) %>%
  bind_rows(summaryComposite(composite3,ageValue=3) ) %>%
  bind_rows(summaryComposite(p_4.1,ageValue=1.1) )

refTable <- data.frame(site=rep(unique(seedBagsData$siteBags),each=6),
                       siteBags=rep(1:20,each=6)) %>% unique

seedSummary <- composite %>%
  dplyr::left_join(refTable,by=c('siteBags'))

ggplot(seedSummary) +
  geom_pointrange(aes(x=age,y=p.med,ymin=ci.lo,ymax=ci.hi)) +
 # geom_smooth(aes(x=age,y=p.med)) +
  ylim(c(0,1)) + theme_bw() + 
  xlim(c(0,4)) +
  facet_wrap(~site) +
  xlab("Age (years)") + ylab("Median predicted intact seeds in October (95% HPDI)") 

ggplot(seedSummary) +
 # geom_pointrange(aes(x=age,y=p.med,ymin=ci.lo,ymax=ci.hi)) +
  geom_smooth(aes(x=age,y=p.med,color=site)) +
  ylim(c(0,1)) + theme_bw() + 
  xlim(c(0,4)) +
 # facet_wrap(~site) +
  xlab("Age (years)") + ylab("Median predicted intact seeds in October (95% HPDI)") 



# In between ---------------------------------------------------------------

p_1 = MCMCchains(mcmcSamples,params="p0_1")
p_2 = MCMCchains(mcmcSamples,params="p0_2")
p_3 = MCMCchains(mcmcSamples,params="p0_3")
p_4 = MCMCchains(mcmcSamples,params="p0_4")


composite0.1 = p_1[,1:20]
composite1 = p_1[,1:20]*(1-p_2[,1:20])*p_3[,1:20]
p_4.1 = p_4[,1:20]
composite0.2 = p_1[,61:80]
composite2 = p_1[,61:80]*(1-p_2[,61:80])*p_3[,61:80]
composite0.3 = p_1[,101:120]
composite3 = p_1[,101:120]*(1-p_2[,101:120])*p_3[,101:120]

summaryComposite<- function(x,ageValue=1){
  
  tmp<-data.frame(siteBags=1:20,
                  age = rep(ageValue,20),
                  p.med=apply(x,2,median),
                  ci.lo=apply(x,2,HPDI,prob=.95)[1,],
                  ci.hi=apply(x,2,HPDI,prob=.95)[2,])
  return(tmp)
}

composite<-summaryComposite(composite0.1,ageValue=1/3) %>%
  bind_rows(summaryComposite(composite1,ageValue=1) ) %>%
  bind_rows(summaryComposite(composite0.2,ageValue=(1+1/3)) ) %>%
  bind_rows(summaryComposite(composite2,ageValue=2) ) %>%
  bind_rows(summaryComposite(composite0.3,ageValue=(2+1/3)) ) %>%
  bind_rows(summaryComposite(composite3,ageValue=3) )# %>%
  #bind_rows(summaryComposite(p_4.1,ageValue=1.1) )

# composite<- summaryComposite(composite1,ageValue=1)  %>%
#   bind_rows(summaryComposite(composite0.2,ageValue=(1+1/3)) ) %>%
#   bind_rows(summaryComposite(composite2,ageValue=2) ) %>%
#   bind_rows(summaryComposite(composite3,ageValue=3) ) %>%
#   bind_rows(summaryComposite(p_4.1,ageValue=1.1) )



refTable <- data.frame(site=rep(unique(seedBagsData$siteBags),each=6),
                       siteBags=rep(1:20,each=6)) %>% unique

seedSummary <- composite %>%
  dplyr::left_join(refTable,by=c('siteBags'))

ggplot(seedSummary) +
  geom_pointrange(aes(x=age,y=p.med,ymin=ci.lo,ymax=ci.hi)) +
  ylim(c(0,1)) + theme_bw() + xlim(c(0,4)) +
  facet_wrap(~site) +
  xlab("Age (years)") + ylab("Median predicted intact seeds in October (95% HPDI)") 


################################################################################
# Spatial data
#################################################################################
position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,northing,elevation) %>%
  dplyr::mutate(easting=easting/1000,northing=northing/1000)

seedSummary <- seedSummary %>%
  dplyr::left_join(position,by='site')


ggplot(seedSummary) +
  geom_pointrange(aes(x=easting,y=p.med,ymin=ci.lo,ymax=ci.hi,color=as.factor(age)), alpha=0.5) +
  ylim(c(0,1)) + theme_bw() + 
  facet_wrap(~age) +
  xlab("Age (years)") + ylab("Median predicted intact seeds in October (95% HPDI)") 
