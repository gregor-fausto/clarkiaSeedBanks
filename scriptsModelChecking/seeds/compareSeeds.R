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

mcmcSamples<-readRDS(file=paste0("/Users/Gregor/Dropbox/dataLibrary/posteriors/modelChecking/seedSamples.RDS"))

#directory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/"
#simFiles <- paste0(directory,list.files(directory))

data <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/viabilityData.RDS")

# Summarize intact seeds -------------------------------------------------------

# summarize mcmc samples for per-bag estimates of viability (round 2 of trials)
summary <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p0_1[siteBags,yearBags]) %>%
  dplyr::group_by(siteBags,yearBags) %>%
  dplyr::summarise(p.med = mean(p0_1)) 

# load viability data passed that gets passed to JAGS
seedBagsData<-readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/seedBagsRawData.RDS")

summary <- data.frame(site=unique(seedBagsData$siteBags),siteBags=1:20) %>%
  dplyr::left_join(summary,by=c('siteBags')) %>%
  dplyr::mutate(yearBags=yearBags+2005) %>%
  dplyr::select(-siteBags) %>%
  dplyr::rename(siteBags=site)

# organize and join data
seedEstimates <- seedBagsData %>%
  dplyr::mutate(yearBags=as.numeric(as.character(yearBags)) ) %>%
  dplyr::left_join(summary,by=c('siteBags','yearBags')) %>%
  dplyr::rename(site=siteBags,round=yearBags,age=ageBags) 

est<-seedEstimates %>%
  dplyr::mutate(p=totalJan/seedStart) %>%
  dplyr::group_by(site,round,p.med) %>%
  dplyr::summarise(p.mu = mean(p),count=n(),s=var(totalJan))

# plot observed vs. predicted estimates of per-bag viability (predicted is from fit)
ggplot() +
  geom_abline(slope=1,intercept=0) +
  geom_point(data=seedEstimates,aes(y=totalJan/seedStart,x=p.med),size=0.5,color="darkgray") +
  geom_point(data=est,aes(y=p.mu,x=p.med)) +
  xlab("Predicted proportion of seeds in January") + ylab("Observed proportion of seeds") +
  theme_bw() +
  xlim(c(0,1)) + ylim(c(0,1))

# More variable samples have more pooling and thus greater difference between the methods
ggplot(data=est) +
  geom_hline(yintercept=0,color='black') +
  geom_point(aes(y=p.med-p.mu,x=s),alpha=0.25) +
  xlab("Sample variance of seeds in January") + ylab("Residual of Bayesian-MLE estimates") +
  theme_bw()


# Summarize seedlings -------------------------------------------------------

# summarize mcmc samples for per-bag estimates of viability (round 2 of trials)
summary <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p0_2[siteBags,yearBags]) %>%
  dplyr::group_by(siteBags,yearBags) %>%
  dplyr::summarise(p.med = mean(p0_2)) 

# load viability data passed that gets passed to JAGS
seedBagsData<-readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/seedBagsRawData.RDS")

summary <- data.frame(site=unique(seedBagsData$siteBags),siteBags=1:20) %>%
  dplyr::left_join(summary,by=c('siteBags')) %>%
  dplyr::mutate(yearBags=yearBags+2005) %>%
  dplyr::select(-siteBags) %>%
  dplyr::rename(siteBags=site)

# organize and join data
seedEstimates <- seedBagsData %>%
  dplyr::mutate(yearBags=as.numeric(as.character(yearBags)) ) %>%
  dplyr::left_join(summary,by=c('siteBags','yearBags')) %>%
  dplyr::rename(site=siteBags,round=yearBags,age=ageBags) 

est<-seedEstimates %>%
  dplyr::mutate(p=seedlingJan/totalJan) %>%
  dplyr::group_by(site,round,p.med) %>%
  dplyr::summarise(p.mu = mean(p),count=n(),s=sum(totalJan))

# plot observed vs. predicted estimates of per-bag viability (predicted is from fit)
ggplot() +
  geom_abline(slope=1,intercept=0) +
  geom_point(data=seedEstimates,aes(y=seedlingJan/totalJan,x=p.med),size=0.5,color="darkgray") +
  geom_point(data=est,aes(y=p.mu,x=p.med)) +
  xlab("Predicted proportion of seedlings") + ylab("Observed proportion of seedlings") +
  theme_bw() +
  xlim(c(0,1)) + ylim(c(0,1))

# More variable samples have more pooling and thus greater difference between the methods
ggplot(data=est) +
  geom_hline(yintercept=0,color='black') +
  geom_point(aes(y=p.med-p.mu,x=s),alpha=0.25) +
  xlab("Sample size variance") + ylab("Residual of Bayesian-MLE estimates") +
  theme_bw()

# Summarize seedlings -------------------------------------------------------

# summarize mcmc samples for per-bag estimates of viability (round 2 of trials)
summary <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p0_4[siteBags,yearBags]) %>%
  dplyr::group_by(siteBags,yearBags) %>%
  dplyr::summarise(p.med = median(p0_4)) 

# load viability data passed that gets passed to JAGS
seedBagsData<-readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/seedBagsRawData.RDS")

summary <- data.frame(site=unique(seedBagsData$siteBags),siteBags=1:20) %>%
  dplyr::left_join(summary,by=c('siteBags')) %>%
  dplyr::mutate(yearBags=yearBags+2005) %>%
  dplyr::select(-siteBags) %>%
  dplyr::rename(siteBags=site)

# organize and join data
seedEstimates <- seedBagsData %>%
  dplyr::mutate(yearBags=as.numeric(as.character(yearBags)) ) %>%
  dplyr::left_join(summary,by=c('siteBags','yearBags')) %>%
  dplyr::rename(site=siteBags,round=yearBags,age=ageBags) 

est<-seedEstimates %>%
  dplyr::mutate(p=seedlingJan/seedStart) %>%
  dplyr::group_by(site,round,p.med) %>%
  dplyr::summarise(p.mu = median(p),count=n(),s=sum(totalJan))

# plot observed vs. predicted estimates of per-bag viability (predicted is from fit)
ggplot() +
  geom_abline(slope=1,intercept=0) +
  geom_point(data=seedEstimates,aes(y=seedlingJan/seedStart,x=p.med),size=0.5,color="darkgray") +
  geom_point(data=est,aes(y=p.mu,x=p.med)) +
  xlab("Predicted proportion of seedlings") + ylab("Observed proportion of seedlings") +
  theme_bw() +
  xlim(c(0,.3)) + ylim(c(0,.6))


ggplot(data=seedBagsData,aes(x=seedlingJan/seedStart)) +
  geom_density(aes(fill=yearBags),alpha=0.25) +
  geom_vline(data=est %>% dplyr::rename(siteBags=site,yearBags=round),
             aes(xintercept=p.med,color=as.factor(yearBags))) +
  facet_wrap(~siteBags,scales='free')


# Summarize October intact -------------------------------------------------------

# summarize mcmc samples for per-bag estimates of viability (round 2 of trials)
summary <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p0_3[siteBags,yearBags]) %>%
  dplyr::group_by(siteBags,yearBags) %>%
  dplyr::summarise(p.med = mean(p0_3)) 

# load viability data passed that gets passed to JAGS
seedBagsData<-readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/seedBagsRawData.RDS")

summary <- data.frame(site=unique(seedBagsData$siteBags),siteBags=1:20) %>%
  dplyr::left_join(summary,by=c('siteBags')) %>%
  dplyr::mutate(yearBags=yearBags+2005) %>%
  dplyr::select(-siteBags) %>%
  dplyr::rename(siteBags=site)

# organize and join data
seedEstimates <- seedBagsData %>%
  dplyr::mutate(yearBags=as.numeric(as.character(yearBags)) ) %>%
  dplyr::left_join(summary,by=c('siteBags','yearBags')) %>%
  dplyr::rename(site=siteBags,round=yearBags,age=ageBags) 

est<-seedEstimates %>%
  dplyr::mutate(p=intactOct/(totalJan-seedlingJan)) %>%
  dplyr::group_by(site,round,p.med) %>%
  dplyr::summarise(p.mu = mean(p),count=n(),s=sum((totalJan-seedlingJan)))

# plot observed vs. predicted estimates of per-bag viability (predicted is from fit)
ggplot() +
  geom_abline(slope=1,intercept=0) +
  geom_point(data=seedEstimates,aes(y=intactOct/(totalJan-seedlingJan),x=p.med),size=0.5,color="darkgray") +
  geom_point(data=est,aes(y=p.mu,x=p.med)) +
  xlab("Predicted proportion of seeds in October") + ylab("Observed proportion of seeds in October") +
  theme_bw() +
  xlim(c(0,1)) + ylim(c(0,1))


# Summarize October intact -------------------------------------------------------

# summarize mcmc samples for per-bag estimates of viability (round 2 of trials)
summary <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p0_5[siteBags,yearBags]) %>%
  dplyr::group_by(siteBags,yearBags) %>%
  dplyr::summarise(p.med = mean(p0_5)) 

# load viability data passed that gets passed to JAGS
seedBagsData<-readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/seedBagsRawData.RDS")

summary <- data.frame(site=unique(seedBagsData$siteBags),siteBags=1:20) %>%
  dplyr::left_join(summary,by=c('siteBags')) %>%
  dplyr::mutate(yearBags=yearBags+2005) %>%
  dplyr::select(-siteBags) %>%
  dplyr::rename(siteBags=site)

# organize and join data
seedEstimates <- seedBagsData %>%
  dplyr::mutate(yearBags=as.numeric(as.character(yearBags)) ) %>%
  dplyr::left_join(summary,by=c('siteBags','yearBags')) %>%
  dplyr::rename(site=siteBags,round=yearBags,age=ageBags) 

est<-seedEstimates %>%
  dplyr::mutate(p=intactOct/seedStart) %>%
  dplyr::group_by(site,round,p.med) %>%
  dplyr::summarise(p.mu = mean(p),count=n(),s=sum(seedStart))

# plot observed vs. predicted estimates of per-bag viability (predicted is from fit)
ggplot() +
  geom_abline(slope=1,intercept=0) +
  geom_point(data=seedEstimates,aes(y=intactOct/seedStart,x=p.med),size=0.5,color="darkgray") +
  geom_point(data=est,aes(y=p.mu,x=p.med)) +
  xlab("Predicted proportion of seeds in October") + ylab("Observed proportion of seeds in October") +
  theme_bw() +
  xlim(c(0,1)) + ylim(c(0,1))


# Comparing derived quantities --------------------------------------------

theta1=MCMCchains(mcmcSamples,params="theta_1")
theta2=MCMCchains(mcmcSamples,params="theta_2")
theta3=MCMCchains(mcmcSamples,params="theta_3")

theta4=MCMCchains(mcmcSamples,params="theta_4")
theta5=MCMCchains(mcmcSamples,params="theta_5")

plot(apply(theta1*theta2,2,median),apply(theta4,2,median))
abline(a=0,b=1)

plot(apply(theta1*(1-theta2)*theta3,2,mean),apply(theta5,2,mean))
abline(a=0,b=1)
     

p1=MCMCchains(mcmcSamples,params="p0_1")
p2=MCMCchains(mcmcSamples,params="p0_2")
p3=MCMCchains(mcmcSamples,params="p0_3")
p4=MCMCchains(mcmcSamples,params="p0_4")
p5=MCMCchains(mcmcSamples,params="p0_5")

plot(apply(p1*p2,2,median),apply(p4,2,median))
abline(a=0,b=1)

plot(apply(p1*(1-p2)*p3,2,median),apply(p5,2,median))
abline(a=0,b=1)
