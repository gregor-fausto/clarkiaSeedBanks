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

mcmcSamples<-readRDS(file=paste0("/Users/Gregor/Dropbox/dataLibrary/posteriors/modelChecking/viabilitySamplesAllAges.RDS"))

#directory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/"
#simFiles <- paste0(directory,list.files(directory))

data <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/viabilityData.RDS")
siteNames = unique(data$siteViab)

position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,northing,elevation) %>%
  dplyr::mutate(easting=easting/1000,northing=northing/1000)

# Summarize per-bag viability -------------------------------------------------------

# summarize mcmc samples for per-bag estimates of viability (round 2 of trials)
summary <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(theta_v[bag]) %>%
  dplyr::group_by(bag) %>%
  dplyr::summarise(p.med = median(theta_v)) 

# load viability data passed that gets passed to JAGS
viabilityRawData<-readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/viabilityRawData.RDS")

# organize and join data
viabilityEstimates <- viabilityRawData %>%
  dplyr::rename(id=bag ) %>%
  tibble::rowid_to_column( "bag") %>%
  dplyr::left_join(summary,by='bag') %>%
  dplyr::rename(site=siteViab,round=yearViab,age=ageViab) %>%
  dplyr::mutate(round=as.double(round),
                age=as.double(age)) 

# plot observed vs. predicted estimates of per-bag viability (predicted is from fit)
ggplot(data=viabilityEstimates) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(y=viabStain/viabStart,x=p.med)) +
  xlab("Predicted per-bag viability") + ylab("Observed per-bag viability (stained/started)") +
  theme_bw()

ggplot(data=viabilityEstimates) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(y=viabStain/viabStart,x=p.med)) +
  xlab("Predicted per-bag viability") + ylab("Observed per-bag viability (stained/started)") +
  theme_bw() +
  facet_wrap(~age) 

# comparison of observed data vs. values predicted by the statistical model
# points scaled to sample size (need different scaling?)
ggplot(data=viabilityEstimates) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(y=viabStain/viabStart,x=p.med,size=viabStart*.5),alpha=0.25) +
  xlab("Predicted per-bag viability") + ylab("Observed per-bag viability (stained/started)") +
  theme_bw()

# PLOT 
ggplot(data=viabilityEstimates) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(y=viabStain/viabStart,x=p.med,size=viabStart),alpha=0.25) +
  xlab("Predicted per-bag viability") + ylab("Observed per-bag viability (stained/started)") +
  theme_bw() + scale_size_continuous(name = "n Trials") +
  facet_wrap(~age) 

# # load estimates from Am Nat paper method
# estimatesAmNat = readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/matList.RDS")
# 
# # organize data
# estimatesAmNat <- estimatesAmNat %>%
#   tidyr::unite(col='id', c(site,bagNo,round,age), sep="-", remove=FALSE)
# 
# estimatesAmNat <- estimatesAmNat %>% 
#   dplyr::rename(bag=bagNo) %>%
#   dplyr::mutate(bag = as.double(bag)) %>%
#   dplyr::mutate(round=as.double(round)) %>%
#   dplyr::mutate(pv = as.numeric(pv),
#                 PVtot = as.numeric(PVtot),
#                 mult = as.numeric(mult),
#                 age = as.numeric(age),
#                 ns = as.numeric(ns),
#                 vs = as.numeric(vs)) %>%
#   dplyr::filter(!is.na(site))
# 
# # join data frames
# viabilityComparison = viabilityEstimates %>%
#   dplyr::left_join(estimatesAmNat,by=c("id",'site','round','age') )
# 
# # plot comparing 
# ggplot(data=viabilityComparison) +
#   geom_abline(slope=1,intercept=0) +
#   geom_point(aes(y=pv,x=p.med)) +
#   xlab("Predicted per-bag viability (Bayesian)") + ylab("Predicted per-bag viability (frequentist)") +
#   theme_bw()
# 
# # Small samples have more pooling and thus greater difference between the methods
# ggplot(data=viabilityComparison) +
#   geom_hline(yintercept=0,color='black') +
#   geom_point(aes(y=p.med-pv,x=viabStart),alpha=0.25) +
#   xlab("Number of trials (seeds starting viability trials)") + ylab("Residual of Bayesian-frequentist models") +
#   theme_bw()

# Summarize per-bag germination -------------------------------------------------------

# summarize mcmc samples for per-bag estimates of viability (round 1 of trials)
summary <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(theta_g[bag]) %>%
  dplyr::group_by(bag) %>%
  dplyr::summarise(p.med = median(theta_g)) 

# organize and join data
viabilityEstimates <- viabilityRawData %>%
  dplyr::rename(id=bag ) %>%
  tibble::rowid_to_column( "bag") %>%
  dplyr::left_join(summary,by='bag') %>%
  dplyr::rename(site=siteViab,round=yearViab,age=ageViab) %>%
  dplyr::mutate(round=as.double(round),
                age=as.double(age)) 

# plot observed vs. predicted estimates of per-bag viability (predicted is from fit)
ggplot(data=viabilityEstimates) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(y=germCount/germStart,x=p.med)) +
  xlab("Predicted per-bag germination") + ylab("Observed per-bag germination (count/start)") +
  theme_bw()

# plot observed vs. predicted estimates of per-bag viability (predicted is from fit)
ggplot(data=viabilityEstimates) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(y=germCount/germStart,x=p.med)) +
  xlab("Predicted per-bag germination") + ylab("Observed per-bag germination (count/start)") +
  theme_bw() +
  facet_wrap(~age)

# comparison of observed data vs. values predicted by the statistical model
# points scaled to sample size (need different scaling?)
ggplot(data=viabilityEstimates) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(y=germCount/germStart,x=p.med,size=germStart*.5),alpha=0.25) +
  xlab("Predicted per-bag germination") + ylab("Observed per-bag germination (count/start)") +
  theme_bw()

ggplot(data=viabilityEstimates) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(y=germCount/germStart,x=p.med,size=germStart*.5),alpha=0.25) +
  xlab("Predicted per-bag germination") + ylab("Observed per-bag germination (count/start)") +
  theme_bw() + scale_size_continuous(name = "n Trials") +
  facet_wrap(~age) 


# Overall viability over time -------------------------------------------------------

# summarize mcmc samples for per-bag estimates of total viability
summary <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(nu0_1[siteViab,matIndex]) %>%
  dplyr::group_by(siteViab,matIndex) %>%
  dplyr::summarise(p.med = median(nu0_1),
                   ci.lo = HPDI(nu0_1,0.95)[1], 
                   ci.hi = HPDI(nu0_1,0.95)[2],
                   ci.lo2 = HPDI(nu0_1,0.5)[1], 
                   ci.hi2 = HPDI(nu0_1,0.5)[2]) 

refTable <- data.frame(site=rep(unique(viabilityRawData$siteViab),each=6),
                       siteViab=rep(1:20,each=6),
                       matIndex=rep(1:6,20),
                       round=rep(c(1,2,3,1,2,1),20),age=rep(c(1,1,1,2,2,3),20))

viabilitySummary <- summary %>%
  dplyr::left_join(refTable,by=c('siteViab','matIndex'))

# PLOT
ggplot(viabilitySummary) +
  geom_pointrange(aes(x=age,y=p.med,
                      ymin=ci.lo,ymax=ci.hi,
                      color=as.factor(round)),alpha=0.5) +
  ylim(c(0,1)) + theme_bw() +facet_wrap(~site) +
  xlab("Age (years)") + ylab("Median predicted viability (95% HPDI)") +
  scale_color_discrete(name="Round")

# summarize mcmc samples for per-bag estimates of total viability
summary <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(nu_1[siteViab,age]) %>%
  dplyr::group_by(siteViab,age) %>%
  dplyr::summarise(p.med = median(nu_1),
                   ci.lo = HPDI(nu_1,0.95)[1], 
                   ci.hi = HPDI(nu_1,0.95)[2],
                   ci.lo2 = HPDI(nu_1,0.5)[1], 
                   ci.hi2 = HPDI(nu_1,0.5)[2]) 

refTable <- data.frame(site=rep(unique(viabilityRawData$siteViab),1),
                       siteViab=rep(1:20,1))

viabilitySummary <- summary %>%
  dplyr::left_join(refTable,by=c('siteViab'))

ggplot(viabilitySummary) +
  geom_pointrange(aes(x=age,y=p.med,
                      ymin=ci.lo,ymax=ci.hi)) +
  # geom_linerange(aes(x=age,
  #                    ymin=ci.lo2,ymax=ci.hi2),size=1.25) +
  ylim(c(0,1)) + theme_bw() +facet_wrap(~site) +
  xlab("Age (years)") + ylab("Median predicted viability (95% HPDI)") 


viabilitySummary <- viabilitySummary %>%
  dplyr::left_join(position,by='site')

ggplot(viabilitySummary) +
  geom_pointrange(aes(x=easting,y=p.med,ymin=ci.lo,ymax=ci.hi,color=as.factor(age)), alpha=0.5) +
  ylim(c(0,1)) + theme_bw() +   facet_wrap(~age) +
  xlab("Easting") + ylab("Median predicted survival viability (95% HPDI)") +
  scale_color_discrete(name="Age")


# ggplot(data=viabilityRawData %>%
#          dplyr::group_by(age,round) %>%
#          dplyr::summarise(mu.g = mean(germCount/germStart),
#                           mu.v = mean(viabStain/viabStart)) %>%
#          tidyr::pivot_longer(cols=3:4)) +
#   geom_point(aes(x=age,y=value)) +
#   ylim(c(0,1)) +
#   facet_wrap(name~round)
# 
# ggplot(data=viabilityRawData %>%
#          dplyr::group_by(age,round) %>%
#          dplyr::summarise(mu.g = mean(germCount/germStart),
#                           mu.v = mean(viabStain/viabStart)) %>%
#          dplyr::mutate(nu = mu.g + mu.v*(1-mu.g)) %>%
#          tidyr::pivot_longer(cols=3:5)) +
#   geom_point(aes(x=round,y=value)) +
#   ylim(c(0,1)) +
#   facet_wrap(name~age)
# 


# Overall viability *per bag* -------------------------------------------------------

# # summarize mcmc samples for per-bag estimates of total viability
# summary <- mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(nu_bag[bag]) %>%
#   dplyr::group_by(bag) %>%
#   dplyr::summarise(p.med = median(nu_bag)) 
# 
# # organize and join data
# viabilityEstimates <- viabilityRawData %>%
#   dplyr::rename(id=bag ) %>%
#   tibble::rowid_to_column( "bag") %>%
#   dplyr::left_join(summary,by='bag') %>%
#   dplyr::rename(site=siteViab,round=yearViab,age=ageViab) %>%
#   dplyr::mutate(round=as.double(round),
#                 age=as.double(age)) 
# 
# # join to am nat estimates
# # viabilityComparison = viabilityEstimates %>%
# #   dplyr::left_join(estimatesAmNat,by=c("id",'site','round','age') )
# 
# # comparison of predictions from both statistical models
# ggplot(data=viabilityEstimates) +
#   geom_abline(slope=1,intercept=0) +
#   geom_point(aes(x=p.med,y=PVtot)) +
#   xlab("Predicted composite viability") + ylab("Predicted PVtot viability") +
#   theme_bw()
# 
# # comparison of predictions from both statistical models
# # accounting for sample size (need better sample size)
# ggplot(data=viabilityComparison) +
#   geom_abline(slope=1,intercept=0,color='red') +
#   geom_point(aes(y=PVtot,x=p.med,size=viabStart+germStart),alpha=0.25) +
#   xlab("Predicted composite viability") + ylab("Predicted PVtot viability") +
#   theme_bw() 
# 
# # Small samples have more pooling and thus greater difference between the methods
# ggplot(data=viabilityComparison) +
#   geom_hline(yintercept=0,color='black') +
#   geom_point(aes(y=p.med-PVtot,x=viabStart+germStart),alpha=0.25) +
#   xlab("Number of trials (seeds starting germ+viability trials)") + ylab("Residual of Bayesian-frequentist estimates") +
#   theme_bw()
# 
# # Effect of pooling is strongest for values at the extremes
# # ie bags with estimated viability of 0 or 1 in the frequentist data
# # as these are likely to be pulled towards the year/site mean
# ggplot(data=viabilityComparison) +
#   geom_hline(yintercept=0,color='red') +
#   geom_point(aes(x=pv,y=abs(p.med-PVtot)),alpha=0.25) +
#   geom_smooth(aes(x=pv,y=abs(p.med-PVtot))) +
#   xlab("Per-bag viability estimate (MLE)") + ylab("Absolute value of residual of \nBayesian-frequentist estimates") +
#   theme_bw()

# Overall viability site*year -------------------------------------------------------

# summarize mcmc samples for pop*year estimates of viability
# note summarizing median and variance of estimates here
summary <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(nu0_1[siteViab,matIndex]) %>%
  dplyr::group_by(siteViab,matIndex) %>%
  dplyr::summarise(p.med = median(nu0_1),
                   p.var = var(nu0_1))  %>%
  dplyr::rename(site=siteViab,matIndex=matIndex)

refTable <- data.frame(matIndex=1:6,round=c(1,2,3,1,2,1),age=c(1,1,1,2,2,3))

viabilitySummary <- summary %>%
  dplyr::left_join(refTable,by=c('matIndex'))

# summarize mcmc samples for pop esimates of viability
summarySite <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(nu_1[siteViab,ageViab]) %>%
  dplyr::group_by(siteViab,ageViab) %>%
  dplyr::summarise(p.med = median(nu_1))  %>%
  dplyr::rename(site=siteViab,age=ageViab)

# organize viability data from JAGS fit
viabilityRawData <- viabilityRawData %>%
  dplyr::rename(id=bag ) %>%
  tibble::rowid_to_column( "bag") %>%
  dplyr::rename(site=siteViab,round=yearViab,age=ageViab) %>%
  dplyr::mutate(round=as.double(round),
                age=as.double(age)) 

# create reference table
refTable <- data.frame(siteName=unique(viabilityRawData$site),site=1)

# join site summary to reference table
viabilityEstimatesSite <- summarySite %>%
  dplyr::left_join(refTable,by="site") %>%
  dplyr::ungroup() %>%
  dplyr::select(-site) %>%
  dplyr::rename(site=siteName) %>%
  dplyr::rename(est_site=p.med) %>%
  dplyr::mutate(site=as.character(site))

# join site*year summary to reference table
viabilityEstimates <- viabilitySummary %>%
  dplyr::left_join(refTable,by="site") %>%
  dplyr::ungroup() %>%
  dplyr::select(-site) %>%
  dplyr::rename(site=siteName) %>%
#  dplyr::mutate(round = as.numeric(year)) %>%
  dplyr::mutate(site=as.character(site))

viabilityComp<-viabilityRawData %>%
  dplyr::left_join(viabilityEstimates,by=c('site','round','age'))

# comparison of predictions from both statistical models
ggplot(data=viabilityComp) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(x=p.med,y=viabStain/viabStart)) + xlim(c(0,1)) + ylim(c(0,1)) +
  xlab("Predicted composite viability") + ylab("Predicted Am Nat viability") +
  theme_bw() +
  facet_wrap(~round)


# load estimates viability for pop*year level from am nat paper
PV1 = readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/pv1List.RDS")

PV1 <- PV1 %>%
  dplyr::mutate(round = as.numeric(round),
                age = as.numeric(age),
                PVmean = as.numeric(PVmean),
                PVwts = as.numeric(PVwts),
                PVratio = as.numeric(PVratio))

# join am nat estimates to JAGS estimates
viabilityPopYear<-viabilityEstimates %>%
  dplyr::left_join(PV1,c("site","round")) 


# comparison of predictions from both statistical models
ggplot(data=viabilityPopYear) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(x=p.med,y=PVratio)) + xlim(c(0,1)) + ylim(c(0,1)) +
  xlab("Predicted composite viability") + ylab("Predicted Am Nat viability") +
  theme_bw()

ggplot(data=viabilityPopYear) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(x=p.med,y=PVratio,size=PVwts,col=as.factor(round)),alpha=0.5) + xlim(c(0,1)) + ylim(c(0,1)) +
  xlab("Predicted composite viability") + ylab("Predicted Am Nat viability") +
  theme_bw()

# Plot variance of estimate against difference between model estimates
ggplot(data=viabilityPopYear) +
  geom_point(aes(y=abs(p.med-PVratio),x=p.var),alpha=0.25) +
  xlab("Variance in composite viability") + ylab("Absolute value of residual of \nBayesian-frequentist estimates") +
  theme_bw()

# Look at variance in the Am Nat estimates
# PV tot is the composite viability per bag
varAmNat <- estimatesAmNat %>%
  dplyr::group_by(site,round) %>%
  dplyr::summarise(var = var(PVtot,na.rm=TRUE),
                   mu = mean(PVtot,na.rm=TRUE))

# join am nat estimates to JAGS estimates
viabilityPopYearAmNat<-viabilityEstimates %>%
  dplyr::left_join(varAmNat,c("site","round")) %>%
  dplyr::left_join(PV1,c("site","round")) 

# Plot variance of PVtot against difference between model estimates
ggplot(data=viabilityPopYearAmNat) +
  geom_point(aes(y=abs(p.med-PVratio),x=var),alpha=0.25) +
  xlab("Sample variance in PVtot") + ylab("Absolute value of residual of \nBayesian-frequentist estimates") +
  theme_bw()

# variance in PVtot positively correlated with variance in composite estimate
ggplot(data=viabilityPopYearAmNat) +
  geom_point(aes(y=var,x=p.var),alpha=0.25) +
  ylab("Variance in PVtot") + xlab("Variance in composite estimate of viability") +
  theme_bw()


ggplot(data=viabilityPopYearAmNat) +
  geom_hline(data=viabilityEstimatesSite,aes(yintercept=est_site),col='red') +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(y=p.med,x=PVratio,size=PVwts/100,alpha=var)) + 
  ylab("Predicted composite viability") + xlab("Predicted Am Nat viability") +
  theme_bw() + facet_wrap(~site,scales='free')

# visualizing shrinkage
shrinkage<-viabilityPopYearAmNat %>%
  dplyr::select(site,round,p.med,PVratio) %>%
  dplyr::rename(est1_freq=PVratio,est2_bayes=p.med) %>%
  tidyr::pivot_longer(cols=3:4) %>%
  dplyr::filter(!is.na(site)) 

ggplot(data=shrinkage) +
  geom_hline(data=viabilityEstimatesSite,aes(yintercept=est_site),col='red') +
  geom_point(aes(x=name,y=value)) +
  geom_line(aes(x=name,y=value,group=round,color=as.factor(round)) )+
  facet_wrap(~site) +
  theme_bw()

ggplot(estimatesAmNat,aes(x=gs,y=PVtot,group=round)) + geom_point() + 
  geom_smooth(method='lm') +facet_wrap(~site,scales='free')

# variance in PVtot positively correlated with variance in composite estimate
estimatesAmNat %>% 
  dplyr::filter(site=="SM" & round==3) %>%
  ggplot() + geom_histogram(aes(x=PVtot)) +
  theme_bw() + xlim(c(0,1))

# looking at another issue
viabilityRawData1 %>% 
  dplyr::group_by(siteViab,yearViab) %>%
  dplyr::summarise(s=sum(germStart)) %>% 
  dplyr::rename(site=siteViab,round=yearViab) %>%
  dplyr::mutate(round=as.double(as.character(round))) %>% 
  dplyr::left_join(viabilityPopYear,by=c('site','round')) %>% 
  dplyr::mutate(del=s-PVwts) %>% 
  dplyr::filter(del!=0) %>% View



viabilityPopYearVsRaw <- viabilityPopYear %>%
  dplyr::left_join(viabilityRawData1 %>%
                     dplyr::rename(site=siteViab,
                                   round=yearViab,
                                   age=ageViab) %>%
                     dplyr::mutate(round=as.double(as.character(round))),by=c("site","round"))

ggplot(data=viabilityPopYearVsRaw) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(x=p.med,y=viabStain/viabStart)) +
  theme_bw()

ggplot(data=viabilityPopYearVsRaw) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(x=PVratio,y=viabStain/viabStart)) +
  theme_bw()

ggplot(data=viabilityPopYearVsRaw) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(x=p.med,y=germCount/germStart)) +
  theme_bw()

ggplot(data=viabilityPopYearVsRaw) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(x=PVratio,y=germCount/germStart)) +
  theme_bw()

ggplot(data=viabilityPopYearVsRaw) +
  geom_abline(slope=1,intercept=0) +
  geom_point(aes(x=PVratio,y=p.med)) +
  theme_bw()


estimatesAmNat = readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/matList.RDS")

estimatesAmNat<-estimatesAmNat %>%
  dplyr::mutate(gs=as.numeric(gs),
                ns = as.numeric(ns),
                vs = as.numeric(vs),
                ng = as.numeric(ng)) %>%
  dplyr::mutate(round=as.numeric(round)) %>%
  dplyr::mutate(PVtot = as.numeric(PVtot))%>%
  dplyr::mutate(age=as.numeric(age)) %>%
  dplyr::mutate(bagNo=as.numeric(bagNo)) 


df <- estimatesAmNat %>%
  dplyr::group_by(site,round) %>%
  dplyr::summarise(var = var(PVtot,na.rm=TRUE),
                   mu = mean(PVtot,na.rm=TRUE))


ggplot(data=df) +
  geom_abline(slope=1,intercept=0,color='red') +
  geom_point(aes(x=mu,y=var),alpha=0.25) +
  theme_bw()

viabilityComparison2<-df %>%
  dplyr::select(site,round,mu,var) %>%
  dplyr::left_join(viabilityComparison,by=c("site","round"))

ggplot(data=viabilityComparison2) +
  geom_abline(slope=1,intercept=0,color='red') +
  geom_point(aes(x=var,y=p.var),alpha=0.25) +
  theme_bw()


# variance in estimate of PVtot is related to residual in p.med-PV1
ggplot(data=viabilityComparison2) +
  geom_abline(slope=1,intercept=0,color='red') +
  geom_point(aes(x=PVratio,y=p.med,size=1/sqrt(var)),alpha=0.25) +
  theme_bw()

ggplot(data=viabilityComparison2) +
  geom_abline(slope=1,intercept=0,color='red') +
  geom_point(aes(y=p.med-PVratio,x=var),alpha=0.25) +
  theme_bw()

ggplot(data=viabilityComparison2) +
  geom_abline(slope=1,intercept=0,color='red') +
  geom_point(aes(y=p.med,x=var,size=PVwts),alpha=0.25) +
  theme_bw() +facet_wrap(~site)


ggplot(data=viabilityComparison2) +
  geom_abline(slope=1,intercept=0,color='red') +
  geom_point(aes(y=PVratio,x=var),alpha=.5) +
  theme_bw()

ggplot(data=viabilityComparison2) +
  geom_boxplot(aes(x=site,y=var),alpha=0.25) +
  theme_bw() 

ggplot(data=viabilityComparison2) +
  geom_boxplot(aes(x=site,y=p.med-PVratio),alpha=0.25) +
  theme_bw() 

cor(viabilityComparison2$var,viabilityComparison2$p.med-viabilityComparison2$PVratio,use="complete")

df<-mat.list %>% dplyr::filter(site=="SM",round==1) 
plot(df$PVwts,type='b',ylim=c(0,300))
lines(df$PVmean,type='b',col='red')

df2<-viabilityComparison2 %>% dplyr::filter(site=="SM",round==1) 
plot(df$PVmean/df$PVwts,type='b',ylim=c(0,1))
lines(df$PVtot,type='b',col='red')
abline(h=df2$PVratio,lty='dotted')
abline(h=df2$p.med,lty='dotdash')



# visualizing shrinkage
shrinkage<-viabilityComparison2 %>%
  dplyr::select(site,round,p.med,PVratio) %>%
  dplyr::rename(est1_freq=PVratio,est2_bayes=p.med) %>%
  tidyr::pivot_longer(cols=3:4) %>%
  dplyr::filter(!is.na(site)) 

ggplot(data=shrinkage) +
  geom_hline(data=viabilityEstimatesSite,aes(yintercept=est_site),col='red') +
  geom_point(aes(x=name,y=value)) +
  geom_line(aes(x=name,y=value,group=round)) +
  facet_wrap(~site) +
  theme_bw()


seedBagMaster <- readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/seedBagMaster.RDS")

df2<-estimatesAmNat %>%
  dplyr::left_join(viabilityComparison %>%
                     dplyr::select(site,round,PVratio,p.med),by=c("site","round")) %>%
  dplyr::left_join(seedBagMaster,by=c("site","round","age","bagNo")) %>%
  dplyr::filter(!is.na(site)) 


ggplot(df2) +
  geom_hline(aes(yintercept=p.med,color=as.factor(round)),linetype='solid') +
  geom_hline(aes(yintercept=PVratio,color=as.factor(round)),linetype='dotted') +
  geom_boxplot(aes(x=transect,y=PVtot,color=as.factor(round))) +
  facet_wrap(~site) + theme_bw()

ggplot(df2 %>% dplyr::filter(site=="BR"&round==1)) +
  geom_hline(aes(yintercept=p.med,color=as.factor(round)),linetype='solid') +
  geom_hline(aes(yintercept=PVratio,color=as.factor(round)),linetype='dotted') +
  geom_point(aes(x=transect,y=PVtot,color=as.factor(round),size=gs)) +
  facet_wrap(~site) + theme_bw()


out=df2 %>% 
  dplyr::left_join(viabilityEstimatesSite,by=c('site'))


ggplot(data=out %>% dplyr::filter(site=="DEM")) +
  geom_hline(aes(yintercept=est_site ),col='orange') +
  geom_hline(aes(yintercept=p.med,color=as.factor(round)),linetype='solid') +
  geom_hline(aes(yintercept=PVratio,color=as.factor(round)),linetype='dotted') +
  geom_point(aes(x=transect,y=PVtot,color=as.factor(round),size=gs)) +
  facet_wrap(round~site) + theme_bw()


viabFinal <- readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/viabilityFinal.RDS")

viabFinal %>% dplyr::filter(age==1) %>%
  dplyr::group_by(round,block) %>%
  # dplyr::summarise(mu.g = mean(germCount/germStart)) %>%
  ggplot() +
  geom_boxplot(aes(x=block,y=germCount/germStart,group=block)) +
  facet_wrap(~round)

viabFinal %>% dplyr::filter(age==1) %>%
  dplyr::group_by(round,block) %>%
  #dplyr::summarise(mu.v = mean(viabStain/viabStart,na.rm=TRUE)) %>%
  ggplot() +
  geom_boxplot(aes(x=block,y=viabStain/viabStart,group=block)) +
  facet_wrap(~round)

df<-viabFinal %>% dplyr::filter(age==1) %>% dplyr::mutate(viabStain/viabStart)
ref<-viabFinal %>%
  dplyr::select(site,bagNo,round,age,block) %>%
  dplyr::rename(bag.y=bagNo)

df<-viabilityComparison %>%
  dplyr::left_join(ref,by=c("site","bag.y","round","age")) 


# difference 1: I estimate a site*year estimate for total viability
# rather than per bag

# calculate per bag viability and germination probability
# combine to calculate site*year total viability as derived quantity
# could do this at per bag level too I guess and compare to PVtot

# but the site*year combo is what's used to combine with germination trial
# and multi-year dataset

# influence of outliers on estimates will depend on product of mult*JanIntact


# block issues

# fewer seeds in later blocks
out=viabFinal %>% dplyr::group_by(round,age,block) %>% dplyr::summarise(mu = mean(viabStain/viabStart,na.rm=TRUE),n=sum(viabStart))  %>% dplyr::filter(age==1)
ggplot(out) + geom_point(aes(x=block,y=n)) + facet_wrap(~round)
ggplot(out) + geom_point(aes(x=block,y=mu)) + facet_wrap(~round)
ggplot(out) + geom_point(aes(x=n,y=mu)) + facet_wrap(~round)

out=viabFinal %>% dplyr::group_by(round,age,block) %>% dplyr::summarise(mu = mean(germCount/germStart,na.rm=TRUE),n=sum(germStart))  %>% dplyr::filter(age==1)
ggplot(out) + geom_point(aes(x=block,y=n)) + facet_wrap(~round)
ggplot(out) + geom_point(aes(x=block,y=mu)) + facet_wrap(~round)
ggplot(out) + geom_point(aes(x=n,y=mu)) + facet_wrap(~round)




