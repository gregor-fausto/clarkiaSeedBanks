rm(list=ls(all=TRUE)) # clear R environment

# load packages

library(MCMCvis)
library(tidybayes)
library(tidyverse)
library(magrittr)
library(bayesplot)
library(rethinking)

# Summary

# Viability is generally high: only SM has nu<0.75
# and especially for the first year, these differences are small
# because it translates to a difference of ~5% in January
# noteworthy that the sites with lowest viability are SM, EC, BR
# all of which are in the southeast corner of the range

# these are all for age 1 viabilities
# Load MCMC samples -------------------------------------------------------

mcmcSamples<-readRDS(file=paste0("/Users/Gregor/Dropbox/dataLibrary/posteriors/modelChecking/viabilitySamples.RDS"))

#directory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/"
#simFiles <- paste0(directory,list.files(directory))

data <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/viabilityData.RDS")
viabilityRawData1<-readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/viabilityRawData1.RDS")
siteNames = unique(viabilityRawData1$siteViab)

siteIndex <- data.frame(siteIndex=siteNames,site=1:20)
yearIndex <- data.frame(yearViab=1:3,year=2006:2008)

nu_py<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(nu0_1[siteViab,yearViab]) %>%
  dplyr::group_by(siteViab,yearViab) %>%
  dplyr::summarise(med = median(nu0_1), 
                   ci.lo = HPDI(nu0_1,.89)[1], 
                   ci.hi = HPDI(nu0_1,.89)[2],
                   ci.lo2 = HPDI(nu0_1,.5)[1],
                   ci.hi2 = HPDI(nu0_1,.5)[2]
  )

nu_py.df<-nu_py %>%
  dplyr::mutate(yearViab = as.integer(yearViab)) %>%
  dplyr::left_join(yearIndex,by="yearViab") %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(yearViab)) %>%
  dplyr::rename(site = siteViab) %>%
  dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)

nu_p<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(nu_1[siteViab]) %>%
  dplyr::group_by(siteViab) %>%
  dplyr::summarise(med = median(nu_1), 
                   ci.lo = HPDI(nu_1,.89)[1], 
                   ci.hi = HPDI(nu_1,.89)[2],
                   ci.lo2 = HPDI(nu_1,.5)[1],
                   ci.hi2 = HPDI(nu_1,.5)[2]
  )

nu_p.df<-nu_p %>%
  dplyr::rename(site=siteViab) %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(site)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::select(site,med,ci.lo,ci.hi,ci.lo2,ci.hi2) %>%
  dplyr::rename(med_p=med,ci.lo_p=ci.lo,ci.hi_p=ci.hi,
                ci.lo2_p=ci.lo2,ci.hi2_p=ci.hi2)


# by site
nu_py.df<-nu_py.df %>%
  dplyr::ungroup() %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::select(-c(site)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::left_join(nu_p.df, by="site") %>%
  dplyr::arrange(med_p) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(desc(-med),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  mutate(id = row_number()) 
  
ggplot(data=nu_py.df,aes(x=site,y=med)) + 

  geom_hline(aes(yintercept=med_p),linetype='dotted') +
  geom_hline(aes(yintercept=ci.lo_p),linetype='solid',alpha=0.5) +
  geom_hline(aes(yintercept=ci.hi_p),linetype='solid',alpha=0.5) +
  geom_pointrange(
    aes(ymin = ci.lo, ymax = ci.hi, group=id, color = year),
    position = position_dodge(.5)
  ) +
  facet_grid(. ~ site, scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  ylab("Probability of composite viability") +
  ylim(c(0,1))



nu_py.df<-nu_py.df %>%
  dplyr::ungroup() %>%
  dplyr::arrange(year,med) %>%
  dplyr::group_by(year) %>%
  mutate(id = row_number()) #%>% 
 
ggplot(data=nu_py.df,aes(x=year,y=med)) + 
  geom_pointrange(
    aes(ymin = ci.lo, ymax = ci.hi, group = id),
    position = position_dodge(1)
  ) +
  facet_grid(. ~ year, scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  ylab("Probability of composite viability") +
  ylim(c(0,1))
