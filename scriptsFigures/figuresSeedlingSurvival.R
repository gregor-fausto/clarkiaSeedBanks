#################################################################################
################################################################################
################################################################################
# Code for figures to compare the following modeling approaches for the seedling survivorship data
#
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

directory = "/Users/Gregor/Dropbox/dataLibrary/posteriors/"
simFiles <- paste0(directory,list.files(directory))

mcmcSamples <- readRDS(simFiles[[12]])
dirFigures = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/"

################################################################################
# Data directory
#################################################################################

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedlingSurvival/"
dataFiles <- paste0(directory,list.files(directory))

data <- readRDS(dataFiles[[2]])
censusSeedlingsFruitingPlants <- readRDS(dataFiles[[1]])

################################################################################
# Spatial data
#################################################################################
position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,northing,elevation) %>%
  dplyr::mutate(easting=easting/1000,northing=northing/1000)

################################################################################
# ...
#################################################################################

siteNames <- unique(censusSeedlingsFruitingPlants$site)

parsToMonitor_1 = c("theta_1","mu0_1","sigma0_1","mu_1","sigma_1")
parsToMonitor_deriv = c("p0_1","p_1")

sigma_p <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p0_1[site]) %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(med = median(p0_1), 
                   ci.lo = HPDI(p0_1,.89)[1], 
                   ci.hi = HPDI(p0_1,.89)[2],
                   ci.lo2 = HPDI(p0_1,.5)[1],
                   ci.hi2 = HPDI(p0_1,.5)[2]
  )

sigma_p<-cbind(sigma_p,siteNames)

sigma_p.med = sigma_p %>% dplyr::select(siteNames,med) %>% 
  dplyr::rename(site=siteNames) 

siteIndex <- data.frame(siteIndex=unique(censusSeedlingsFruitingPlants$site),site=1:20)
yearIndex <- data.frame(yearIndex=as.numeric(unique(censusSeedlingsFruitingPlants$year)),
                        year=unique(unique(data$year)))

sigma_py<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p_1[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(med = median(p_1), 
                   ci.lo = HPDI(p_1,.89)[1], 
                   ci.hi = HPDI(p_1,.89)[2],
                   ci.lo2 = HPDI(p_1,.5)[1],
                   ci.hi2 = HPDI(p_1,.5)[2]
  )

sigma_py<-sigma_py %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(site,year)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)


interannualSigmaDF<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p_1[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(lambda.med = median(p_1), 
                   ci.lo = HPDI(p_1,.89)[1], 
                   ci.hi = HPDI(p_1,.89)[2]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  dplyr::select(-c(site,year)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::left_join(sigma_p.med, by="site") %>%
  dplyr::arrange(med) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(desc(lambda.med),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  mutate(id = row_number())
  
interannualSigma <- interannualSigmaDF %>%
  ggplot(aes(x = id , y = lambda.med)) + 
  geom_hline(aes(yintercept=med),linetype='dotted') +
  geom_point(aes(color=year)) +
  geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +
  coord_flip() +
  facet_grid(site ~ ., scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ylab("Probability of seedling survival to fruiting") +
  ylim(c(0,1)) +
  labs(color="Year") 

spatialSigma <-sigma_p %>%
  dplyr::select(-site) %>%
  dplyr::rename(site=siteNames) %>%
  dplyr::left_join(position,by="site") %>%
  ggplot(aes(x = easting , y = med)) + 
  geom_point() +
  geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  ylab("Probability of seedling survival to fruiting") +
  xlab("Easting (km)") +
  ylim(c(0,1))

ggsave(filename=paste0(dirFigures,"interannualSigma.pdf"),
       plot=interannualSigma,width=6,height=12)

ggsave(filename=paste0(dirFigures,"spatialSigma.pdf"),
       plot=spatialSigma,width=6,height=12)


allSurvival = interannualSigmaDF %>%
  dplyr::ungroup() %>%
  dplyr::left_join(siteIndex %>% 
                     dplyr::rename(siteNumber=site,site=siteIndex),by="site") %>%
  dplyr::arrange((siteNumber)) %>%
  dplyr::mutate(year=as.integer(as.character(year))) %>%
  ggplot(aes(x=year,y=lambda.med,group=site,ymin=ci.lo,ymax=ci.hi)) +
  geom_ribbon(alpha=0.1,fill="cornflowerblue") +
  geom_line(color="black") +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  ylab("Probability of seedling survival to fruiting") +
  scale_x_continuous(breaks = c(2006, 2008, 2010, 2012, 2014, 2016, 2018),
                     labels = c(2006, "", 2010, "", 2014, "", 2018)) +
  scale_size_continuous(range = c(.5,2.5)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none") +
  xlab("Year")

popSurvival <- interannualSigmaDF %>%
  dplyr::ungroup() %>%
  dplyr::left_join(siteIndex %>% 
                     dplyr::rename(siteNumber=site,site=siteIndex),by="site") %>%
  dplyr::arrange((siteNumber)) %>%
  dplyr::mutate(year=as.integer(as.character(year))) %>%
  ggplot(aes(x=year,y=lambda.med,group=site,ymin=ci.lo,ymax=ci.hi)) +
  geom_rect(aes(xmin = 2013, xmax = 2017, ymin = -Inf, ymax = 0),fill='gold') +
  geom_ribbon(alpha=0.25,fill="cornflowerblue") +
  geom_line(color="black") +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  ylab("Probability of seedling survival to fruiting") +
  scale_x_continuous(breaks = c(2006, 2008, 2010, 2012, 2014, 2016, 2018),
                     labels = c(2006, "", 2010, "", 2014, "", 2018)) +
  scale_size_continuous(range = c(.5,2.5)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none") +
  xlab("Year") + facet_wrap(~site)

ggsave(filename=paste0("~/Dropbox/chapter-3/figures/all-survival.pdf"),
       plot=allSurvival,width=12,height=6)

ggsave(filename=paste0("~/Dropbox/chapter-3/figures/pop-survival.pdf"),
       plot=popSurvival,width=12,height=12)


sampleSize = censusSeedlingsFruitingPlants %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(year = as.double(year)) 


trialNumber = censusSeedlingsFruitingPlants %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(trials = sum(seedlingNumber),
                   trialsVar = var(seedlingNumber)) %>%
  dplyr::mutate(year = as.double(year)) 


df <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p_1[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(lambda.med = median(p_1), 
                   ci.lo = HPDI(p_1,.89)[1], 
                   ci.hi = HPDI(p_1,.89)[2]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  dplyr::select(-c(site,year)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::left_join(sampleSize,by=c("site","year")) %>%
  dplyr::left_join(trialNumber,by=c("site","year"))

df<- df %>%
  dplyr::left_join(sigma_p.med, by="site") %>%
  dplyr::arrange(med) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(desc(lambda.med),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  mutate(id = row_number())

df %>%
  ggplot(aes(x = id , y = lambda.med)) + 
  geom_hline(aes(yintercept=med),linetype='dotted') +
  geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25,alpha=0.25) +
  geom_point(aes(color=year,size=trials),alpha=0.5) +
  coord_flip() +
  facet_grid(site ~ ., scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ylab("Probability of seedling survival to fruiting") +
  ylim(c(0,1)) +
  labs(color="Year") 

df %>%
  dplyr::filter(site=="BG")

censusSeedlingsFruitingPlants %>%
  dplyr::filter(site=="BG",year==2015) %>% View
