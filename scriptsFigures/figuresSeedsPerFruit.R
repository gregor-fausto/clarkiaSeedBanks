#################################################################################
################################################################################
################################################################################
# Code for figures to compare the following modeling approaches for the seedling survivorship data
#
# Scripts by Gregor Siegmund
# fausto.siegmund@gmail.com
# last updated 04-22-2020
#################################################################################
#################################################################################
#################################################################################
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)

dirFigures = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/"

library(MCMCvis)
library(tidybayes)
library(tidyverse)
library(magrittr)
library(bayesplot)
library(HDInterval)
library(coda)
library(rethinking)

directory = "/Users/Gregor/Dropbox/dataLibrary/posteriors/"
posteriorFiles <- paste0(directory,list.files(directory))

mcmcSamples <- readRDS(posteriorFiles[[7]])

################################################################################
# Read in raw data
#################################################################################

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/"
dataFiles <- paste0(directory,list.files(directory))

data <- readRDS(dataFiles[[5]])

################################################################################
# Spatial data
#################################################################################
position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,northing,elevation) %>%
  dplyr::mutate(easting=easting/1000,northing=northing/1000)


################################################################################
# Posteriors for undamaged fruits (2013-)
#################################################################################

countSeedPerFruit <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countSeedPerFruit.RDS")
countSeedPerFruit <- countSeedPerFruit %>% dplyr::filter(demography==1)
siteNames = unique(countSeedPerFruit$site)

# recover site indices and year indices
siteIndex <- data.frame(siteIndex=unique(countSeedPerFruit$site),site=unique(data$site3))
yearIndex <- data.frame(yearIndex=unique(countSeedPerFruit$year),
                        year=unique(data$year3)) 

# extract population level distributions for undamaged plants
seeds_p <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(mu_p_seeds[site]) %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(med = median(mu_p_seeds,0), 
                   ci.lo = HPDI(mu_p_seeds,0.95)[1], 
                   ci.hi = HPDI(mu_p_seeds,0.95)[2],
                   ci.lo2 = HPDI(mu_p_seeds,0.5)[1], 
                   ci.hi2 = HPDI(mu_p_seeds,0.5)[2]
  )

seeds_p <- cbind(seeds_p,siteIndex=siteNames)

seeds_py<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(mu_py_seeds[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(med = median(mu_py_seeds), 
                   ci.lo = HPDI(mu_py_seeds,0.95)[1], 
                   ci.hi = HPDI(mu_py_seeds,0.95)[2],
                   ci.lo2 = HPDI(mu_py_seeds,0.5)[1], 
                   ci.hi2 = HPDI(mu_py_seeds,0.5)[2]
  )

seeds_py$med[seeds_py$site=="DLW"&seeds_py$year==2014]

seeds_py <- seeds_py %>%
  dplyr::mutate(year = as.integer(year)) %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(year)) %>%
  dplyr::select(-c(site)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)

# this needs to be fixed
seeds_p.med = seeds_p %>%
  dplyr::rename(lambda.med=med) %>%
  dplyr::select(site,lambda.med) %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::select(-site) %>%
  dplyr::rename(site = siteIndex)

# need to figure out ci.hi
interannualDF<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(mu_py_seeds[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(med = median(mu_py_seeds),
                   ci.lo = HPDI(mu_py_seeds,0.89)[1],
                   ci.hi = HPDI(mu_py_seeds,0.89)[2]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  dplyr::select(-c(site,year)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::left_join(seeds_p.med, by="site") %>%
  dplyr::arrange(lambda.med) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(desc(med),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  mutate(id = row_number())

interannualSeeds <- interannualDF %>% 
  ggplot(aes(x = id , y = med)) + 
  geom_hline(aes(yintercept=lambda.med),linetype='dotted') +
  geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +  
  
  geom_point(aes(color=year),size=1) +
  coord_flip() +
  facet_grid(site ~ ., scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ylab("Seeds per fruit (undamaged fruits)") +
  ylim(c(0,200))

spatialSeeds <-seeds_p %>%
  dplyr::select(-site) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::left_join(position,by="site") %>%
  ggplot(aes(x = easting , y = med)) + 
  geom_point() +
  geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  ylab("Seeds per fruit (undamaged fruits)") +
  xlab("Easting (km)") +
  labs(title="Seeds per fruit (2006-2018)") 

# ribbon plot
GeomRibbon$handle_na <- function(data, params) {  data }

# time plot
ribbonSeeds <-   interannualDF %>%
  dplyr::left_join(position,by="site") %>%
  ggplot(aes(x=year,y=med,group=site,ymin=ci.lo,ymax=ci.hi)) +
  geom_ribbon(alpha=0.5,colour="grey50") +
  geom_line() +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  xlab("Year") +
  ylab("Seeds per fruit (undamaged fruits)") +
  facet_wrap(~site,scales='free_y') +
  xlab("Easting (km)") +
  labs(title="Seeds per fruit (2006-2018)") 

## Save Figures
ggsave(filename=paste0(dirFigures,"interannualSeeds.pdf"),
       plot=interannualSeeds,width=6,height=16)

ggsave(filename=paste0(dirFigures,"spatialSeeds.pdf"),
       plot=spatialSeeds,width=4,height=4)

ggsave(filename=paste0(dirFigures,"ribbonSeeds.pdf"),
       plot=ribbonSeeds,width=12,height=6)
