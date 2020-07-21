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

library(MCMCvis)
library(tidybayes)
library(tidyverse)
library(magrittr)
library(bayesplot)

directory = "/Users/Gregor/Dropbox/dataLibrary/posteriors/"
posteriorFiles <- paste0(directory,list.files(directory))

mcmcSamples <- readRDS(posteriorFiles[[4]])


#saveRDS(data,file=paste0(fileDirectory,"data.rds"))
#saveRDS(countFruitsPerPlantAllPlots,file=paste0(fileDirectory,"countFruitsPerPlantAllPlots.rds"))
#saveRDS(countUndamagedDamagedFruitsPerPlantAllPlots,file=paste0(fileDirectory,"countUndamagedDamagedFruitsPerPlantAllPlots.rds"))
#saveRDS(countSeedPerUndamagedFruit,file=paste0(fileDirectory,"countSeedPerUndamagedFruit.rds"))
#saveRDS(countSeedPerDamagedFruit,file=paste0(fileDirectory,"countSeedPerDamagedFruit.rds"))


################################################################################
# Read in raw data
#################################################################################

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/"
simFiles <- paste0(directory,list.files(directory))

data <- readRDS(simFiles[[5]])

countFruitsPerPlantAllPlots <- readRDS(simFiles[[1]])
countUndamagedDamagedFruitsPerPlantAllPlots <- readRDS(simFiles[[4]])
countSeedPerUndamagedFruit <- readRDS(simFiles[[3]])

siteNames = unique(countFruitsPerPlantAllPlots$site)

################################################################################
# Spatial data
#################################################################################
position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,northing,elevation) %>%
  dplyr::mutate(easting=easting/1000,northing=northing/1000)

################################################################################
# Posteriors for total fruit equivalents (2006-2012)
#################################################################################

# recover site indices and year indices
siteIndex <- data.frame(siteIndex=unique(countFruitsPerPlantAllPlots$site),site=1:20)
yearIndex <- data.frame(year=1:data$n_year,yearIndex=unique(countFruitsPerPlantAllPlots$year))

# extract population level distributions
tfeSummaryPop <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p0_tfe[siteNames]) %>%
  dplyr::summarise(med = median(p0_tfe), 
                   ci.lo = quantile(p0_tfe,probs=0.025), 
                   ci.hi = quantile(p0_tfe,probs=0.975),
                   ci.lo2 = quantile(p0_tfe,probs=0.25), 
                   ci.hi2 = quantile(p0_tfe,probs=0.75)
  )

tfeSummaryPopDF <- cbind(tfeSummaryPop,site=siteNames) %>%
  dplyr::left_join(position,by="site")

ggplot(data=tfeSummaryPopDF) +
  geom_point(aes(x=easting,y=med),size=1) +
  geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.5) +
  geom_linerange(aes(x=easting,ymin=ci.lo2,ymax=ci.hi2),size=1) +
  theme_bw() +
  xlab("Easting (km)") +
  ylab("Total fruit equivalents per plant") +
  labs(title="Total fruit equivalents per plant (2006-2012)")

tfeSummaryPopYear<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(tfe[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(med = median(tfe), 
                   ci.lo = quantile(tfe,probs=0.025), 
                   ci.hi = quantile(tfe,probs=0.975),
                   ci.lo2 = quantile(tfe,probs=0.25), 
                   ci.hi2 = quantile(tfe,probs=0.75)
  )

tfeSummaryPopYearDF <- tfeSummaryPopYear %>%
  dplyr::mutate(year = as.integer(year)) %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(year)) %>%
  dplyr::select(-c(site)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)

tfeSummaryPopMed = tfeSummaryPopDF %>% dplyr::select(site,med)

 mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(tfe[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(lambda.med = median(tfe), 
                   ci.lo = quantile(tfe,probs=0.025), 
                   ci.hi = quantile(tfe,probs=0.975),
                   ci.lo2 = quantile(tfe,probs=0.25), 
                   ci.hi2 = quantile(tfe,probs=0.75)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  dplyr::select(-c(site,year)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::left_join(tfeSummaryPopMed, by="site") %>%
  dplyr::arrange(med) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(desc(lambda.med),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  mutate(id = row_number()) %>% 
  ggplot(aes(x = id , y = lambda.med)) + 
  geom_hline(aes(yintercept=med),linetype='dotted') +
   geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +
   geom_linerange(aes(x=id,ymin=ci.lo2,ymax=ci.hi2),size=1) +  
  geom_point(aes(color=year),size=1, shape=1) +

  coord_flip() +
  facet_grid(site ~ ., scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ylab("Fruits per plant (total fruit equivalents)") 

 ################################################################################
 # Posteriors for undamaged fruits (2013-)
 #################################################################################
 
 # recover site indices and year indices
 siteIndex <- data.frame(siteIndex=unique(countUndamagedDamagedFruitsPerPlantAllPlots$site2),site=1:20)
 yearIndex <- data.frame(year=1:data$n_year2,yearIndex=unique(countUndamagedDamagedFruitsPerPlantAllPlots$year2))
 
 # extract population level distributions for undamaged plants
 ufSummaryPop <- mcmcSamples %>%
   tidybayes::recover_types(data) %>%
   tidybayes::spread_draws(p0_und[siteNames]) %>%
   dplyr::summarise(med = median(p0_und), 
                    ci.lo = quantile(p0_und,probs=0.025), 
                    ci.hi = quantile(p0_und,probs=0.975),
                    ci.lo2 = quantile(p0_und,probs=0.25), 
                    ci.hi2 = quantile(p0_und,probs=0.75)
   )
 
 ufSummaryPopDF <- cbind(ufSummaryPop,site=siteNames) %>%
   dplyr::left_join(position,by="site")
 
 ggplot(data=ufSummaryPopDF) +
   geom_point(aes(x=easting,y=med),size=1) +
   geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.5) +
   geom_linerange(aes(x=easting,ymin=ci.lo2,ymax=ci.hi2),size=1) +
   theme_bw() +
   xlab("Easting (km)") +
   ylab("Undamaged fruits per plant") +
   labs(title="Undamaged fruits per plant (2013-2018)")
 
 ufSummaryPopYear<-mcmcSamples %>%
   tidybayes::recover_types(data) %>%
   tidybayes::spread_draws(p_und[site,year]) %>%
   dplyr::group_by(site,year) %>%
   dplyr::summarise(med = median(p_und), 
                    ci.lo = quantile(p_und,probs=0.025), 
                    ci.hi = quantile(p_und,probs=0.975),
                    ci.lo2 = quantile(p_und,probs=0.25), 
                    ci.hi2 = quantile(p_und,probs=0.75)
   )
 
 ufSummaryPopYearDF <- ufSummaryPopYear %>%
   dplyr::mutate(year = as.integer(year)) %>%
   dplyr::left_join(siteIndex,by="site") %>%
   dplyr::left_join(yearIndex,by="year") %>%
   dplyr::ungroup() %>%
   dplyr::select(-c(year)) %>%
   dplyr::select(-c(site)) %>%
   dplyr::rename(site = siteIndex) %>%
   dplyr::rename(year = yearIndex) %>%
   dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)
 
 ufSummaryPopMed = ufSummaryPopDF %>% dplyr::select(site,med)
 
 mcmcSamples %>%
   tidybayes::recover_types(data) %>%
   tidybayes::spread_draws(p_und[site,year]) %>%
   dplyr::group_by(site,year) %>%
   dplyr::summarise(lambda.med = median(p_und), 
                    ci.lo = quantile(p_und,probs=0.025), 
                    ci.hi = quantile(p_und,probs=0.975),
                    ci.lo2 = quantile(p_und,probs=0.25), 
                    ci.hi2 = quantile(p_und,probs=0.75)
   ) %>%
   dplyr::ungroup() %>%
   dplyr::left_join(siteIndex,by="site") %>%
   dplyr::left_join(yearIndex,by="year") %>%
   dplyr::select(-c(site,year)) %>%
   dplyr::rename(site = siteIndex) %>%
   dplyr::rename(year = yearIndex) %>%
   dplyr::left_join(ufSummaryPopMed, by="site") %>%
   dplyr::arrange(med) %>% 
   dplyr::mutate(site=factor(site,levels=unique(site))) %>%
   dplyr::group_by(site) %>%
   dplyr::arrange(desc(lambda.med),.by_group = TRUE) %>%
   dplyr::mutate(year=factor(year)) %>%
   mutate(id = row_number()) %>% 
   ggplot(aes(x = id , y = lambda.med)) + 
   geom_hline(aes(yintercept=med),linetype='dotted') +
   geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +
   geom_linerange(aes(x=id,ymin=ci.lo2,ymax=ci.hi2),size=1) +  
   geom_point(aes(color=year),size=1, shape=1) +
   
   coord_flip() +
   facet_grid(site ~ ., scales="free_x", space="free_x") +
   theme_bw() +
   theme(panel.spacing=unit(0,"pt"), 
         panel.border=element_rect(colour="grey50", fill=NA)) +
   theme(axis.title.y=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank()) +
   ylab("Fruits per plant (undamaged fruits)")
 
 
 ################################################################################
 # Posteriors for damaged fruits (2013-)
 #################################################################################
 
 # recover site indices and year indices
 siteIndex <- data.frame(siteIndex=unique(countUndamagedDamagedFruitsPerPlantAllPlots$site2),site=1:20)
 yearIndex <- data.frame(year=1:data$n_year2,yearIndex=unique(countUndamagedDamagedFruitsPerPlantAllPlots$year2))
 
 # extract population level distributions for undamaged plants
 dfSummaryPop <- mcmcSamples %>%
   tidybayes::recover_types(data) %>%
   tidybayes::spread_draws(p0_dam[siteNames]) %>%
   dplyr::summarise(med = median(p0_dam), 
                    ci.lo = quantile(p0_dam,probs=0.025), 
                    ci.hi = quantile(p0_dam,probs=0.975),
                    ci.lo2 = quantile(p0_dam,probs=0.25), 
                    ci.hi2 = quantile(p0_dam,probs=0.75)
   )
 
 dfSummaryPopDF <- cbind(dfSummaryPop,site=siteNames) %>%
   dplyr::left_join(position,by="site")
 
 ggplot(data=dfSummaryPopDF) +
   geom_point(aes(x=easting,y=med),size=1) +
   geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.5) +
   geom_linerange(aes(x=easting,ymin=ci.lo2,ymax=ci.hi2),size=1) +
   theme_bw() +
   xlab("Easting (km)") +
   ylab("Damaged fruits per plant") +
   labs(title="Damaged fruits per plant (2013-2018)") 
 
 dfSummaryPopYear<-mcmcSamples %>%
   tidybayes::recover_types(data) %>%
   tidybayes::spread_draws(p_dam[site,year]) %>%
   dplyr::group_by(site,year) %>%
   dplyr::summarise(med = median(p_dam), 
                    ci.lo = quantile(p_dam,probs=0.025), 
                    ci.hi = quantile(p_dam,probs=0.975),
                    ci.lo2 = quantile(p_dam,probs=0.25), 
                    ci.hi2 = quantile(p_dam,probs=0.75)
   )
 
 dfSummaryPopYearDF <- dfSummaryPopYear %>%
   dplyr::mutate(year = as.integer(year)) %>%
   dplyr::left_join(siteIndex,by="site") %>%
   dplyr::left_join(yearIndex,by="year") %>%
   dplyr::ungroup() %>%
   dplyr::select(-c(year)) %>%
   dplyr::select(-c(site)) %>%
   dplyr::rename(site = siteIndex) %>%
   dplyr::rename(year = yearIndex) %>%
   dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)
 
 dfSummaryPopMed = dfSummaryPopDF %>% dplyr::select(site,med)
 
 mcmcSamples %>%
   tidybayes::recover_types(data) %>%
   tidybayes::spread_draws(p_dam[site,year]) %>%
   dplyr::group_by(site,year) %>%
   dplyr::summarise(lambda.med = median(p_dam), 
                    ci.lo = quantile(p_dam,probs=0.025), 
                    ci.hi = quantile(p_dam,probs=0.975),
                    ci.lo2 = quantile(p_dam,probs=0.25), 
                    ci.hi2 = quantile(p_dam,probs=0.75)
   ) %>%
   dplyr::ungroup() %>%
   dplyr::left_join(siteIndex,by="site") %>%
   dplyr::left_join(yearIndex,by="year") %>%
   dplyr::select(-c(site,year)) %>%
   dplyr::rename(site = siteIndex) %>%
   dplyr::rename(year = yearIndex) %>%
   dplyr::left_join(dfSummaryPopMed, by="site") %>%
   dplyr::arrange(med) %>% 
   dplyr::mutate(site=factor(site,levels=unique(site))) %>%
   dplyr::group_by(site) %>%
   dplyr::arrange(desc(lambda.med),.by_group = TRUE) %>%
   dplyr::mutate(year=factor(year)) %>%
   mutate(id = row_number()) %>% 
   ggplot(aes(x = id , y = lambda.med)) + 
   geom_hline(aes(yintercept=med),linetype='dotted') +
   geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +
   geom_linerange(aes(x=id,ymin=ci.lo2,ymax=ci.hi2),size=1) +  
   geom_point(aes(color=year),size=1, shape=1) +
   
   coord_flip() +
   facet_grid(site ~ ., scales="free_x", space="free_x") +
   theme_bw() +
   theme(panel.spacing=unit(0,"pt"), 
         panel.border=element_rect(colour="grey50", fill=NA)) +
   theme(axis.title.y=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank()) +
   ylab("Fruits per plant (undamaged fruits)") +
   ylim(c(0,100))
 
################################################################################
# Posteriors for TFE composite (2013-2018)
#################################################################################
 
 # recover site indices and year indices
 siteIndex <- data.frame(siteIndex=unique(countSeedPerUndamagedFruit$site3),site=1:20)
 yearIndex <- data.frame(year=1:data$n_year3,yearIndex=unique(countSeedPerUndamagedFruit$year3))
 
 # extract population level distributions for undamaged plants
 # usSummaryPop <- mcmcSamples %>%
 #   tidybayes::recover_types(data) %>%
 #   tidybayes::spread_draws(tfe_comp[siteNames]) %>%
 #   dplyr::summarise(med = median(p0_seeds), 
 #                    ci.lo = quantile(p0_seeds,probs=0.025), 
 #                    ci.hi = quantile(p0_seeds,probs=0.975),
 #                    ci.lo2 = quantile(p0_seeds,probs=0.25), 
 #                    ci.hi2 = quantile(p0_seeds,probs=0.75)
 #   )
 # 
 # usSummaryPopDF <- cbind(usSummaryPop,site=siteNames) %>%
 #   dplyr::left_join(position,by="site")
 
 # ggplot(data=usSummaryPopDF) +
 #   geom_point(aes(x=easting,y=med),size=1) +
 #   geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.5) +
 #   geom_linerange(aes(x=easting,ymin=ci.lo2,ymax=ci.hi2),size=1) +
 #   theme_bw() +
 #   xlab("Easting (km)") +
 #   ylab("Seeds per undamaged fruit") +
 #   labs(title="Seeds per undamaged fruit (2006-2018)") 
 # 
 tfecompSummaryPopYear<-mcmcSamples %>%
   tidybayes::recover_types(data) %>%
   tidybayes::spread_draws(tfe_comp[site,year]) %>%
   dplyr::group_by(site,year) %>%
   dplyr::summarise(med = median(tfe_comp),
                    ci.lo = quantile(tfe_comp,probs=0.025),
                    ci.hi = quantile(tfe_comp,probs=0.975),
                    ci.lo2 = quantile(tfe_comp,probs=0.25),
                    ci.hi2 = quantile(tfe_comp,probs=0.75)
   )

 tfecompSummaryPopYearDF <-  tfecompSummaryPopYear %>%
   dplyr::mutate(year = as.integer(year)) %>%
   dplyr::left_join(siteIndex,by="site") %>%
   dplyr::left_join(yearIndex,by="year") %>%
   dplyr::ungroup() %>%
   dplyr::select(-c(year)) %>%
   dplyr::select(-c(site)) %>%
   dplyr::rename(site = siteIndex) %>%
   dplyr::rename(year = yearIndex) %>%
   dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)

 #tfecompSummaryPopMed =  tfecompSummaryPopDF %>% dplyr::select(site,med)
 tfecompSummaryPopMed =  tfecompSummaryPopYearDF %>% 
   dplyr::group_by(site) %>%
   dplyr::summarise(med=mean(med)) %>%
   dplyr::select(site,med)
 
 
 mcmcSamples %>%
   tidybayes::recover_types(data) %>%
   tidybayes::spread_draws(tfe_comp[site,year]) %>%
   dplyr::group_by(site,year) %>%
   dplyr::summarise(lambda.med = median(tfe_comp),
                    ci.lo = quantile(tfe_comp,probs=0.025),
                    ci.hi = quantile(tfe_comp,probs=0.975),
                    ci.lo2 = quantile(tfe_comp,probs=0.25),
                    ci.hi2 = quantile(tfe_comp,probs=0.75)
   ) %>%
   dplyr::ungroup() %>%
   dplyr::left_join(siteIndex,by="site") %>%
   dplyr::left_join(yearIndex,by="year") %>%
   dplyr::select(-c(site,year)) %>%
   dplyr::rename(site = siteIndex) %>%
   dplyr::rename(year = yearIndex) %>%
   dplyr::left_join(tfecompSummaryPopMed, by="site") %>%
   dplyr::arrange(med) %>%
   dplyr::mutate(site=factor(site,levels=unique(site))) %>%
   dplyr::group_by(site) %>%
   dplyr::arrange(desc(lambda.med),.by_group = TRUE) %>%
   dplyr::mutate(year=factor(year)) %>%
   mutate(id = row_number()) %>%
   ggplot(aes(x = id , y = lambda.med)) +
   geom_hline(aes(yintercept=med),linetype='dotted') +
   geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +
   geom_linerange(aes(x=id,ymin=ci.lo2,ymax=ci.hi2),size=1) +
   geom_point(aes(color=year),size=1, shape=1) +

   coord_flip() +
   facet_grid(site ~ ., scales="free_x", space="free_x") +
   theme_bw() +
   theme(panel.spacing=unit(0,"pt"),
         panel.border=element_rect(colour="grey50", fill=NA)) +
   theme(axis.title.y=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank()) +
   ylab("Fruits per plant (undamaged fruits)") +
   ylim(c(0,100))
 