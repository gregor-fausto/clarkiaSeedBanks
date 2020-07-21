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


# # # -------------------------------------------------------------------
# # # -------------------------------------------------------------------
# # # For testing
# # # -------------------------------------------------------------------
# # # -------------------------------------------------------------------
# keep<-data$site2==12
# data$site2<-data$site2[keep]-11
# data$year2<-data$year2[keep]
# data$y_und<-data$y_und[keep]
# data$n2 = length(data$y_und)
# data$n_site2 = 1

################################################################################
# Posteriors
#################################################################################
# 
#  MCMCsummary(mcmcSamples,params="nu_tfe")
# MCMCsummary(mcmcSamples,params="nu_und")
# MCMCsummary(mcmcSamples,params="nu_dam")
# MCMCsummary(mcmcSamples,params="nu_seeds")
# MCMCsummary(mcmcSamples,params="nu_dam_seeds")

# 
# par(mfrow=c(1,1))
# hist(rgamma(1000,1,1),breaks=100)
# hist(MCMCchains(mcmcSamples,params="nu_und")[,2],breaks=100,add=TRUE,col='red')
# hist(MCMCchains(mcmcSamples,params="tau0_und")[,4],breaks=100)
# 
# par(mfrow=c(2,3))
# for(i in 1:6){
#    hist(MCMCchains(mcmcSamples,params="mu_und")[,i],breaks=100)
# }
# 
# for(i in 1:6){
#    hist(MCMCchains(mcmcSamples,params="tau_und")[,i],breaks=100)
# }
# 
# par(mfrow=c(1,1))
# hist(MCMCchains(mcmcSamples,params="mu_p"),breaks=100)
# 
# par(mfrow=c(2,3))
# for(i in c(1,18)){
#    hist(MCMCchains(mcmcSamples,params="mu_py")[,i],breaks=100)
# }
# 
# par(mfrow=c(1,1))
# hist(MCMCchains(mcmcSamples,params="mean_p")[,1],breaks=100)
# abline(v=mean(data$y_und),col="red",lwd=2)
# 
# 
# par(mfrow=c(2,3))
# for(i in 1:6){
#    hist(MCMCchains(mcmcSamples,params="mean_py")[,i],breaks=100)
#    tmp <- data$y_und[data$year2==i]
#    abline(v=mean(tmp),col="red",lwd=2)
# }
# 
# 
# hist(data$y_und,breaks=100)
# par(mfrow=c(2,3))
# for(i in 1:6){
#    tmp <- data$y_und[data$year2==i]
#    hist(tmp,breaks=50)
# }



################################################################################
# Spatial data
#################################################################################
position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,northing,elevation) %>%
  dplyr::mutate(easting=easting/1000,northing=northing/1000)


 ################################################################################
 # Posteriors for undamaged fruits (2013-)
 #################################################################################
 
countUndamagedDamagedFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countUndamagedDamagedFruitsPerPlantAllPlots.RDS")

siteNames = unique(countUndamagedDamagedFruitsPerPlantAllPlots$site)
 
# recover site indices and year indices
 siteIndex <- data.frame(siteIndex=unique(countUndamagedDamagedFruitsPerPlantAllPlots$site),site=1:20)
 yearIndex <- data.frame(year=1:6,yearIndex=unique(countUndamagedDamagedFruitsPerPlantAllPlots$year))
 
 siteIndex <- data.frame(siteIndex=unique(data$site2),siteBags=1:20)
 yearIndex <- data.frame(yearBags=1:3,yearIndex=2006:2008)
 
 # extract population level distributions for undamaged plants
 ufSummaryPop <- mcmcSamples %>%
   tidybayes::recover_types(data) %>%
   tidybayes::spread_draws(mu_p_und[siteNames]) %>%
    dplyr::group_by(siteNames) %>%
   dplyr::summarise(med = median(mu_p_und,0), 
                    ci.lo = HPDI(mu_p_und,0.95)[1], 
                    ci.hi = HPDI(mu_p_und,0.95)[2],
                    ci.lo2 = HPDI(mu_p_und,0.5)[1], 
                    ci.hi2 = HPDI(mu_p_und,0.5)[2]
   )
 
 ufSummaryPopDF <- cbind(ufSummaryPop,site=siteNames) %>%
   dplyr::left_join(position,by="site")

 ufSummaryPopYear<-mcmcSamples %>%
   tidybayes::recover_types(data) %>%
   tidybayes::spread_draws(mu_py_und[site,year]) %>%
   dplyr::group_by(site,year) %>%
   dplyr::summarise(med = median(mu_py_und), 
                    ci.lo = HPDI(mu_py_und,0.95)[1], 
                    ci.hi = HPDI(mu_py_und,0.95)[2],
                    ci.lo2 = HPDI(mu_py_und,0.5)[1], 
                    ci.hi2 = HPDI(mu_py_und,0.5)[2]
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
 
   # this needs to be fixed
 ufSummaryPopMed = ufSummaryPopDF %>%
    dplyr::rename(lambda.med=med) %>%
    dplyr::select(site,lambda.med)
 
interannualUF <- mcmcSamples %>%
   tidybayes::recover_types(data) %>%
   tidybayes::spread_draws(mu_py_und[site,year]) %>%
   dplyr::group_by(site,year) %>%
    dplyr::summarise(med = median(mu_py_und),
                     ci.lo = HPDI(mu_py_und,0.89)[1],
                     ci.hi = HPDI(mu_py_und,0.89)[2]
    ) %>%
   dplyr::ungroup() %>%
   dplyr::left_join(siteIndex,by="site") %>%
   dplyr::left_join(yearIndex,by="year") %>%
   dplyr::select(-c(site,year)) %>%
   dplyr::rename(site = siteIndex) %>%
   dplyr::rename(year = yearIndex) %>%
    
   dplyr::left_join(ufSummaryPopMed, by="site") %>%
   dplyr::arrange(lambda.med) %>% 
   dplyr::mutate(site=factor(site,levels=unique(site))) %>%
   dplyr::group_by(site) %>%
   dplyr::arrange(desc(med),.by_group = TRUE) %>%
   dplyr::mutate(year=factor(year)) %>%
   mutate(id = row_number()) %>% 
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
 