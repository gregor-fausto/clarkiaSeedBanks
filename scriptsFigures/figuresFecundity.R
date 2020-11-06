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
 siteIndex <- data.frame(siteIndex=unique(countUndamagedDamagedFruitsPerPlantAllPlots$site),site=unique(data$site2))
 yearIndex <- data.frame(yearIndex=unique(countUndamagedDamagedFruitsPerPlantAllPlots$year),
                         year=unique(data$year2)) 
 
 # extract population level distributions for undamaged plants
 uf_p <- mcmcSamples %>%
   tidybayes::recover_types(data) %>%
   tidybayes::spread_draws(mu_p_und[site]) %>%
    dplyr::group_by(site) %>%
   dplyr::summarise(med = median(mu_p_und,0), 
                    ci.lo = HPDI(mu_p_und,0.95)[1], 
                    ci.hi = HPDI(mu_p_und,0.95)[2],
                    ci.lo2 = HPDI(mu_p_und,0.5)[1], 
                    ci.hi2 = HPDI(mu_p_und,0.5)[2]
   )
 
 uf_p <- cbind(uf_p,siteIndex=siteNames)

 uf_py<-mcmcSamples %>%
   tidybayes::recover_types(data) %>%
   tidybayes::spread_draws(mu_py_und[site,year]) %>%
   dplyr::group_by(site,year) %>%
   dplyr::summarise(med = median(mu_py_und), 
                    ci.lo = HPDI(mu_py_und,0.95)[1], 
                    ci.hi = HPDI(mu_py_und,0.95)[2],
                    ci.lo2 = HPDI(mu_py_und,0.5)[1], 
                    ci.hi2 = HPDI(mu_py_und,0.5)[2]
   )
 
 uf_py <- uf_py %>%
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
 uf_p.med = uf_p %>%
    dplyr::rename(lambda.med=med) %>%
    dplyr::select(site,lambda.med) %>%
    dplyr::left_join(siteIndex,by="site") %>%
    dplyr::select(-site) %>%
    dplyr::rename(site = siteIndex)
    
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
   dplyr::left_join(uf_p.med, by="site") %>%
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
 
spatialUF <-uf_p %>%
   dplyr::select(-site) %>%
   dplyr::rename(site = siteIndex) %>%
   dplyr::left_join(position,by="site") %>%
   ggplot(aes(x = easting , y = med)) + 
   geom_point() +
   geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +
   theme_bw() +
   theme(panel.spacing=unit(0,"pt"), 
         panel.border=element_rect(colour="grey50", fill=NA)) +
   ylab("Fruits per plant (undamaged fruits)") +
   xlab("Easting (km)") +
   labs(title="Undamaged fruits per plant (2013-2018)") 

 ################################################################################
 # Posteriors for damaged fruits (2013-)
 #################################################################################

countUndamagedDamagedFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countUndamagedDamagedFruitsPerPlantAllPlots.RDS")

siteNames = unique(countUndamagedDamagedFruitsPerPlantAllPlots$site)

# recover site indices and year indices
siteIndex <- data.frame(siteIndex=unique(countUndamagedDamagedFruitsPerPlantAllPlots$site),site=unique(data$site2))
yearIndex <- data.frame(yearIndex=unique(countUndamagedDamagedFruitsPerPlantAllPlots$year),
                        year=unique(data$year2)) 

# extract population level distributions for undamaged plants
df_p <- mcmcSamples %>%
   tidybayes::recover_types(data) %>%
   tidybayes::spread_draws(mu_p_dam[site]) %>%
   dplyr::group_by(site) %>%
   dplyr::summarise(med = median(mu_p_dam,0), 
                    ci.lo = HPDI(mu_p_dam,0.89)[1], 
                    ci.hi = HPDI(mu_p_dam,0.89)[2],
                    ci.lo2 = HPDI(mu_p_dam,0.5)[1], 
                    ci.hi2 = HPDI(mu_p_dam,0.5)[2]
   )

df_p <- cbind(df_p,siteIndex=siteNames)

# df_py<-mcmcSamples %>%
#    tidybayes::recover_types(data) %>%
#    tidybayes::spread_draws(mu_py_dam[site,year]) %>%
#    dplyr::group_by(site,year) %>%
#    dplyr::summarise(med = median(mu_py_dam), 
#                     ci.lo = HPDI(mu_py_dam,0.95)[1], 
#                     ci.hi = HPDI(mu_py_dam,0.95)[2],
#                     ci.lo2 = HPDI(mu_py_dam,0.5)[1], 
#                     ci.hi2 = HPDI(mu_py_dam,0.5)[2]
#    )
# 
# df_py <- df_py %>%
#    dplyr::mutate(year = as.integer(year)) %>%
#    dplyr::left_join(siteIndex,by="site") %>%
#    dplyr::left_join(yearIndex,by="year") %>%
#    dplyr::ungroup() %>%
#    dplyr::select(-c(year)) %>%
#    dplyr::select(-c(site)) %>%
#    dplyr::rename(site = siteIndex) %>%
#    dplyr::rename(year = yearIndex) %>%
#    dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)

# this needs to be fixed
df_p.med = df_p %>%
   dplyr::rename(lambda.med=med) %>%
   dplyr::select(site,lambda.med) %>%
   dplyr::left_join(siteIndex,by="site") %>%
   dplyr::select(-site) %>%
   dplyr::rename(site = siteIndex)

interannualDF <- mcmcSamples %>%
   tidybayes::recover_types(data) %>%
   tidybayes::spread_draws(mu_py_dam[site,year]) %>%
   dplyr::group_by(site,year) %>%
   dplyr::summarise(med = median(mu_py_dam),
                    ci.lo = HPDI(mu_py_dam,0.89)[1],
                    ci.hi = HPDI(mu_py_dam,0.89)[2]
   ) %>%
   dplyr::ungroup() %>%
   dplyr::left_join(siteIndex,by="site") %>%
   dplyr::left_join(yearIndex,by="year") %>%
   dplyr::select(-c(site,year)) %>%
   dplyr::rename(site = siteIndex) %>%
   dplyr::rename(year = yearIndex) %>%
   dplyr::left_join(df_p.med, by="site") %>%
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
   ylab("Fruits per plant (damaged fruits)") 

spatialDF <-df_p %>%
   dplyr::select(-site) %>%
   dplyr::rename(site = siteIndex) %>%
   dplyr::left_join(position,by="site") %>%
   ggplot(aes(x = easting , y = med)) + 
   geom_point() +
   geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +
   theme_bw() +
   theme(panel.spacing=unit(0,"pt"), 
         panel.border=element_rect(colour="grey50", fill=NA)) +
   ylab("Fruits per plant (damaged fruits)") +
   xlab("Easting (km)") +
   labs(title="Damaged fruits per plant (2013-2018)") 

################################################################################
# Total fruit equivalents (2006-2012)
#################################################################################

countFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countFruitsPerPlantAllPlots.RDS")

siteNames = unique(countFruitsPerPlantAllPlots$site)

# recover site indices and year indices
siteIndex <- data.frame(siteIndex=unique(countFruitsPerPlantAllPlots$site),site=unique(data$site))
yearIndex <- data.frame(yearIndex=unique(countFruitsPerPlantAllPlots$year),
                        year=unique(data$year)) 

# extract population level distributions for undamaged plants
tfe_p <- mcmcSamples %>%
   tidybayes::recover_types(data) %>%
   tidybayes::spread_draws(mu_p_tfe[site]) %>%
   dplyr::group_by(site) %>%
   dplyr::summarise(med = median(mu_p_tfe,0), 
                    ci.lo = HPDI(mu_p_tfe,0.89)[1], 
                    ci.hi = HPDI(mu_p_tfe,0.89)[2],
                    ci.lo2 = HPDI(mu_p_tfe,0.5)[1], 
                    ci.hi2 = HPDI(mu_p_tfe,0.5)[2]
   )

tfe_p <- cbind(tfe_p,siteIndex=siteNames)

# df_py<-mcmcSamples %>%
#    tidybayes::recover_types(data) %>%
#    tidybayes::spread_draws(mu_py_dam[site,year]) %>%
#    dplyr::group_by(site,year) %>%
#    dplyr::summarise(med = median(mu_py_dam), 
#                     ci.lo = HPDI(mu_py_dam,0.95)[1], 
#                     ci.hi = HPDI(mu_py_dam,0.95)[2],
#                     ci.lo2 = HPDI(mu_py_dam,0.5)[1], 
#                     ci.hi2 = HPDI(mu_py_dam,0.5)[2]
#    )
# 
# df_py <- df_py %>%
#    dplyr::mutate(year = as.integer(year)) %>%
#    dplyr::left_join(siteIndex,by="site") %>%
#    dplyr::left_join(yearIndex,by="year") %>%
#    dplyr::ungroup() %>%
#    dplyr::select(-c(year)) %>%
#    dplyr::select(-c(site)) %>%
#    dplyr::rename(site = siteIndex) %>%
#    dplyr::rename(year = yearIndex) %>%
#    dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)

# this needs to be fixed
tfe_p.med = tfe_p %>%
   dplyr::rename(lambda.med=med) %>%
   dplyr::select(site,lambda.med) %>%
   dplyr::left_join(siteIndex,by="site") %>%
   dplyr::select(-site) %>%
   dplyr::rename(site = siteIndex)


interannualTFE <- mcmcSamples %>%
   tidybayes::recover_types(data) %>%
   tidybayes::spread_draws(mu_py_tfe[site,year]) %>%
   dplyr::group_by(site,year) %>%
   dplyr::summarise(med = median(mu_py_tfe),
                    ci.lo = HPDI(mu_py_tfe,0.89)[1],
                    ci.hi = HPDI(mu_py_tfe,0.89)[2]
   ) %>%
   dplyr::ungroup() %>%
   dplyr::left_join(siteIndex,by="site") %>%
   dplyr::left_join(yearIndex,by="year") %>%
   dplyr::select(-c(site,year)) %>%
   dplyr::rename(site = siteIndex) %>%
   dplyr::rename(year = yearIndex) %>%
   dplyr::left_join(tfe_p.med, by="site") %>%
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
   ylab("Fruits per plant (total fruit equivalents)") 

spatialTFE <-tfe_p %>%
   dplyr::select(-site) %>%
   dplyr::rename(site = siteIndex) %>%
   dplyr::left_join(position,by="site") %>%
   ggplot(aes(x = easting , y = med)) + 
   geom_point() +
   geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +
   theme_bw() +
   theme(panel.spacing=unit(0,"pt"), 
         panel.border=element_rect(colour="grey50", fill=NA)) +
   ylab("Fruits per plant (total fruit equivalents)") +
   xlab("Easting (km)") +
   labs(title="Total fruit equivalents per plant (2006-2012)") 
 
################################################################################
# Total fruit equivalent composites (2013-2018)
#################################################################################

countSeedPerDamagedFruit <- readRDS(dataFiles[[2]])

siteNames = unique(countSeedPerDamagedFruit$site4) 

# recover site indices and year indices
siteIndex <- data.frame(siteIndex=unique(countSeedPerDamagedFruit$site4),site=unique(data$site4))
yearIndex <- data.frame(yearIndex=unique(countSeedPerDamagedFruit$year4),
                        year=unique(data$year4)) 

# extract population level distributions for undamaged plants
# tfe_p <- mcmcSamples %>%
#    tidybayes::recover_types(data) %>%
#    tidybayes::spread_draws(mu_p_tfe[site]) %>%
#    dplyr::group_by(site) %>%
#    dplyr::summarise(med = median(mu_p_tfe,0), 
#                     ci.lo = HPDI(mu_p_tfe,0.89)[1], 
#                     ci.hi = HPDI(mu_p_tfe,0.89)[2],
#                     ci.lo2 = HPDI(mu_p_tfe,0.5)[1], 
#                     ci.hi2 = HPDI(mu_p_tfe,0.5)[2]
#    )
# 
# tfe_p <- cbind(tfe_p,siteIndex=siteNames)

tfecomp_py<-mcmcSamples %>%
   tidybayes::recover_types(data) %>%
   tidybayes::spread_draws(mu_py_tfe_comp[site,year]) %>%
   dplyr::group_by(site,year) %>%
   dplyr::summarise(med = median(mu_py_tfe_comp),
                    ci.lo = HPDI(mu_py_tfe_comp,0.95)[1],
                    ci.hi = HPDI(mu_py_tfe_comp,0.95)[2],
                    ci.lo2 = HPDI(mu_py_tfe_comp,0.5)[1],
                    ci.hi2 = HPDI(mu_py_tfe_comp,0.5)[2]
   )

tfecomp_py <- tfecomp_py %>%
   dplyr::mutate(year = as.integer(year)) %>%
   dplyr::left_join(siteIndex,by="site") %>%
   dplyr::left_join(yearIndex,by="year") %>%
   dplyr::ungroup() %>%
   dplyr::select(-c(year)) %>%
   dplyr::select(-c(site)) %>%
   dplyr::rename(site = siteIndex) %>%
   dplyr::rename(year = yearIndex) %>%
   dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)


   
   tfeCompDF<- mcmcSamples %>%
   tidybayes::recover_types(data) %>%
   tidybayes::spread_draws(mu_py_tfe_comp[site,year]) %>%
   dplyr::group_by(site,year) %>%
   dplyr::summarise(med = median(mu_py_tfe_comp),
                    ci.lo = HPDI(mu_py_tfe_comp,0.89)[1],
                    ci.hi = HPDI(mu_py_tfe_comp,0.89)[2]
   ) %>%
   dplyr::ungroup() %>%
   dplyr::left_join(siteIndex,by="site") %>%
   dplyr::left_join(yearIndex,by="year") %>%
   dplyr::select(-c(site,year)) %>%
   dplyr::rename(site = siteIndex) %>%
   dplyr::rename(year = yearIndex) %>%
   # dplyr::left_join(tfe_p.med, by="site") %>%
   # dplyr::arrange(lambda.med) %>% 
   dplyr::mutate(site=factor(site,levels=unique(site))) %>%
   dplyr::group_by(site) %>%
   dplyr::arrange(desc(med),.by_group = TRUE) %>%
   dplyr::mutate(year=factor(year)) %>%
   mutate(id = row_number()) 

# manually filter out S22 2013
tfeCompDF<- tfeCompDF %>% 
   dplyr::filter(site!="S22"|year!=2013) %>%
   dplyr::filter(site!="OSR"|year!=2013) %>%
   dplyr::filter(site!="LCE"|year!=2013) %>%
   dplyr::filter(site!="GCN"|year!=2013)
   
# plot   
interannualTFEcomp <-  ggplot(tfeCompDF,aes(x = id , y = med)) + 
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
   ylab("Fruits per plant (total fruit equivalents)") 

################################################################################
# ALL TOTAL FRUIT EQUIVALENTS
#################################################################################

countFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countFruitsPerPlantAllPlots.RDS")

siteNames = unique(countFruitsPerPlantAllPlots$site)

# recover site indices and year indices
siteIndex <- data.frame(siteIndex=unique(countFruitsPerPlantAllPlots$site),site=unique(data$site))
yearIndex <- data.frame(yearIndex=unique(countFruitsPerPlantAllPlots$year),
                        year=unique(data$year)) 

tfeDF<-mcmcSamples %>%
   tidybayes::recover_types(data) %>%
   tidybayes::spread_draws(mu_py_tfe[site,year]) %>%
   dplyr::group_by(site,year) %>%
   dplyr::summarise(med = median(mu_py_tfe),
                    ci.lo = HPDI(mu_py_tfe,0.89)[1],
                    ci.hi = HPDI(mu_py_tfe,0.89)[2]
   ) %>%
   dplyr::ungroup() %>%
   dplyr::left_join(siteIndex,by="site") %>%
   dplyr::left_join(yearIndex,by="year") %>%
   dplyr::select(-c(site,year)) %>%
   dplyr::rename(site = siteIndex) %>%
   dplyr::rename(year = yearIndex)


countSeedPerDamagedFruit <- readRDS(dataFiles[[2]])

siteNames = unique(countSeedPerDamagedFruit$site4) 

# recover site indices and year indices
siteIndex <- data.frame(siteIndex=unique(countSeedPerDamagedFruit$site4),site=unique(data$site4))
yearIndex <- data.frame(yearIndex=unique(countSeedPerDamagedFruit$year4),
                        year=unique(data$year4)) 

tfeCompDF<- mcmcSamples %>%
   tidybayes::recover_types(data) %>%
   tidybayes::spread_draws(mu_py_tfe_comp[site,year]) %>%
   dplyr::group_by(site,year) %>%
   dplyr::summarise(med = median(mu_py_tfe_comp),
                    ci.lo = HPDI(mu_py_tfe_comp,0.89)[1],
                    ci.hi = HPDI(mu_py_tfe_comp,0.89)[2]
   ) %>%
   dplyr::ungroup() %>%
   dplyr::left_join(siteIndex,by="site") %>%
   dplyr::left_join(yearIndex,by="year") %>%
   dplyr::select(-c(site,year)) %>%
   dplyr::rename(site = siteIndex) %>%
   dplyr::rename(year = yearIndex) 

# manually filter out S22 2013
tfeCompDF<- tfeCompDF %>% 
   dplyr::filter(site!="S22"|year!=2013) %>%
   dplyr::filter(site!="OSR"|year!=2013) %>%
   dplyr::filter(site!="LCE"|year!=2013) %>%
   dplyr::filter(site!="GCN"|year!=2013)

reference<-tfeCompDF %>% dplyr::select(site,year) %>% unique()
reference$year <- as.integer(reference$year)
tfeCompDF$year<-as.integer(tfeCompDF$year)

# undamaged
countFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/countUndamagedDamagedFruitsPerPlantAllPlots.rds")
siteIndex <- data.frame(siteIndex=unique(countFruitsPerPlantAllPlots$site),site=unique(data$site))
yearIndex <- data.frame(yearIndex=unique(countFruitsPerPlantAllPlots$year),
                        year=unique(data$year)) 
ufDF<-mcmcSamples %>%
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
   dplyr::rename(year = yearIndex)

ufDF$year<-as.integer(ufDF$year)

referenceUF<-ufDF %>% dplyr::select(site,year) %>% unique()

ufOnly<-dplyr::setdiff(data.frame(referenceUF),data.frame(reference)) %>%
   dplyr::mutate(site.year=paste0(site,year))

# countUndamagedDamagedFruitsPerPlantAllPlots %>%
#    dplyr::mutate(site.year=paste0(site,year)) %>%
#    dplyr::filter(site.year %in%ufOnly$site.year) 

ufDF_missing<-ufDF %>%
   dplyr::mutate(site.year=paste0(site,year)) %>%
   dplyr::filter(site.year %in%ufOnly$site.year) %>% 
   dplyr::select(-site.year)

# bind all datasets
tfeFullDF<-bind_rows(tfeDF,tfeCompDF,ufDF_missing)

# fill in with NAs
referenceDF <- expand.grid(site=unique(tfeFullDF$site),year=unique(tfeFullDF$year))
tfeFullDF<-tfeFullDF %>%
   dplyr::full_join(referenceDF,by=c("site","year"))

GeomRibbon$handle_na <- function(data, params) {  data }

# time plot
ribbonTFEfull <-   tfeFullDF %>%
   ggplot(aes(x=year,y=med,group=site,ymin=ci.lo,ymax=ci.hi)) +
   geom_ribbon(alpha=0.5,colour="grey50") +
   geom_line() +
   theme_bw() +
   theme(panel.spacing=unit(0,"pt"), 
         panel.border=element_rect(colour="grey50", fill=NA)) +
   xlab("Year") +
   ylab("Total fruit equivalents per plant") +
   facet_wrap(~site,scales='free_y') +
   labs(title="Total fruit equivalents per plant, composite (2006-2018)") 
   

# time plot
   tfeFullDFsorted <-  tfeFullDF %>%
      dplyr::group_by(site) %>%
      dplyr::arrange(desc(med),.by_group = TRUE) %>%
      dplyr::mutate(year=factor(year)) %>%
      mutate(id = row_number())

   interannualTFEfull <-   ggplot(tfeFullDFsorted,aes(x = id , y = med)) + 
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
   ylab("Fruits per plant (total fruit equivalents)")  +
      labs(color="Year") 

   # spatial plot
   tfeFullDFsorted <-  tfeFullDFsorted %>%
      dplyr::left_join(position,by="site")
   
spatialTFEFull <-   ggplot(tfeFullDFsorted,aes(x = easting , y = med)) + 
      geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +  
      geom_point(size=1) +
      theme_bw() +
      facet_wrap( ~ year, scales="free_y",nrow=2) +
      theme(panel.spacing=unit(0,"pt"), 
            panel.border=element_rect(colour="grey50", fill=NA)) +
      xlab("Easting") +
      ylab("Fruits per plant (total fruit equivalents)")  +
      labs(color="Year")
     
## Save Figures
# ggsave(filename=paste0(dirFigures,"interannualUF.pdf"),
#        plot=interannualUF,width=6,height=16)
# ggsave(filename=paste0(dirFigures,"interannualDF.pdf"),
#        plot=interannualDF,width=6,height=12)
# ggsave(filename=paste0(dirFigures,"interannualTFE.pdf"),
#        plot=interannualTFE,width=6,height=12)
# ggsave(filename=paste0(dirFigures,"interannualTFEcomp.pdf"),
#        plot=interannualTFEcomp,width=6,height=12)
# ggsave(filename=paste0(dirFigures,"interannualTFEFull.pdf"),
#        plot=interannualTFEfull,width=6,height=12)
# 
# ggsave(filename=paste0(dirFigures,"spatialUF.pdf"),
#        plot=spatialUF,width=4,height=4)
# ggsave(filename=paste0(dirFigures,"spatialDF.pdf"),
#        plot=spatialDF,width=4,height=4)
# ggsave(filename=paste0(dirFigures,"spatialTFE.pdf"),
#        plot=spatialTFE,width=4,height=4)
# # ggsave(filename=paste0(dirFigures,"interannualTFEcomp.pdf"),
# #        plot=interannualTFEcomp,width=6,height=12)
# ggsave(filename=paste0(dirFigures,"spatialTFEFull.pdf"),
#        plot=spatialTFEFull,width=12,height=6)
# 
# ggsave(filename=paste0(dirFigures,"ribbonTFEfull.pdf"),
#        plot=ribbonTFEfull,width=12,height=6)
