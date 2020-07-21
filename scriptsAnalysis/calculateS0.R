# Calculate s0

library(tidyverse)

# get the file with fruits per plot from transects
countFruitsPerPlantTransects<-readRDS("~/Dropbox/dataLibrary/postProcessingData/countFruitsPerPlantTransects.RDS")

transectCounts <- countFruitsPerPlantTransects %>%
   dplyr::mutate(position=as.character(position)) %>%
   dplyr::group_by(site,year,transect,position) %>%
   dplyr::summarise(total=sum(countFruitsPerPlant))
 
 # get file with median number of seeds per fruit
phiSummary<-read.csv("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/parameterSummary/seedsSummary.csv")[,-1]

# join the two frames for fitness
seedOutput=transectCounts %>%
   dplyr::filter(year<2009) %>%
   dplyr::left_join(phiSummary,by=c("site","year")) %>%
   dplyr::mutate(seed_output = total*mu) %>%
   dplyr::select(site,year,transect,position,seed_output)

# get seedling counts
censusSeedlingsFruitingPlants <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedlingSurvival/censusSeedlingsFruitingPlants.rds")
censusSeedlingsFruitingPlants <- censusSeedlingsFruitingPlants %>%
   dplyr::mutate(year = as.integer(year))

# join seedling counts to fitness data frame
yearOverYear<-seedOutput %>%
   dplyr::left_join(censusSeedlingsFruitingPlants,by=c("site","transect","position","year")) 

#get belowgroundrates
s1Summary<-read.csv("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/parameterSummary/s1Summary.csv")[,-1]
g1Summary<-read.csv("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/parameterSummary/g1Summary.csv")[,-1]
s2Summary<-read.csv("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/parameterSummary/s2Summary.csv")[,-1]

s1Summary<-s1Summary %>%
   dplyr::select(site,mu) %>%
   dplyr::rename(mu_s1 = mu)

g1Summary<-g1Summary %>%
   dplyr::select(site,mu) %>%
   dplyr::rename(mu_g1 = mu)

tmp<-yearOverYear %>%
   dplyr::left_join(s1Summary %>% dplyr::select(site,mu_s1),by=c("site")) %>%
   dplyr::left_join(g1Summary %>% dplyr::select(site,mu_g1) ,by=c("site"))

s0<-tmp %>% 
   dplyr::ungroup() %>%
   dplyr::mutate(s0=seedlingNumber/(seed_output*mu_s1*mu_g1)) %>%
   dplyr::group_by(site) %>%
   dplyr::summarise(mu_s0 = mean(s0,na.rm=TRUE))

plot(s0$mu_s0);abline(h=1)

s0<-s2Summary %>% dplyr::select(site,mu) %>%
   dplyr::mutate(s0=mu^(3/8)) %>%
   dplyr::select(site,s0)

write.csv(s0,file=paste0(dirPars,"s0Summary.csv"))

# adjusted
   dplyr::left_join(s0,by="site") %>%
   dplyr::mutate(s0_adj = s0*mu_s0)

################################################################################
# Germination G1
#################################################################################
# directory = "/Users/Gregor/Dropbox/dataLibrary/posteriors/"
# simFiles <- paste0(directory,list.files(directory))
# 
# mcmcSamples <- readRDS(simFiles[[8]])
# 
# directory = "/Users/Gregor/Dropbox/dataLibrary/workflow/tidyData/"
# simFiles <- paste0(directory,list.files(directory))
# 
# data <- readRDS(simFiles[[1]])
# siteNames = unique(data$siteBags)
# 
# g1_py<-mcmcSamples %>%
#    tidybayes::recover_types(data) %>%
#    tidybayes::spread_draws(g1.0[siteBags,yearBags]) %>%
#    dplyr::group_by(siteBags,yearBags) %>%
#    dplyr::summarise(mu = mean(g1.0),
#                     ci.lo95 = quantile(g1.0,probs=0.025), 
#                     ci.hi95 = quantile(g1.0,probs=0.975),
#                     med = median(g1.0), 
#                     hpdi.lo95 = HPDI(g1.0,.95)[1], 
#                     hpdi.hi95 = HPDI(g1.0,.95)[2]
#    )
# 
# siteIndex <- data.frame(siteIndex=unique(data$siteBags),siteBags=1:20)
# yearIndex <- data.frame(yearBags=1:3,yearIndex=2006:2008)
# 
# g1_py<-g1_py %>%
#    dplyr::mutate(yearBags = as.integer(yearBags)) %>%
#    dplyr::left_join(yearIndex,by="yearBags") %>%
#    dplyr::ungroup() %>%
#    dplyr::select(-c(yearBags)) %>%
#    dplyr::rename(site = siteBags) %>%
#    dplyr::rename(year = yearIndex) 
# 
# 
# 
# ################################################################################
# # Winter seed survival (s1)
# #################################################################################
# s1_py<-mcmcSamples %>%
#    tidybayes::recover_types(data) %>%
#    tidybayes::spread_draws(s1.0[siteBags,yearBags]) %>%
#    dplyr::group_by(siteBags,yearBags) %>%
#    dplyr::summarise(mu = mean(s1.0),
#                     ci.lo95 = quantile(s1.0,probs=0.025), 
#                     ci.hi95 = quantile(s1.0,probs=0.975),
#                     med = median(s1.0), 
#                     hpdi.lo95 = HPDI(s1.0,.95)[1], 
#                     hpdi.hi95 = HPDI(s1.0,.95)[2]
#    )
# 
# s1_py<-s1_py %>%
#    dplyr::mutate(yearBags = as.integer(yearBags)) %>%
#    dplyr::left_join(yearIndex,by="yearBags") %>%
#    dplyr::ungroup() %>%
#    dplyr::select(-c(yearBags)) %>%
#    dplyr::rename(site = siteBags) %>%
#    dplyr::rename(year = yearIndex)
# 
# 
# ## calculate annual
# 
# tmp<-yearOverYear %>%
#    dplyr::left_join(s1_py %>% dplyr::select(site,year,mu) %>% dplyr::rename(mu_s1 = mu),by=c("site","year")) %>%
#    dplyr::left_join(g1_py %>% dplyr::select(site,year,mu) %>% dplyr::rename(mu_g1 = mu),by=c("site","year"))
# 
# s0<-tmp %>% 
#    dplyr::ungroup() %>%
#    dplyr::mutate(s0=seedlingNumber/(seed_output*mu_s1*mu_g1)) %>%
#    dplyr::group_by(site,year) %>%
#    dplyr::summarise(mu_s0 = mean(s0,na.rm=TRUE))
# 
# plot(s0$mu_s0);abline(h=1)
