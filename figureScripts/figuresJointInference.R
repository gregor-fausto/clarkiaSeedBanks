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

mcmcSamples <- readRDS(simFiles[[8]])
dirFigures = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/"

################################################################################
# Germination
#################################################################################

directory = "/Users/Gregor/Dropbox/dataLibrary/workflow/tidyData/"
simFiles <- paste0(directory,list.files(directory))

data <- readRDS(simFiles[[1]])
siteNames = unique(data$siteBags)

################################################################################
# Spatial data
#################################################################################
position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,northing,elevation) %>%
  dplyr::mutate(easting=easting/1000,northing=northing/1000)


################################################################################
# Germination 1
#################################################################################

g1_p <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(g1[siteBags]) %>%
  dplyr::group_by(siteBags) %>%
  dplyr::summarise(med = median(g1), 
                   ci.lo = HPDI(g1,.89)[1], 
                   ci.hi = HPDI(g1,.89)[2],
                   ci.lo2 = HPDI(g1,.5)[1],
                   ci.hi2 = HPDI(g1,.5)[2]
  )

g1_p<-cbind(g1_p,site=siteNames)

g1_p.med = g1_p %>% dplyr::select(site,med)

siteIndex <- data.frame(siteIndex=unique(data$siteBags),siteBags=1:20)
yearIndex <- data.frame(yearBags=1:3,yearIndex=2006:2008)

g1_py<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(g1.0[siteBags,yearBags]) %>%
  dplyr::group_by(siteBags,yearBags) %>%
  dplyr::summarise(med = median(g1.0), 
                   ci.lo = HPDI(g1.0,.89)[1], 
                   ci.hi = HPDI(g1.0,.89)[2],
                   ci.lo2 = HPDI(g1.0,.5)[1],
                   ci.hi2 = HPDI(g1.0,.5)[2]
  )

g1_py<-g1_py %>%
  dplyr::mutate(yearBags = as.integer(yearBags)) %>%
  dplyr::left_join(yearIndex,by="yearBags") %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(yearBags)) %>%
  dplyr::rename(site = siteBags) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)


siteIndex <- data.frame(siteIndex=unique(data$siteBags),site=1:20)
yearIndex <- data.frame(year=1:3,yearIndex=2006:2008)

interannualG1<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(g1.0[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(lambda.med = median(g1.0), 
                   ci.lo = HPDI(g1.0,.89)[1], 
                   ci.hi = HPDI(g1.0,.89)[2]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  dplyr::select(-c(site,year)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::left_join(g1_p.med, by="site") %>%
  dplyr::arrange(med) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(desc(lambda.med),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  mutate(id = row_number()) %>% 
  ggplot(aes(x = id , y = lambda.med)) + 
  geom_point(aes(color=year)) +
  geom_hline(aes(yintercept=med),linetype='dotted') +
  geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +
  coord_flip() +
  facet_grid(site ~ ., scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ylab("Probability of germination [P(G)]") +
  ylim(c(0,1))

spatialG1 <-g1_p %>%
  dplyr::left_join(position,by="site") %>%
  ggplot(aes(x = easting , y = med)) + 
  geom_point() +
  geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  ylab("[P(G)]") +
  xlab("Easting (km)") +
  ylim(c(0,1))


################################################################################
# Winter seed survival (s1)
#################################################################################

s1_p <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(s1[siteBags]) %>%
  dplyr::group_by(siteBags) %>%
  dplyr::summarise(med = median(s1), 
                   ci.lo = HPDI(s1,.89)[1], 
                   ci.hi = HPDI(s1,.89)[2],
                   ci.lo2 = HPDI(s1,.5)[1],
                   ci.hi2 = HPDI(s1,.5)[2]
  )

s1_p<-cbind(s1_p,site=siteNames)

s1_p.med = s1_p %>% dplyr::select(site,med)

siteIndex <- data.frame(siteIndex=unique(data$siteBags),siteBags=1:20)
yearIndex <- data.frame(yearBags=1:3,yearIndex=2006:2008)

s1_py<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(s1.0[siteBags,yearBags]) %>%
  dplyr::group_by(siteBags,yearBags) %>%
  dplyr::summarise(med = median(s1.0), 
                   ci.lo = HPDI(s1.0,.89)[1], 
                   ci.hi = HPDI(s1.0,.89)[2],
                   ci.lo2 = HPDI(s1.0,.5)[1],
                   ci.hi2 = HPDI(s1.0,.5)[2]
  )

s1_py<-s1_py %>%
  dplyr::mutate(yearBags = as.integer(yearBags)) %>%
  dplyr::left_join(yearIndex,by="yearBags") %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(yearBags)) %>%
  dplyr::rename(site = siteBags) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)


siteIndex <- data.frame(siteIndex=unique(data$siteBags),site=1:20)
yearIndex <- data.frame(year=1:3,yearIndex=2006:2008)

interannualS1<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(s1.0[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(lambda.med = median(s1.0), 
                   ci.lo = HPDI(s1.0,.89)[1], 
                   ci.hi = HPDI(s1.0,.89)[2]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  dplyr::select(-c(site,year)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::left_join(s1_p.med, by="site") %>%
  dplyr::arrange(med) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(desc(lambda.med),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  mutate(id = row_number()) %>% 
  ggplot(aes(x = id , y = lambda.med)) + 
  geom_point(aes(color=year)) +
  geom_hline(aes(yintercept=med),linetype='dotted') +
  geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +
  coord_flip() +
  facet_grid(site ~ ., scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ylab("Probability of seed survival in first winter [P(S1)]") +
  ylim(c(0,1))

spatialS1 <-s1_p %>%
  dplyr::left_join(position,by="site") %>%
  ggplot(aes(x = easting , y = med)) + 
  geom_point() +
  geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  ylab("[P(S1)]") +
  xlab("Easting (km)")+
  ylim(c(0,1))

# 
# ################################################################################
# # Summer seed survival (s2)
# #################################################################################

s2_p <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(s2[siteBags]) %>%
  dplyr::group_by(siteBags) %>%
  dplyr::summarise(med = median(s2), 
                   ci.lo = HPDI(s2,.89)[1], 
                   ci.hi = HPDI(s2,.89)[2],
                   ci.lo2 = HPDI(s2,.5)[1],
                   ci.hi2 = HPDI(s2,.5)[2]
  )

s2_p<-cbind(s2_p,site=siteNames)

s2_p.med = s2_p %>% dplyr::select(site,med)

siteIndex <- data.frame(siteIndex=unique(data$siteBags),siteBags=1:20)
yearIndex <- data.frame(yearBags=1:3,yearIndex=2006:2008)

s2_py<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(s2.0[siteBags,yearBags]) %>%
  dplyr::group_by(siteBags,yearBags) %>%
  dplyr::summarise(med = median(s2.0), 
                   ci.lo = HPDI(s2.0,.89)[1], 
                   ci.hi = HPDI(s2.0,.89)[2],
                   ci.lo2 = HPDI(s2.0,.5)[1],
                   ci.hi2 = HPDI(s2.0,.5)[2]
  )

s2_py<-s2_py %>%
  dplyr::mutate(yearBags = as.integer(yearBags)) %>%
  dplyr::left_join(yearIndex,by="yearBags") %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(yearBags)) %>%
  dplyr::rename(site = siteBags) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)


siteIndex <- data.frame(siteIndex=unique(data$siteBags),site=1:20)
yearIndex <- data.frame(year=1:3,yearIndex=2006:2008)

interannualS2<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(s2.0[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(lambda.med = median(s2.0), 
                   ci.lo = HPDI(s2.0,.89)[1], 
                   ci.hi = HPDI(s2.0,.89)[2]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  dplyr::select(-c(site,year)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::left_join(s2_p.med, by="site") %>%
  dplyr::arrange(med) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(desc(lambda.med),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  mutate(id = row_number()) %>% 
  ggplot(aes(x = id , y = lambda.med)) + 
  geom_point(aes(color=year)) +
  geom_hline(aes(yintercept=med),linetype='dotted') +
  geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +
  coord_flip() +
  facet_grid(site ~ ., scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ylab("Probability of seed survival in first winter [P(S2)]") +
  ylim(c(0,1))

spatialS2 <-s2_p %>%
  dplyr::left_join(position,by="site") %>%
  ggplot(aes(x = easting , y = med)) + 
  geom_point() +
  geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  ylab("[P(S2)]") +
  xlab("Easting (km)") +
  ylim(c(0,1))
# ################################################################################
# # Winter seed survival, year 2 (s3)
# #################################################################################

s3_p <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(s3[siteBags]) %>%
  dplyr::group_by(siteBags) %>%
  dplyr::summarise(med = median(s3), 
                   ci.lo = HPDI(s3,.89)[1], 
                   ci.hi = HPDI(s3,.89)[2],
                   ci.lo2 = HPDI(s3,.5)[1],
                   ci.hi2 = HPDI(s3,.5)[2]
  )

s3_p<-cbind(s3_p,site=siteNames)

s3_p.med = s3_p %>% dplyr::select(site,med)

siteIndex <- data.frame(siteIndex=unique(data$siteBags),siteBags=1:20)
yearIndex <- data.frame(yearBags=1:2,yearIndex=2006:2007)

s3_py<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(s3.0[siteBags,yearBags]) %>%
  dplyr::group_by(siteBags,yearBags) %>%
  dplyr::summarise(med = median(s3.0), 
                   ci.lo = HPDI(s3.0,.89)[1], 
                   ci.hi = HPDI(s3.0,.89)[2],
                   ci.lo2 = HPDI(s3.0,.5)[1],
                   ci.hi2 = HPDI(s3.0,.5)[2]
  )

s3_py<-s3_py %>%
  dplyr::mutate(yearBags = as.integer(yearBags)) %>%
  dplyr::left_join(yearIndex,by="yearBags") %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(yearBags)) %>%
  dplyr::rename(site = siteBags) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)


siteIndex <- data.frame(siteIndex=unique(data$siteBags),site=1:20)
yearIndex <- data.frame(year=1:2,yearIndex=2006:2007)

interannualS3<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(s3.0[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(lambda.med = median(s3.0), 
                   ci.lo = HPDI(s3.0,.89)[1], 
                   ci.hi = HPDI(s3.0,.89)[2]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  dplyr::select(-c(site,year)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::left_join(s3_p.med, by="site") %>%
  dplyr::arrange(med) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(desc(lambda.med),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  mutate(id = row_number()) %>% 
  ggplot(aes(x = id , y = lambda.med)) + 
  geom_point(aes(color=year)) +
  geom_hline(aes(yintercept=med),linetype='dotted') +
  geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +
  coord_flip() +
  facet_grid(site ~ ., scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ylab("Probability of seed survival in first winter [P(S3)]") +
  ylim(c(0,1))

spatialS3 <-s3_p %>%
  dplyr::left_join(position,by="site") %>%
  ggplot(aes(x = easting , y = med)) + 
  geom_point() +
  geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  ylab("[P(S3)]") +
  xlab("Easting (km)") +
  ylim(c(0,1))

# combine plots



interannualS1.j <- interannualS1 + theme(legend.position="none") + scale_fill_manual(values=c("#1b9e77", "#d95f02", "#7570b3"))
interannualG1.j <- interannualG1 + theme(legend.position="none") + scale_fill_manual(values=c("#1b9e77", "#d95f02", "#7570b3"))
interannualS2.j <- interannualS2 + theme(legend.position="none") + scale_fill_manual(values=c("#1b9e77", "#d95f02", "#7570b3"))
interannualS3.j <- interannualS3 + theme(legend.position="none") + scale_fill_manual(values=c("#1b9e77", "#d95f02", "#7570b3"))

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend= get_legend(interannualS1 + scale_fill_manual(values=c("#1b9e77", "#d95f02", "#7570b3"),name="Year"))

interannualPlot <- gridExtra::grid.arrange(interannualS1.j,interannualG1.j,interannualS2.j,interannualS3.j,nrow=1,legend)

ggsave(filename=paste0(dirFigures,"interannualBelowground.pdf"),
       plot=interannualPlot,width=20,height=10)



spatialS1.j <- spatialS1 + theme(legend.position="none") +theme(axis.title.x=element_blank())
spatialG1.j <- spatialG1 + theme(legend.position="none") +theme(axis.title.x=element_blank())
spatialS2.j <- spatialS2 + theme(legend.position="none") +theme(axis.title.x=element_blank())
spatialS3.j <- spatialS3 + theme(legend.position="none")

spatialPlot <- gridExtra::grid.arrange(spatialS1.j,spatialG1.j,spatialS2.j,spatialS3.j,ncol=1)

ggsave(filename=paste0(dirFigures,"spatialBelowground.pdf"),
       plot=spatialPlot,width=5,height=10)
