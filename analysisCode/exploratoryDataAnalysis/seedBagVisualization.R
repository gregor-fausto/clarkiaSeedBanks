# script to follow data visualization workflow from Gabry
# with seed bag data
# rainclouds from Allen et al. 2019

source("/Users/Gregor/Dropbox/clarkiaSeedBanks/library/geomFlatViolin.R")
load(file="/Users/Gregor/Dropbox/clarkiaSeedBanks/library/dataFromWorkFlowFile/seedBagsData.rda")

library(ggplot2)

age1 <- seedBags %>%
  dplyr::filter(age==1)

ggplot(data=seedBags) +
  geom_point(aes(x=site,y=seedlingJan))

library(cowplot) 
w=6
h=3

# Seedlings in January
seedBagsSummary <- seedBags %>%
  dplyr::group_by(site,yearData) %>%
  dplyr::summarise(seedlingsMean = mean(seedlingJan,na.rm=TRUE),
                   se = sd(seedlingJan,na.rm=TRUE)/sqrt(n())) %>%
  dplyr::mutate(yearData = as.factor(yearData))

globalMean <- mean(seedBags$seedlingJan,na.rm=TRUE)

pSeedlingsJan <- ggplot(seedBags %>%   dplyr::mutate(yearData = as.factor(yearData)),aes(x=yearData,y=seedlingJan,fill=yearData,color=yearData))+
  facet_wrap(~site) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim=TRUE)+
  geom_point(position = position_jitter(width = .15), size = .25)+
  geom_boxplot(aes(x = as.numeric(as.factor(yearData))+0.25, y = seedlingJan),outlier.shape = NA,
               alpha = 0.3, width = .1, colour = "BLACK") +
  guides(fill=FALSE, color=FALSE) +
  ylab('seedlings')+xlab('site')+theme_cowplot()+
  ggtitle('Figure: Seedlings in January, Age 1') +
  geom_hline(yintercept=globalMean,linetype="dashed")
pSeedlingsJan

pSeedlingsJan2 <- ggplot(seedBags,aes(x=site,y=seedlingJan,fill=site,color=site))+
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim=TRUE)+
  geom_point(position = position_jitter(width = .15), size = .25)+
  geom_boxplot(aes(x = as.numeric(as.factor(site))+0.25, y = seedlingJan),outlier.shape = NA,
               alpha = 0.3, width = .1, colour = "BLACK") +
  coord_flip() + guides(fill=FALSE, color=FALSE) +
  ylab('seedlings')+xlab('site')+theme_cowplot()+
  ggtitle('Figure: Seedlings in January, Age 1') +
  geom_hline(yintercept=globalMean,linetype="dashed")
pSeedlingsJan2

# seeds in January
seedBagsSummary <- seedBags %>%
  dplyr::group_by(site,yearData) %>%
  dplyr::summarise(seedsMean = mean(intactJan,na.rm=TRUE),
                   se = sd(intactJan,na.rm=TRUE)/sqrt(n())) 

globalMean <- mean(seedBags$intactJan,na.rm=TRUE)

pIntactJan <- ggplot(seedBags %>%   dplyr::mutate(yearData = as.factor(yearData)),aes(x=yearData,y=intactJan,fill=yearData,color=yearData))+
  facet_wrap(~site) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim=TRUE)+
  geom_point(position = position_jitter(width = .15), size = .25)+
  geom_boxplot(aes(x = as.numeric(as.factor(yearData))+0.25, y = intactJan),outlier.shape = NA,
               alpha = 0.3, width = .1, colour = "BLACK") +
  guides(fill=FALSE, color=FALSE) +
  ylab('seedlings')+xlab('site')+theme_cowplot()+
  ggtitle('Figure: Intact seeds in January, Age 1') +
  geom_hline(yintercept=globalMean,linetype="dashed")
pIntactJan

pIntactJan2 <- ggplot(seedBags,aes(x=site,y=intactJan,fill=site,color=site))+
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim=TRUE)+
  geom_point(position = position_jitter(width = .15), size = .25)+
  geom_boxplot(aes(x = as.numeric(as.factor(site))+0.25, y = intactJan),outlier.shape = NA,
               alpha = 0.3, width = .1, colour = "BLACK") +
  coord_flip() + guides(fill=FALSE, color=FALSE) +
  ylab('seedlings')+xlab('site')+theme_cowplot()+
  ggtitle('Figure: Intact seeds in January, Age 1') +
  geom_hline(yintercept=globalMean,linetype="dashed")
pIntactJan2


# seeds in October
seedBagsSummary <- seedBags %>%
  dplyr::group_by(site,yearData) %>%
  dplyr::summarise(seedsMean = mean(intactOct,na.rm=TRUE),
                   se = sd(intactOct,na.rm=TRUE)/sqrt(n())) 

globalMean <- mean(seedBags$intactOct,na.rm=TRUE)

pIntactOct <- ggplot(seedBags %>%   dplyr::mutate(yearData = as.factor(yearData)),aes(x=yearData,y=intactOct,fill=yearData,color=yearData))+
  facet_wrap(~site) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim=TRUE)+
  geom_point(position = position_jitter(width = .15), size = .25)+
  geom_boxplot(aes(x = as.numeric(as.factor(yearData))+0.25, y = intactOct),outlier.shape = NA,
               alpha = 0.3, width = .1, colour = "BLACK") +
  guides(fill=FALSE, color=FALSE) +
  ylab('seedlings')+xlab('site')+theme_cowplot()+
  ggtitle('Figure: Intact seeds in October, Age 1') +
  geom_hline(yintercept=globalMean,linetype="dashed")
pIntactOct

pIntactOct2 <- ggplot(seedBags,aes(x=site,y=intactOct,fill=site,color=site))+
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim=TRUE)+
  geom_point(position = position_jitter(width = .15), size = .25)+
  geom_boxplot(aes(x = as.numeric(as.factor(site))+0.25, y = intactOct),outlier.shape = NA,
               alpha = 0.3, width = .1, colour = "BLACK") +
  coord_flip() + guides(fill=FALSE, color=FALSE) +
  ylab('seedlings')+xlab('site')+theme_cowplot()+
  ggtitle('Figure: Intact seeds in October, Age 1') +
  geom_hline(yintercept=globalMean,linetype="dashed")
pIntactOct2
# 
# 
# p1 <- ggplot(seedBagsSummary, aes(x = site, y = seedlingsMean, fill = site))+
#   geom_bar(stat = "identity", width = .8)+
#   geom_errorbar(aes(ymin = seedlingsMean - se, ymax = seedlingsMean+se), width = .2)+
#   guides(fill=FALSE)+
#   ylim(0, 20)+
#   ylab('seedlings')+xlab('Site')+theme_cowplot()+
#   ggtitle("Figure R1: Barplot +/- SEM")
# 
# p1
# 
# p2 <- ggplot(seedBags,aes(x=site,y=seedlingJan))+
#   geom_flat_violin(position = position_nudge(x = .2, y = 0),adjust =2)+
#   geom_point(position = position_jitter(width = .15), size = .25)+
#   ylab('seedlings')+xlab('site')+theme_cowplot()+
#   ggtitle('Figure R2: Basic Rainclouds or Little Prince Plot')
# p2
# 
# p3 <- ggplot(seedBags,aes(x=site,y=seedlingJan,fill=site))+
#   geom_flat_violin(position = position_nudge(x = .2, y = 0),adjust =2)+
#   geom_point(position = position_jitter(width = .15), size = .25)+
#   coord_flip() + guides(fill=FALSE) +
#   ylab('seedlings')+xlab('site')+theme_cowplot()+
#   ggtitle('Figure R4: Basic Rainclouds with Color')
# p3
# 
# p4 <- ggplot(seedBags,aes(x=site,y=seedlingJan,fill=site))+
#   geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2)+
#   geom_point(position = position_jitter(width = .15), size = .25)+
#   geom_boxplot(aes(x = as.numeric(as.factor(site))+0.25, y = seedlingJan),outlier.shape = NA,
#                alpha = 0.3, width = .1, colour = "BLACK") +
#   coord_flip() + guides(fill=FALSE, color=FALSE) +
#   ylab('seedlings')+xlab('site')+theme_cowplot()+
#   ggtitle('Figure R4:  Rainclouds Plot with Boxplots')
# p4
# 
# p5 <- ggplot(seedBags,aes(x=site,y=seedlingJan,fill=site,color=site))+
#   geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim=FALSE)+
#   geom_point(position = position_jitter(width = .15), size = .25)+
#   geom_boxplot(aes(x = as.numeric(as.factor(site))+0.25, y = seedlingJan),outlier.shape = NA,
#                alpha = 0.3, width = .1, colour = "BLACK") +
#   coord_flip() + guides(fill=FALSE, color=FALSE) +
#   ylab('seedlings')+xlab('site')+theme_cowplot()+
#   ggtitle('Figure R5:  Rainclouds Plot with Boxplots') +
#   geom_hline(yintercept=globalMean,linetype="dashed")
# p5
