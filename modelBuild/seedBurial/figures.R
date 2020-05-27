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

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBurial/"
simFiles <- paste0(directory,list.files(directory))
dirFigures = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/")

mcmcSamples <- readRDS(simFiles[[3]])

################################################################################
# Checks
#################################################################################

parsToMonitor_1 = c("theta_1","mu0_1","sigma0_1","mu_1","sigma_1","p_1")
parsToMonitor_2 = c("theta_2","mu0_2","sigma0_2","mu_2","sigma_2","p_2")
parsToMonitor_3 = c("theta_3","mu0_3","sigma0_3","mu_3","sigma_3","p_3")
parsToMonitor_g = c("theta_g","mu0_g","sigma0_g","mu_g","sigma_g","p_g")
parsToMonitor_v = c("theta_v","mu0_v","sigma0_v","mu_v","sigma_v","p_v")
parsToMonitor_deriv = c("nu_1","s1","g1","s2")


library(bayesplot)

################################################################################
# Germination
#################################################################################

data <- readRDS(simFiles[[1]])
seedBagExperiment <- readRDS(simFiles[[2]])
siteNames = unique(seedBagExperiment$siteBags)

position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,northing,elevation) %>%
  dplyr::mutate(easting=easting/1000,northing=northing/1000)

################################################################################
# Germination
#################################################################################

mcmcSummary<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(g1[siteBags]) %>%
  dplyr::group_by(siteBags) %>%
  dplyr::summarise(med = median(g1), 
                   ci.lo = quantile(g1,probs=0.025), 
                   ci.hi = quantile(g1,probs=0.975),
                   ci.lo2 = quantile(g1,probs=0.25), 
                   ci.hi2 = quantile(g1,probs=0.75)
  )
mcmcSummary<-cbind(mcmcSummary,site=siteNames)

g1Summary <- mcmcSummary %>%
  dplyr::rename(ci.lo95 = ci.lo,
                ci.hi95 = ci.hi,
                ci.lo50 = ci.lo2,
                ci.hi50 = ci.hi2)

fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/parameterSummary/")
dir.create(file.path(fileDirectory), showWarnings = FALSE)
saveRDS(g1Summary,file=paste0(fileDirectory,"g1Summary.RDS"))


vr = position %>% 
  dplyr::left_join(mcmcSummary,by="site")

pdf(paste0(dirFigures,"spatial-g1.pdf"), width=8, height=6)

par(mfrow=c(1,1))
plot(vr$easting,vr$med,ylim=c(0.0,1),
     ylab="Germination probability [P(G)]",
     xlab="Easting (km)", 
     cex.lab = 1.5, cex.axis = 1.5,
     pch=16,type='n')
# abline(h=.05,lwd=.5,col='gray')
# abline(h=.1,lwd=.5,col='gray')
# abline(h=.15,lwd=.5,col='gray')
# abline(h=.2,lwd=.5,col='gray')
# abline(h=.25,lwd=.5,col='gray')
points(vr$easting,vr$med,ylim=c(0.0,1),
     pch=16)

segments(x0=vr$easting, y0=vr$ci.lo, y1=vr$ci.hi)
segments(x0=vr$easting, y0=vr$ci.lo2, y1=vr$ci.hi2,lwd=2)

dev.off()

g1<-MCMCchains(mcmcSamples,params="g1")

plot(density(g1),type='n',ylim=c(0,12))
for(i in 1:20){
  lines(density(g1[,i]),lwd=.5)
}


siteIndex <- data.frame(siteIndex=unique(seedBagExperiment$siteBags),siteBags=1:20)
yearIndex <- data.frame(yearBags=1:3,yearIndex=2006:2008)

mcmcSummary<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(g1.0[siteBags,yearBags]) %>%
  dplyr::group_by(siteBags,yearBags) %>%
  dplyr::summarise(med = median(g1.0), 
                   ci.lo = quantile(g1.0,probs=0.025), 
                   ci.hi = quantile(g1.0,probs=0.975),
                   ci.lo2 = quantile(g1.0,probs=0.25), 
                   ci.hi2 = quantile(g1.0,probs=0.75)
  )
mcmcSummary<-mcmcSummary %>%
  dplyr::left_join(siteIndex,by="siteBags") %>%
  dplyr::left_join(yearIndex,by="yearBags") %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(siteBags,yearBags)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)

vr = position %>% 
  dplyr::left_join(mcmcSummary,by="site")


tmp <- vr %>% 
  dplyr::select(site,year,med,ci.lo,ci.hi)

g1Summary <- tmp
write.csv(g1Summary , 
          file = "~/Dropbox/clarkiaSeedBanks/products/dataFiles/g1Summary.csv",
          row.names=FALSE)

g1 <- ggplot(data=vr) +
  geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +
  geom_linerange(aes(x=easting,ymin=ci.lo2,ymax=ci.hi2),size=1) +
  geom_point(aes(x=easting,y=med)) +
  facet_wrap(~year,nrow=1) +
  theme_bw() +
  ylab("Germination probability [P(G)]") +
  xlab("Easting (km)") 

g1

ggsave(filename=paste0(dirFigures,"annual-g1.pdf"),
       plot=g1,width=20,height=4)

library(pdftools)
pdf_combine(c(paste0(dirFigures,"spatial-g1.pdf"),
              paste0(dirFigures,"annual-g1.pdf")), 
            output = paste0(dirFigures,"summary-g1.pdf"))


################################################################################
# Winter seed survival (s1)
#################################################################################

mcmcSummary<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(s1[siteBags]) %>%
  dplyr::group_by(siteBags) %>%
  dplyr::summarise(med = median(s1), 
                   ci.lo = quantile(s1,probs=0.025), 
                   ci.hi = quantile(s1,probs=0.975),
                   ci.lo2 = quantile(s1,probs=0.25), 
                   ci.hi2 = quantile(s1,probs=0.75)
  )
mcmcSummary<-cbind(mcmcSummary,site=siteNames)


s1Summary <- mcmcSummary %>%
  dplyr::rename(ci.lo95 = ci.lo,
                ci.hi95 = ci.hi,
                ci.lo50 = ci.lo2,
                ci.hi50 = ci.hi2)

fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/parameterSummary/")
dir.create(file.path(fileDirectory), showWarnings = FALSE)
saveRDS(s1Summary,file=paste0(fileDirectory,"s1Summary.RDS"))


vr = position %>% 
  dplyr::left_join(mcmcSummary,by="site")

pdf(paste0(dirFigures,"spatial-s1.pdf"), width=8, height=6)

par(mfrow=c(1,1))
plot(vr$easting,vr$med,ylim=c(0.0,1),
     ylab="Seed survival probability in first winter [P(S1)]",
     xlab="Easting (km)", 
     cex.lab = 1.5, cex.axis = 1.5,
     pch=16,type='n')

points(vr$easting,vr$med,ylim=c(0.0,1),
       pch=16)

segments(x0=vr$easting, y0=vr$ci.lo, y1=vr$ci.hi)
segments(x0=vr$easting, y0=vr$ci.lo2, y1=vr$ci.hi2,lwd=2)

dev.off()

s1<-MCMCchains(mcmcSamples,params="s1")

plot(density(s1),type='n',ylim=c(0,8))
for(i in 1:20){
  lines(density(s1[,i]),lwd=.5)
}



siteIndex <- data.frame(siteIndex=unique(seedBagExperiment$siteBags),siteBags=1:20)
yearIndex <- data.frame(yearBags=1:3,yearIndex=2006:2008)

mcmcSummary<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(s1.0[siteBags,yearBags]) %>%
  dplyr::group_by(siteBags,yearBags) %>%
  dplyr::summarise(med = median(s1.0), 
                   ci.lo = quantile(s1.0,probs=0.025), 
                   ci.hi = quantile(s1.0,probs=0.975),
                   ci.lo2 = quantile(s1.0,probs=0.25), 
                   ci.hi2 = quantile(s1.0,probs=0.75)
  )
mcmcSummary<-mcmcSummary %>%
  dplyr::left_join(siteIndex,by="siteBags") %>%
  dplyr::left_join(yearIndex,by="yearBags") %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(siteBags,yearBags)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)

vr = position %>% 
  dplyr::left_join(mcmcSummary,by="site")


tmp <- vr %>% 
  dplyr::select(site,year,med,ci.lo,ci.hi)

s1Summary <- tmp
write.csv(s1Summary , 
          file = "~/Dropbox/clarkiaSeedBanks/products/dataFiles/s1Summary.csv",
          row.names=FALSE)

g1 <- ggplot(data=vr) +
  geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +
  geom_linerange(aes(x=easting,ymin=ci.lo2,ymax=ci.hi2),size=1) +
  geom_point(aes(x=easting,y=med)) +
  facet_wrap(~year,nrow=1) +
  theme_bw() +
  ylab("Seed survival probability in first winter [P(S1)]") +
  xlab("Easting (km)") 


ggsave(filename=paste0(dirFigures,"annual-s1.pdf"),
       plot=g1,width=20,height=4)

library(pdftools)
pdf_combine(c(paste0(dirFigures,"spatial-s1.pdf"),
              paste0(dirFigures,"annual-s1.pdf")), 
            output = paste0(dirFigures,"summary-s1.pdf"))


################################################################################
# Summer seed survival (s2)
#################################################################################

mcmcSummary<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(s2[siteBags]) %>%
  dplyr::group_by(siteBags) %>%
  dplyr::summarise(med = median(s2), 
                   ci.lo = quantile(s2,probs=0.025), 
                   ci.hi = quantile(s2,probs=0.975),
                   ci.lo2 = quantile(s2,probs=0.25), 
                   ci.hi2 = quantile(s2,probs=0.75)
  )
mcmcSummary<-cbind(mcmcSummary,site=siteNames)


s2Summary <- mcmcSummary %>%
  dplyr::rename(ci.lo95 = ci.lo,
                ci.hi95 = ci.hi,
                ci.lo50 = ci.lo2,
                ci.hi50 = ci.hi2)

fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/parameterSummary/")
dir.create(file.path(fileDirectory), showWarnings = FALSE)
saveRDS(s2Summary,file=paste0(fileDirectory,"s2Summary.RDS"))

vr = position %>% 
  dplyr::left_join(mcmcSummary,by="site")

pdf(paste0(dirFigures,"spatial-s2.pdf"), width=8, height=6)

par(mfrow=c(1,1))
plot(vr$easting,vr$med,ylim=c(0.0,1),
     ylab="Seed survival probability in first summer [P(S2)]",
     xlab="Easting (km)", 
     cex.lab = 1.5, cex.axis = 1.5,
     pch=16,type='n')

points(vr$easting,vr$med,ylim=c(0.0,1),
       pch=16)

segments(x0=vr$easting, y0=vr$ci.lo, y1=vr$ci.hi)
segments(x0=vr$easting, y0=vr$ci.lo2, y1=vr$ci.hi2,lwd=2)

dev.off()

s2<-MCMCchains(mcmcSamples,params="s2")

plot(density(s2),type='n',ylim=c(0,8))
for(i in 1:20){
  lines(density(s2[,i]),lwd=.5)
}


siteIndex <- data.frame(siteIndex=unique(seedBagExperiment$siteBags),siteBags=1:20)
yearIndex <- data.frame(yearBags=1:3,yearIndex=2006:2008)

mcmcSummary<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(s2.0[siteBags,yearBags]) %>%
  dplyr::group_by(siteBags,yearBags) %>%
  dplyr::summarise(med = median(s2.0), 
                   ci.lo = quantile(s2.0,probs=0.025), 
                   ci.hi = quantile(s2.0,probs=0.975),
                   ci.lo2 = quantile(s2.0,probs=0.25), 
                   ci.hi2 = quantile(s2.0,probs=0.75)
  )
mcmcSummary<-mcmcSummary %>%
  dplyr::left_join(siteIndex,by="siteBags") %>%
  dplyr::left_join(yearIndex,by="yearBags") %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(siteBags,yearBags)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)

vr = position %>% 
  dplyr::left_join(mcmcSummary,by="site")


tmp <- vr %>% 
  dplyr::select(site,year,med,ci.lo,ci.hi)

s2Summary <- tmp
write.csv(s1Summary , 
          file = "~/Dropbox/clarkiaSeedBanks/products/dataFiles/s2Summary.csv",
          row.names=FALSE)

g1 <- ggplot(data=vr) +
  geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +
  geom_linerange(aes(x=easting,ymin=ci.lo2,ymax=ci.hi2),size=1) +
  geom_point(aes(x=easting,y=med)) +
  facet_wrap(~year,nrow=1) +
  theme_bw() +
  ylab("Seed survival probability in first summer [P(S2)]") +
  xlab("Easting (km)") 


ggsave(filename=paste0(dirFigures,"annual-s2.pdf"),
       plot=g1,width=20,height=4)

library(pdftools)
pdf_combine(c(paste0(dirFigures,"spatial-s2.pdf"),
              paste0(dirFigures,"annual-s2.pdf")), 
            output = paste0(dirFigures,"summary-s2.pdf"))


# ################################################################################
# # ???
# #################################################################################
# ref<-seedBagExperiment %>% 
#   dplyr::select(yearBags,siteBags) %>%
#   tidyr::unite(yearSite,c(yearBags,siteBags))
# 
# nYears = length(unique(seedBagExperiment$yearBags))
# 
# sampleYear <-sample(1:nYears,1)
# par(mfrow=c(1,3))
# for(i in sampleYear){
#   # filter data to year
#   tmp = seedBagExperiment %>% 
#     dplyr::filter(yearBags==unique(seedBagExperiment$yearBags)[i]) %>%
#     dplyr::mutate(observed.p = totalJan/seedStart) %>%
#     dplyr::arrange(observed.p)
#   
#   mcmcSummary<-mcmcSamples %>%
#     tidybayes::spread_draws(theta_1[yearBags]) %>%
#     dplyr::group_by(yearBags) %>%
#     dplyr::summarise(med = median(theta_1), 
#                      ci.lo = quantile(theta_1,probs=0.025), 
#                      ci.hi = quantile(theta_1,probs=0.975))
#   
# tmp<-tmp %>%
#   dplyr::left_join(mcmcSummary,by="yearBags")
# 
# mle = sum(tmp$fruitingPlantNumber)/sum(tmp$seedlingNumber)
# 
# # plot hierarchical population-level estimate
# 
# # toggle on jitter to avoid overplotting
# jitter <- rnorm(length(tmp$observed.p),mean=0,sd=0.01)
# 
# plot(tmp$observed.p + jitter,tmp$med,
#      pch=16, cex = 1,
#      col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#      xlim=c(0,1),ylim=c(0,1),
#      xlab="Observed proportion, y[i]/n[i]",
#      ylab="theta[i] and 95% credible interval",
#      main=yearNames[i],type='n')
# 
# segments(x0=tmp$observed.p + jitter,x1=tmp$observed.p + jitter,
#          y0=tmp$ci.lo,y1=tmp$ci.hi,lwd=0.75,
#          col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5))
# points(tmp$observed.p + jitter,tmp$med,
#        pch=16, cex = 1,
#        col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5))
# 
# # one-to-one line
# abline(a=0,b=1)
# 
# # add mle estimate for the year
# abline(h = mle, lty='dashed')
# }
# 
# # compare hierarchical to non-hierarchical
# # for posterior of p
# for(i in sampleYear){
#   
#   # filter data to year
#   tmp = simData %>% 
#     dplyr::filter(year==i) %>%
#     dplyr::mutate(observed.p = fruitingPlantNumber/seedlingNumber) %>%
#     dplyr::arrange(observed.p) %>%
#     tidyr::unite(yearPlot,c(year,plot))
#   
#   mcmcSummary<-samplesBinLinkBetaPriorPartialMean %>%
#     tidybayes::spread_draws(theta[year,plot]) %>%
#     tidyr::unite(yearPlot,c(year,plot)) %>%
#     dplyr::filter(yearPlot %in% tmp$yearPlot) %>%
#     dplyr::group_by(yearPlot) %>%
#     dplyr::summarise(med = median(theta), 
#                      ci.lo = quantile(theta,probs=0.025), 
#                      ci.hi = quantile(theta,probs=0.975))
#   
#   tmp<-tmp %>%
#     dplyr::left_join(mcmcSummary,by="yearPlot")
#   
#   mle = sum(tmp$fruitingPlantNumber)/sum(tmp$seedlingNumber)
#   
#   # plot hierarchical population-level estimate
#   
#   # toggle on jitter to avoid overplotting
#   jitter <- rnorm(length(tmp$observed.p),mean=0,sd=0.01)
#   
#   plot(tmp$observed.p + jitter,tmp$med,
#        pch=16, cex = 1,
#        col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5),
#        xlim=c(0,1),ylim=c(0,1),
#        xlab="Observed proportion, y[i]/n[i]",
#        ylab="theta[i] and 95% credible interval",
#        main=yearNames[i],type='n')
#   
#   segments(x0=tmp$observed.p + jitter,x1=tmp$observed.p + jitter,
#            y0=tmp$ci.lo,y1=tmp$ci.hi,lwd=0.75,
#            col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5))
#   points(tmp$observed.p + jitter,tmp$med,
#          pch=16, cex = 1,
#          col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5))
#   
#   # one-to-one line
#   abline(a=0,b=1)
#   
#   # add mle estimate for the year
#   abline(h = mle, lty='dashed')
#   
# }
# 
# # compare hierarchical to non-hierarchical
# # for posterior of p
# for(i in sampleYear){
#   
#   # filter data to year
#   tmp = simData %>% 
#     dplyr::filter(year==i) %>%
#     dplyr::mutate(observed.p = fruitingPlantNumber/seedlingNumber) %>%
#     dplyr::arrange(observed.p) %>%
#     tidyr::unite(yearPlot,c(year,plot))
#   
#   mcmcSummary<-samplesBinLinkBetaPriorPartialMeanHier %>%
#     tidybayes::spread_draws(theta[site,year,plot]) %>%
#     tidyr::unite(yearPlot,c(year,plot)) %>%
#     dplyr::filter(yearPlot %in% tmp$yearPlot) %>%
#     dplyr::group_by(yearPlot) %>%
#     dplyr::summarise(med = median(theta), 
#                      ci.lo = quantile(theta,probs=0.025), 
#                      ci.hi = quantile(theta,probs=0.975))
#   
#   tmp<-tmp %>%
#     dplyr::left_join(mcmcSummary,by="yearPlot")
#   
#   mle = sum(tmp$fruitingPlantNumber)/sum(tmp$seedlingNumber)
#   
#   # plot hierarchical population-level estimate
#   
#   # toggle on jitter to avoid overplotting
#   jitter <- rnorm(length(tmp$observed.p),mean=0,sd=0.01)
#   
#   plot(tmp$observed.p + jitter,tmp$med,
#        pch=16, cex = 1,
#        col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        xlim=c(0,1),ylim=c(0,1),
#        xlab="Observed proportion, y[i]/n[i]",
#        ylab="theta[i] and 95% credible interval",
#        main=yearNames[i],type='n')
#   
#   segments(x0=tmp$observed.p + jitter,x1=tmp$observed.p + jitter,
#            y0=tmp$ci.lo,y1=tmp$ci.hi,lwd=0.75,
#            col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
#   points(tmp$observed.p + jitter,tmp$med,
#          pch=16, cex = 1,
#          col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
#   
#   # one-to-one line
#   abline(a=0,b=1)
#   
#   # add mle estimate for the year
#   abline(h = mle, lty='dashed')
#   
# }
# 
# dev.off()
# 
# 
# pdf(file=paste0(dirFigures,"appendix-x-figure54.pdf"), width=8, height=6)
# # compare both hierarchical models
# # purple is year-level parameters
# # gray is year- and population-level parameters
# par(mfrow=c(2,5))
# # compare hierarchical to non-hierarchical
# # for posterior of p
# for(i in 1:nYears){
#   
#   # filter data to year
#   tmp = simData %>% 
#     dplyr::filter(year==i) %>%
#     dplyr::mutate(observed.p = fruitingPlantNumber/seedlingNumber) %>%
#     dplyr::arrange(observed.p) %>%
#     tidyr::unite(yearPlot,c(year,plot))
#   
#   mcmcSummary<-samplesBinLinkBetaPriorPartialMeanHier %>%
#     tidybayes::spread_draws(theta[site,year,plot]) %>%
#     tidyr::unite(yearPlot,c(year,plot)) %>%
#     dplyr::filter(yearPlot %in% tmp$yearPlot) %>%
#     dplyr::group_by(yearPlot) %>%
#     dplyr::summarise(med.h = median(theta), 
#                      ci.lo.h = quantile(theta,probs=0.025), 
#                      ci.hi.h = quantile(theta,probs=0.975))
#   
#   mcmcSummary2<-samplesBinLinkBetaPriorPartialMean %>%
#     tidybayes::spread_draws(theta[year,plot]) %>%
#     tidyr::unite(yearPlot,c(year,plot)) %>%
#     dplyr::filter(yearPlot %in% tmp$yearPlot) %>%
#     dplyr::group_by(yearPlot) %>%
#     dplyr::summarise(med = median(theta), 
#                      ci.lo = quantile(theta,probs=0.025), 
#                      ci.hi = quantile(theta,probs=0.975))
#   
#   tmp<-tmp %>%
#     dplyr::left_join(mcmcSummary,by="yearPlot") %>%
#     dplyr::left_join(mcmcSummary2,by="yearPlot")
#   
#   mle = sum(tmp$fruitingPlantNumber)/sum(tmp$seedlingNumber)
#   
#   # plot hierarchical population-level estimate
#   
#   # toggle on jitter to avoid overplotting
#   jitter <- rnorm(length(tmp$observed.p),mean=0,sd=0.01)
#   
#   plot(tmp$observed.p + jitter,tmp$med.h,
#        pch=16, cex = 1,
#        col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        xlim=c(0,1),ylim=c(0,1),
#        xlab="Observed proportion, y[i]/n[i]",
#        ylab="theta[i] and 95% credible interval",
#        main=yearNames[i],type='n')
#   
#   segments(x0=tmp$observed.p + jitter,x1=tmp$observed.p + jitter,
#            y0=tmp$ci.lo.h,y1=tmp$ci.hi.h,lwd=0.75,
#            col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
#   points(tmp$observed.p + jitter,tmp$med.h,
#          pch=16, cex = 1,
#          col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
#  
#   # toggle on jitter to avoid overplotting
#   jitter <- rnorm(length(tmp$observed.p),mean=0,sd=0.01)
#   
#   # lower hierarchy
#   segments(x0=tmp$observed.p + jitter,x1=tmp$observed.p + jitter,
#            y0=tmp$ci.lo,y1=tmp$ci.hi,lwd=0.75,
#            col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5))
#   points(tmp$observed.p + jitter,tmp$med,
#          pch=16, cex = 1,
#          col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5))
#   
#   # one-to-one line
#   abline(a=0,b=1)
#   
#   # add mle estimate for the year
#   abline(h = mle, lty='dashed')
#   
# }
# dev.off()
# 
# # compare pooling for upper level parameters
# #pdf(file=paste0(dirFigures,"appendix-x-figure54-2.pdf"), width=8, height=6)
# 
# par(mfrow=c(1,1))
#   
#   tmp = simData %>% 
#     dplyr::group_by(year) %>%
#     dplyr::summarise(y = sum(fruitingPlantNumber), n = sum(seedlingNumber)) %>%
#     dplyr::mutate(mle.year = y/n) %>%
#     dplyr::arrange(mle.year)
#   
#   mcmcSummary<-samplesBinLinkBetaPriorPartialMeanHier %>%
#     tidybayes::spread_draws(phi[site,year]) %>%
#     dplyr::group_by(year) %>%
#     dplyr::summarise(med.h = median(phi), 
#                      ci.lo.h = quantile(phi,probs=0.025), 
#                      ci.hi.h = quantile(phi,probs=0.975))
#   
#   mcmcSummary2<-samplesBinLinkBetaPriorPartialMean %>%
#     tidybayes::spread_draws(phi[year]) %>%
#     dplyr::group_by(year) %>%
#     dplyr::summarise(med = median(phi), 
#                      ci.lo = quantile(phi,probs=0.025), 
#                      ci.hi = quantile(phi,probs=0.975))
#   
#   tmp<-tmp %>%
#     dplyr::left_join(mcmcSummary,by="year") %>%
#     dplyr::left_join(mcmcSummary2,by="year")
#   
#   # filter data to year
#   mle.pop = sum(simData$fruitingPlantNumber)/sum(simData$seedlingNumber)
#   
#   # plot hierarchical population-level estimate
#   
# 
#   # toggle on jitter to avoid overplotting
#   jitter <- rnorm(length(tmp$mle.year),mean=0,sd=0.01)
#   
#   plot(tmp$mle.year + jitter,tmp$med,
#        pch=16, cex = 1,
#        col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5),
#        xlim=c(0,1),ylim=c(0,1),
#        xlab="Year-level MLE, sum(y[i])/sum(n[i])",
#        ylab="phi[i] and 95% credible interval",type='n')
# 
#   # toggle on jitter to avoid overplotting
#   jitter <- rnorm(length(tmp$mle.year),mean=0,sd=0.01)
#   
#   # 2-level hierarchical model
#   segments(x0=tmp$mle.year + jitter,x1=tmp$mle.year + jitter,
#            y0=tmp$ci.lo.h,y1=tmp$ci.hi.h,lwd=0.75,
#            col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
#   points(tmp$mle.year + jitter,tmp$med.h,
#          pch=16, cex = 1,
#          col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
#   
#   # 1-level hierarchical model
#   segments(x0=tmp$mle.year - jitter,x1=tmp$mle.year - jitter,
#            y0=tmp$ci.lo,y1=tmp$ci.hi,lwd=0.75,
#            col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5))
#   points(tmp$mle.year - jitter,tmp$med,
#          pch=16, cex = 1,
#          col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5))
#   
#   # one-to-one line
#   abline(a=0,b=1)
#   
#   # add mle estimate for the population
#   abline(h = mle.pop, lty='dashed')
#   
#   
#   mcmcSummaryPop<-samplesBinLinkBetaPriorPartialMeanHier %>%
#     tidybayes::spread_draws(phi0) %>%
#     dplyr::summarise(med = median(phi0), 
#                      ci.lo = quantile(phi0,probs=0.025), 
#                      ci.hi = quantile(phi0,probs=0.975))
#   
#   # add mle estimate for the population
#   abline(h = mcmcSummaryPop$med, lty=,
#          col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
#  
#   out=seq(-.2,1.2,.1)
#   yhi = rep(mcmcSummaryPop$ci.hi,length(out))
#   ylo = rep(mcmcSummaryPop$ci.lo,length(out))
# 
#   abline(h = mcmcSummaryPop$med, 
#          col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
# 
#   polygon(x=c(out,rev(out)),y=c(ylo,rev(yhi)),
#           col=rgb(red = 0, green = 0, blue = 0, alpha = 0.1),
#           border=NA)
# 
# 
#   # -------------------------------------------------------------------
#   # -------------------------------------------------------------------
# # Compare parameterizations
#   # -------------------------------------------------------------------
#   # -------------------------------------------------------------------
#   samplesBinLinkBetaPriorPartialMean %<>% tidybayes::recover_types(simData)
#   samplesBinLinkBetaPriorPartialMeanLogit %<>% tidybayes::recover_types(simData)
#   samplesBinLinkBetaPriorPartialMeanLogitNC %<>% tidybayes::recover_types(simData)
#   
#   ref<-simData %>% 
#     dplyr::select(year,plot) %>%
#     tidyr::unite(yearPlot,c(year,plot))
#   
#   
#   tmp = simData %>% 
#    # dplyr::filter(year==i) %>%
#      dplyr::mutate(observed.p = fruitingPlantNumber/seedlingNumber) %>%
#     dplyr::arrange(observed.p) %>%
#     tidyr::unite(yearPlot,c(year,plot))
#   
#   mcmcSummary<-samplesBinLinkBetaPriorPartialMean %>%
#     tidybayes::spread_draws(theta[year,plot]) %>%
#     tidyr::unite(yearPlot,c(year,plot)) %>%
#    # dplyr::filter(yearPlot %in% tmp$yearPlot) %>%
#     dplyr::group_by(yearPlot) %>%
#     dplyr::summarise(med = median(theta), 
#                      ci.lo = quantile(theta,probs=0.025), 
#                      ci.hi = quantile(theta,probs=0.975))
#   
#   mcmcSummary2<-samplesBinLinkBetaPriorPartialMeanLogit %>%
#     tidybayes::spread_draws(alpha[year,plot]) %>%
#     tidyr::unite(yearPlot,c(year,plot)) %>%
#     dplyr::mutate(theta=boot::inv.logit(alpha)) %>%
#     #dplyr::filter(yearPlot %in% tmp$yearPlot) %>% 
#     dplyr::group_by(yearPlot) %>%
#     dplyr::summarise(med.logit = median(theta), 
#                      ci.lo.logit = quantile(theta,probs=0.025), 
#                      ci.hi.logit = quantile(theta,probs=0.975))
#   
#   mcmcSummary3<-samplesBinLinkBetaPriorPartialMeanLogitNC %>%
#     tidybayes::spread_draws(alpha[year,plot]) %>%
#     tidyr::unite(yearPlot,c(year,plot)) %>%
#     dplyr::mutate(theta=boot::inv.logit(alpha)) %>%
#     #dplyr::filter(yearPlot %in% tmp$yearPlot) %>% 
#     dplyr::group_by(yearPlot) %>%
#     dplyr::summarise(med.logit.nc = median(theta), 
#                      ci.lo.logit.nc = quantile(theta,probs=0.025), 
#                      ci.hi.logit.nc = quantile(theta,probs=0.975))
#   
#   tmp<-tmp %>%
#     dplyr::left_join(mcmcSummary,by="yearPlot") %>%
#     dplyr::left_join(mcmcSummary2,by="yearPlot") %>%
#     dplyr::left_join(mcmcSummary3,by="yearPlot")
#   
#   mle = sum(tmp$fruitingPlantNumber)/sum(tmp$seedlingNumber)
#   
#   # plot hierarchical population-level estimate
#   
#   # toggle on jitter to avoid overplotting
#   jitter <- rnorm(length(tmp$observed.p),mean=0,sd=0.01)
#   
#   par(mfrow=c(1,2))
#   plot(tmp$observed.p + jitter,tmp$med,
#        pch=16, cex = 1,
#        col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5),
#        xlim=c(0,1),ylim=c(0,1),
#        xlab="Observed proportion, y[i]/n[i]",
#        ylab="theta[i] and 95% credible interval",
#        main="",type='n')
#   
#   segments(x0=tmp$observed.p + jitter,x1=tmp$observed.p + jitter,
#            y0=tmp$ci.lo,y1=tmp$ci.hi,lwd=0.75,
#            col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5))
#   points(tmp$observed.p + jitter,tmp$med,
#          pch=16, cex = 1,
#          col=rgb(red = .5, green = 0, blue = .5, alpha = 1))
#   
#   # one-to-one line
#   abline(a=0,b=1)
#   
#   # add mle estimate for the year
#   abline(h = mle, lty='dashed')
#   
#   plot(tmp$observed.p + jitter,tmp$med.logit,
#        pch=16, cex = 1,
#        col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5),
#        xlim=c(0,1),ylim=c(0,1),
#        xlab="Observed proportion, y[i]/n[i]",
#        ylab="theta[i] and 95% credible interval",
#        main="",type='n')
#   
#   segments(x0=tmp$observed.p + jitter,x1=tmp$observed.p + jitter,
#            y0=tmp$ci.lo.logit,y1=tmp$ci.hi.logit,lwd=0.75,
#            col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5))
#   points(tmp$observed.p + jitter,tmp$med.logit,
#          pch=16, cex = 1,
#          col=rgb(red = .5, green = 0, blue = .5, alpha = 1))
#   
#   # one-to-one line
#   abline(a=0,b=1)
#   
#   # add mle estimate for the year
#   abline(h = mle, lty='dashed')
# 
#   par(mfrow=c(1,3))
#   # compare  
#   plot(tmp$med ,tmp$med.logit,
#        pch=16,cex=.75,
#        xlim=c(0,1),ylim=c(0,1),
#        xlab="theta[i] and 95% credible interval\n direct parameterization",
#        ylab="theta[i] and 95% credible interval, logit-parameterization")
#   # one-to-one line
#   abline(a=0,b=1)
#   
#   #compare
#   plot(tmp$med ,tmp$med.logit.nc,
#        pch=16,cex=.75,
#        xlim=c(0,1),ylim=c(0,1),
#        xlab="theta[i] and 95% credible interval\n direct parameterization",
#        ylab="theta[i] and 95% credible interval, non-centered logit-parameterization")
#   # one-to-one line
#   abline(a=0,b=1)
#   
#   #compare
#   plot(tmp$med.logit ,tmp$med.logit.nc,
#        pch=16,cex=.75,
#        xlim=c(0,1),ylim=c(0,1),
#        xlab="theta[i] and 95% credible interval\n logit-parameterization",
#        ylab="theta[i] and 95% credible interval, non-centered logit-parameterization")
#   # one-to-one line
#   abline(a=0,b=1)
#   
#   
#   
# dim(MCMCchains(samplesBinLinkBetaPriorPartialMean,params="theta"))
# dim(MCMCchains(samplesBinLinkBetaPriorPartialMeanLogit,params='alpha'))
# theta <- MCMCchains(samplesBinLinkBetaPriorPartialMean,params="theta") %>%
#   recover_types(simData) %>%
#   spread_draws(theta[year,plot]) %>%
#   tidyr::unite(yearPlot, c(year,plot),remove=FALSE) %>%
#   dplyr::filter(yearPlot %in% ref$yearPlot) %>% 
#   unspread_draws(theta[year,plot]) %>% 
#   dplyr::select(-c('.chain','.iteration','.draw'))
#   
# theta.logit<- MCMCchains(samplesBinLinkBetaPriorPartialMeanLogit,params="alpha") %>%
#   recover_types(simData) %>%
#   spread_draws(alpha[year,plot]) %>%
#   dplyr::mutate(theta = boot::inv.logit(alpha)) %>%
#   dplyr::select(-alpha) %>%
#   tidyr::unite(yearPlot, c(year,plot),remove=FALSE) %>%
#   dplyr::filter(yearPlot %in% ref$yearPlot) %>% 
#   unspread_draws(theta[year,plot]) %>% 
#   dplyr::select(-c('.chain','.iteration','.draw'))
#   
# delta=theta.logit-theta
# delta.summary<-apply(delta,2,quantile,probs=c(.025,.5,.975))  
# 
# test<-delta.summary %>%
#   tidybayes::spread_draws(theta[year,plot]) %>%
#   dplyr::filter(.draw==2) %>%
#   dplyr::left_join(simData,by=c("year","plot"))
# 
# par(mfrow=c(1,1))
# plot(test$seedlingNumber,test$theta)
# 
# 
# par(mfrow=c(1,1))
# 
# tmp = simData %>% 
#   dplyr::group_by(year) %>%
#   dplyr::summarise(y = sum(fruitingPlantNumber), n = sum(seedlingNumber)) %>%
#   dplyr::mutate(mle.year = y/n) %>%
#   dplyr::arrange(mle.year)
# 
# mcmcSummary<-samplesBinLinkBetaPriorPartialMean %>%
#   tidybayes::spread_draws(phi[year]) %>%
#   dplyr::group_by(year) %>%
#   dplyr::summarise(med = median(phi), 
#                    ci.lo = quantile(phi,probs=0.025), 
#                    ci.hi = quantile(phi,probs=0.975))
# 
# mcmcSummary2<-samplesBinLinkBetaPriorPartialMeanLogit %>%
#   tidybayes::spread_draws(mu[year]) %>%
#   dplyr::group_by(year) %>%
#   dplyr::summarise(med.logit = median(boot::inv.logit(mu)), 
#                    ci.lo.logit = quantile(boot::inv.logit(mu),probs=0.025), 
#                    ci.hi.logit = quantile(boot::inv.logit(mu),probs=0.975))
# 
# tmp<-tmp %>%
#   dplyr::left_join(mcmcSummary,by="year") %>%
#   dplyr::left_join(mcmcSummary2,by="year")
# 
# # plot hierarchical population-level estimate
# 
# 
# # toggle on jitter to avoid overplotting
# jitter <- rnorm(length(tmp$mle.year),mean=0,sd=0.01)
# 
# plot(tmp$mle.year + jitter,tmp$med,
#      pch=16, cex = 1,
#      col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5),
#      xlim=c(0,1),ylim=c(0,1),
#      xlab="Year-level MLE, sum(y[i])/sum(n[i])",
#      ylab="phi[i] and 95% credible interval",type='n')
# 
# # toggle on jitter to avoid overplotting
# jitter <- rnorm(length(tmp$mle.year),mean=0,sd=0.01)
# 
# # 1-level hierarchical model, direct param
# segments(x0=tmp$mle.year + jitter,x1=tmp$mle.year + jitter,
#          y0=tmp$ci.lo,y1=tmp$ci.hi,lwd=0.75,
#          col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
# points(tmp$mle.year + jitter,tmp$med,
#        pch=16, cex = 1,
#        col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
# 
# # 1-level hierarchical model, logit param
# segments(x0=tmp$mle.year - jitter,x1=tmp$mle.year - jitter,
#          y0=tmp$ci.lo.logit,y1=tmp$ci.hi.logit,lwd=0.75,
#          col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5))
# points(tmp$mle.year - jitter,tmp$med.logit,
#        pch=16, cex = 1,
#        col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5))
# 
# # one-to-one line
# abline(a=0,b=1)
# 
# # add mle estimate for the population
# abline(h = mle.pop, lty='dashed')
# 
# 
# mcmcSummaryPop<-samplesBinLinkBetaPriorPartialMeanHier %>%
#   tidybayes::spread_draws(phi0) %>%
#   dplyr::summarise(med = median(phi0), 
#                    ci.lo = quantile(phi0,probs=0.025), 
#                    ci.hi = quantile(phi0,probs=0.975))
# 
# # add mle estimate for the population
# abline(h = mcmcSummaryPop$med, lty=,
#        col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
# 
# out=seq(-.2,1.2,.1)
# yhi = rep(mcmcSummaryPop$ci.hi,length(out))
# ylo = rep(mcmcSummaryPop$ci.lo,length(out))
# 
# abline(h = mcmcSummaryPop$med, 
#        col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
# 
# polygon(x=c(out,rev(out)),y=c(ylo,rev(yhi)),
#         col=rgb(red = 0, green = 0, blue = 0, alpha = 0.1),
#         border=NA)
# 
# 
# 
# tmp = simData %>% 
#   # dplyr::filter(year==i) %>%
#   dplyr::mutate(observed.p = fruitingPlantNumber/seedlingNumber) %>%
#   dplyr::arrange(observed.p) %>%
#   tidyr::unite(yearPlot,c(year,plot))
# 
#   samplesBinLinkBetaPriorComplete <- readRDS(simFiles[[5]])
#   
#   samplesBinLinkBetaPriorPartialMean <- readRDS(simFiles[[7]])
#   samplesBinLinkBetaPriorPartialMeanLogit <- readRDS(simFiles[[10]])
#   samplesBinLinkBetaPriorPartialMeanLogitNC <- readRDS(simFiles[[11]])
#   
#   samplesBinLinkBetaPriorPartialMeanHier <- readRDS(simFiles[[6]])
#   samplesBinLinkBetaPriorPartialMeanHierLogit <- readRDS(simFiles[[8]])
#   samplesBinLinkBetaPriorPartialMeanHierLogitNC <- readRDS(simFiles[[9]])
#   
# 
#   
#   # Comparison shows that the distribution of rhat for
#   # the logit parameterization is <1.05
#   # plot comparing Rhat values 
#   par(mfrow=c(1,3))  
#   hist(MCMCsummary(samplesBinLinkBetaPriorPartialMean,params=c("theta"))$Rhat,breaks=40)
#   hist(MCMCsummary(samplesBinLinkBetaPriorPartialMeanLogit,params=c("theta"))$Rhat,breaks=40)
#   hist(MCMCsummary(samplesBinLinkBetaPriorPartialMeanLogitNC,params=c("theta"))$Rhat,breaks=40)
# 
#   par(mfrow=c(1,3))  
#   MCMCsummary(samplesBinLinkBetaPriorPartialMean,params=c("phi"))
#   MCMCsummary(samplesBinLinkBetaPriorPartialMeanLogit,params=c("mu"))
#   MCMCsummary(samplesBinLinkBetaPriorPartialMeanLogitNC,params=c("mu"))
#   
#   par(mfrow=c(1,3))  
#   hist(MCMCsummary(samplesBinLinkBetaPriorPartialMeanHier,params=c("theta"))$Rhat,breaks=40)
#   hist(MCMCsummary(samplesBinLinkBetaPriorPartialMeanHierLogit,params=c("theta"))$Rhat,breaks=40)
#   hist(MCMCsummary(samplesBinLinkBetaPriorPartialMeanHierLogitNC,params=c("theta"))$Rhat,breaks=40)
#   
#   par(mfrow=c(1,3))  
#   MCMCsummary(samplesBinLinkBetaPriorPartialMeanHier,params=c("phi"))
#   MCMCsummary(samplesBinLinkBetaPriorPartialMeanHierLogit,params=c("mu"))
#   MCMCsummary(samplesBinLinkBetaPriorPartialMeanHierLogitNC,params=c("mu"))
#     
#   # look at posteriors
# phi<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params=c("phi"))
# kappa<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params=c("kappa"))
# theta<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params=c("theta"))
# 
# # see neal's funnel
# par(mfrow=c(2,5))
# for(i in 1:nYears){
# plot(boot::logit(phi[,i]),log(kappa[,i]))
# }
# 
# par(mfrow=c(1,1))
# # also the relationship of theta and kappa
# plot(boot::logit(theta[,1]),log(kappa[,1]))
# 
# 
# # look at posteriors
# mu<-MCMCchains(samplesBinLinkBetaPriorPartialMeanLogit,params=c("mu"))
# sigma<-MCMCchains(samplesBinLinkBetaPriorPartialMeanLogit,params=c("sigma"))
# alpha<-MCMCchains(samplesBinLinkBetaPriorPartialMeanLogit,params=c("alpha"))
# theta<-MCMCchains(samplesBinLinkBetaPriorPartialMeanLogit,params=c("theta"))
# 
# par(mfrow=c(1,1))
# # also the relationship of log-odds of success (alpha) and log population scale
# plot(alpha[,1],log(sigma[,1]))
# 
# 
# # look at posteriors non-centered parm
# mu<-MCMCchains(samplesBinLinkBetaPriorPartialMeanLogitNC,params=c("mu"))
# sigma<-MCMCchains(samplesBinLinkBetaPriorPartialMeanLogitNC,params=c("sigma"))
# alpha.std<-MCMCchains(samplesBinLinkBetaPriorPartialMeanLogitNC,params=c("alpha.std"))
# theta<-MCMCchains(samplesBinLinkBetaPriorPartialMeanLogitNC,params=c("theta"))
# 
# par(mfrow=c(1,1))
# plot(mu[,1],log(sigma[,1]))
# 
# # also the relationship of log-odds of success (alpha) and log population scale
# plot(alpha.std[,1],log(sigma[,1]))
# 
# # -------------------------------------------------------------------
# # -------------------------------------------------------------------
# # Posterior Predictive Checks
# 
# # notes: check group = interaction(dat$site,dat$yearStart)
# # for stat density grouped plots
# # -------------------------------------------------------------------
# # -------------------------------------------------------------------
# library(gridExtra)
# library(bayesplot)
# library(dplyr)
# 
# 
# iter<-dim(MCMCchains(mcmcSamples,params="fruitingPlantNumberSim"))[1]
# 
# #stat:skew
# skew <- psych::skew
# 
# # pdf(
# #   "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/products/seedBagsCompletePoolingPPC.pdf",
# #   onefile=TRUE,
# #   paper="USr",
# #   height = 5.5, width = 10)
# 
# # ppc_dens_overlay(simData$fruitingPlantNumber, MCMCchains(samplesBinLinkBetaPriorPartialMeanHier,params="fruitingPlantNumberSim")[sample(iter,1000), ]) +
# #   theme_bw() +
# #   labs(title="Posterior predictive checks for number of fruiting plants", 
# #   caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")
# 
# y <- simData$fruitingPlantNumber
# yrep.nh <- MCMCchains(samplesBinLinkBetaPriorComplete,params="fruitingPlantNumberSim")[sample(iter,1000), ]
# group <- as.factor(simData$year)
# 
# color_scheme_set("brightblue")
# 
# nh.mean<-ppc_stat_grouped(y, yrep.nh,
#                    group=group,facet_args=list(nrow=2)) +
#   theme_bw() +
#   labs(title="Posterior predictive checks for number of fruiting plants; non-hierarchical model", 
#        caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")
# 
# nh.var<-ppc_stat_grouped(y, yrep.nh, group=group,stat="var",facet_args=list(nrow=2)) +
#   theme_bw() +
#   labs(title="Posterior predictive checks for number of fruiting plants; non-hierarchical model", 
#        caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")
# 
# 
# yrep.hyear = MCMCchains(samplesBinLinkBetaPriorPartialMean,params="fruitingPlantNumberSim")[sample(iter,1000), ]
# 
# color_scheme_set("purple")
# 
# hieryear.mean<-ppc_stat_grouped(y, yrep.hyear,group=group,facet_args=list(nrow=2)) +
#   theme_bw() +
#   labs(title="Posterior predictive checks for number of fruiting plants; hierarchical model with year-level parameters", 
#        caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")
# 
# hieryear.var<-ppc_stat_grouped(y, yrep.hyear,group=group,stat='var',facet_args=list(nrow=2)) +
#   theme_bw()
#   labs(title="Posterior predictive checks for number of fruiting plants; hierarchical model with year-level parameters", 
#        caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")
# 
# 
#  yrep.hyearpop= MCMCchains(samplesBinLinkBetaPriorPartialMeanHier,params="fruitingPlantNumberSim")[sample(iter,1000), ]
# 
#  color_scheme_set("darkgray")
# 
#  
#  hieryearpop.mean<-ppc_stat_grouped(y, yrep.hyearpop,group=group,facet_args=list(nrow=2)) +
#   theme_bw() +
#   labs(title="Posterior predictive checks for number of fruiting plants; hierarchical model with year- and population-level parameters", 
#        caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")
# 
# hieryearpop.var<-ppc_stat_grouped(y, yrep.hyearpop,group=group,stat='var',facet_args=list(nrow=2)) +
#   theme_bw() +
#   labs(title="Posterior predictive checks for number of fruiting plants; hierarchical model with year- and population-level parameters", 
#        caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")
# 
# l = list(nh.mean,hieryear.mean,hieryearpop.mean,nh.var,hieryear.var,hieryearpop.var)
# 
# ggsave(paste0(dirFigures,"appendix-x-ppc.pdf"), marrangeGrob(grobs = l,nrow=1,ncol=1))
# 
# 
# 
# 
# par(mfrow=c(2,5))
# # compare hierarchical to non-hierarchical
# # for posterior of p
# for(i in 1:nYears){
#   # plot hierarchical population-level estimate
#   plot(density(hyearpop.p_site),
#        xlim=c(0,1),ylim=c(0,30),
#        xlab="Probability of survivorship",
#        main=yearNames[i])
#   
#   # plot hierarchical year-level estimates
#   lines(density(hyearpop.p[,i]),lty=1,lwd=1,col='red')
#   
#   lines(density(hyear.p[,i]),lty=1,lwd=1,col='blue')
#   
#   tmp = simData %>% dplyr::filter(year==i)
#   points(tmp$fruitingPlantNumber/tmp$seedlingNumber,
#          rep(0,dim(tmp)[1]),
#          pch=16,
#          col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
# }
# 
# hyear.p.sum<-apply(hyear.p,2,quantile,probs=c(.025,.5,.975))
# hyearpop.p.sum<-apply(hyearpop.p,2,quantile,probs=c(.025,.5,.975))
# 
# # show summary of posteriors with 95% credible intervals
# par(mfrow=c(1,1))
# plot(1:9-0.25,t(hyear.p.sum)[,2],pch=1,ylim=c(0,1),xlim=c(0,10))
# segments(x0=1:9-0.25,y0=t(hyear.p.sum)[,1],
#          x1=1:9-0.25,y1=t(hyear.p.sum)[,3])
# points(1:9+0.25,t(hyearpop.p.sum)[,2],pch=3)
# segments(x0=1:9+0.25,y0=t(hyearpop.p.sum)[,1],
#          x1=1:9+0.25,y1=t(hyearpop.p.sum)[,3])
# 
# yearNames = levels(simData$year)
# 
# # par(mfrow=c(2,5))
# # for(i in 1:nYears){
# #   plot(density(p[,i]),xlim=c(-.1,1.1),main=yearNames[i])
# #   #lines(density(omega[,i]),lty='dotted')
# #   lines(density(phi[,i]),lty='dashed')
# # }
# # 
# # simData$year <- as.numeric(simData$year)
# # ggplot() +
# #   stat_density(data=samplesBinLinkBetaPriorComplete %>%
# #                  tidybayes::spread_draws(p[year]),aes(p),
# #                geom='line') +
# #   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
# #                  tidybayes::spread_draws(phi[year]),aes(phi),color="red",
# #                geom='line') +
# #   facet_wrap(~year,scales="free",nrow=2) +
# #   geom_jitter(data=simData ,aes(x=fruitingPlantNumber/seedlingNumber,y=0)) +
# #   theme_minimal()
# 
# # posteriors for the complete pooling model, the mean of the beta, the mode of the beta
# # it shows that the complete pooling model (one p per year, black)
# # has the smallest credible intervals which generally
# # overlap those from the model parameterized via the mean (red)
# # the models parameterized via the mean are slightly broader credible intervals
# 
# # ggplot() +
# #   stat_density(data=samplesBinLinkBetaPriorComplete %>%
# #                  tidybayes::spread_draws(p[year]),aes(p),geom='line') +
# #   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
# #                  tidybayes::spread_draws(phi[year]),aes(phi),color="red",geom='line') +
# #   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
# #                  tidybayes::spread_draws(p[year]),aes(p),color="red",geom='line',linetype='dotted') +
# #   stat_density(data=samplesBinLinkBetaPriorPartialMode %>%
# #                  tidybayes::spread_draws(omega[year]),aes(omega),color="blue",geom='line') +
# #   stat_density(data=samplesBinLinkBetaPriorPartialMode %>%
# #                  tidybayes::spread_draws(p[year]),aes(p),color="blue",geom='line',linetype='dotted') +
# #   facet_wrap(~year,scales="free",nrow=2) +
# #   geom_jitter(data=simData ,aes(x=fruitingPlantNumber/seedlingNumber,y=0)) +
# #   theme_minimal()
# 
# # the posterior distribution for probability of success in year p (marginalized over all plots)
# # is broader than the posterior for the mean probability of success in year p
# # parameterization via the mean and mode give similar results once marginalized
# # though the mode parameterization seems to give broader posteriors
# 
# simData$year <- as.numeric(as.factor(simData$year))
# simData$plot <- as.numeric(as.factor(simData$plot))
# # ggplot() +
# #   stat_density(data=samplesBinLinkBetaPriorComplete %>%
# #                  tidybayes::spread_draws(p[year]) %>%
# #                  dplyr::filter(year==5),aes(p),lwd=1, geom='line') +
# #   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
# #                  tidybayes::spread_draws(phi[year]) %>% 
# #                  dplyr::filter(year==5),aes(phi),color="red",lwd=1, geom='line') +
# #   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
# #                  tidybayes::spread_draws(p[year]) %>% 
# #                  dplyr::filter(year==5),aes(p),
# #                color="red",linetype='dashed',lwd=0.5, geom="line") +
# #   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
# #                  tidybayes::spread_draws(theta[year,plot]) %>% 
# #                  dplyr::filter(year==5) ,aes(theta,group=plot),
# #                color="black",lwd=0.25, geom="line") +
# #   geom_point(data=simData  %>% 
# #                dplyr::filter(year==5),
# #              aes(x=fruitingPlantNumber/seedlingNumber,y=0),position=ggstance::position_dodgev(height=0.3)) +
# #   theme_minimal() +
# #   facet_wrap(~plot, scales = 'free')
# 
# 
# samplesBinLinkBetaPriorPartialMeanHier %<>% recover_types(simData)
# # plot 2
# ggplot() +
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(phi0),aes(phi0),color="red",geom='line',size=1) +
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(phi[site,year]),aes(phi,group=year),geom='line',size=0.25, position='identity') +
#   geom_point(data=simData,
#              aes(x=fruitingPlantNumber/seedlingNumber,y=0),
#              position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
#   # facet_wrap(~year,scales="free_y",nrow=2) +
#   theme_minimal()
# 
# ggplot() +
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(p_site),aes(p_site),color="red",geom='line',size=1) +
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(p[site,year]),aes(p,group=year),geom='line',size=0.25, position='identity') +
#   geom_point(data=simData,
#              aes(x=fruitingPlantNumber/seedlingNumber,y=0),
#              position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
#   facet_wrap(~year,scales="free_y",nrow=2) +
#   theme_minimal()
# 
# # compare levels of the hierarchy
# ggplot() +
#   # posterior distribution of the mean probability of success at the site
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(phi0),aes(phi0),color="red",geom='line',size=.5) +
#   # posterior distribution of the mean probability of success in each year
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(phi[site,year]),aes(phi,group=year),geom='line',size=0.25, position='identity') +
#   # posterior distribution for the probability of success at the site
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(p_site),aes(p_site),color="red",geom='line',size=1) +
#   # posterior distribution for the probability of success in each year 
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(p[site,year]),aes(p,group=year),geom='line',size=1, position='identity') +
#   geom_point(data=simData,
#              aes(x=fruitingPlantNumber/seedlingNumber,y=0),
#              position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
#   facet_wrap(~year,scales="free_y",nrow=2) +
#   theme_minimal()
# 
# # key takeaway from the plot above: 
# # the posterior for probability of success at each level (site, year-within-site) is broader than the mean probability of success at that level (compare thin to solid lines)
# 
# ggplot() +
#   # posterior distribution of the mean probability of success at the site
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(phi0),aes(phi0),color="red",geom='line',size=1) +
#   # posterior distribution of the mean probability of success in each year
#   # stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#   #                tidybayes::spread_draws(phi[site,year]),aes(phi,group=as.factor(year)),geom='line',size=0.25, position='identity') +
#   # posterior distribution for the probability of success at the site
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(p_site),aes(p_site),color="red",geom='line',size=.5) +
#   # posterior distribution for the probability of success in each year 
#   # stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#   #                tidybayes::spread_draws(p[site,year]),aes(p,group=year),geom='line',size=1, position='identity') +
#   geom_point(data=simData,
#              aes(x=fruitingPlantNumber/seedlingNumber,y=0),
#              position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
#   #facet_wrap(~year,scales="free_y",nrow=2) +
#   theme_minimal()
# 
# ggplot() +
#   # posterior distribution of the mean probability of success at the site
#   # stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#   #                tidybayes::spread_draws(phi0),aes(phi0),color="red",geom='line',size=1) +
#   # posterior distribution of the mean probability of success in each year
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(phi[site,year]),aes(phi,group=as.factor(year)),geom='line',size=0.25, position='identity') +
#   # posterior distribution for the probability of success at the site
#   # stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#   #                tidybayes::spread_draws(p_site),aes(p_site),color="red",geom='line',size=.5) +
#   # posterior distribution for the probability of success in each year 
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(p[site,year]),aes(p,group=year),geom='line',size=1, position='identity') +
#   geom_point(data=simData,
#              aes(x=fruitingPlantNumber/seedlingNumber,y=0),
#              position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
#   facet_wrap(~year,scales="free_y",nrow=2) +
#   theme_minimal()
# 
# # effect of hierarchy on estimate of means - estimates pool towards hierarchy with little data
# ggplot() +
#   # posterior distribution of the mean probability of success at the site
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(phi0),aes(phi0),color="red",geom='line',size=1) +
#   # posterior distribution of the mean probability of success in each year
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(phi[site,year]),aes(phi,group=as.factor(year)),geom='line',size=0.25, position='identity',color='red') +
#   # posterior distribution of the mean probability of success in each year
#   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
#                  tidybayes::spread_draws(phi[year]),aes(phi,group=as.factor(year)),geom='line',size=0.25, position='identity') +
#   geom_point(data=simData,
#              aes(x=fruitingPlantNumber/seedlingNumber,y=0),
#              position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
#   facet_wrap(~year,scales="free_y",nrow=2) +
#   theme_minimal()
# 
# 
# 
# 
# 
# # compare the three different fits (no hierarchy, hierarchy at 1 level, hierarchy at 2 levels)
# # difference is greatest with small datasets (year = 9 only has 1 data point)
# ggplot() +
#   stat_density(data=samplesBinLinkBetaPriorComplete %>%
#                  tidybayes::spread_draws(p[year]),aes(p),color="black",geom='line',size=1) +
#   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
#                  tidybayes::spread_draws(phi[year]),aes(phi),color="red",geom='line') +
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(phi[site,year]),aes(phi),color="blue",geom='line') +
#   geom_point(data=simData,
#              aes(x=fruitingPlantNumber/seedlingNumber,y=0),
#              position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
#   facet_wrap(~year,scales="free_y",nrow=2) +
#   theme_minimal()
# 
# # compare distribution of p
# ggplot() +
#   stat_density(data=samplesBinLinkBetaPriorComplete %>%
#                  tidybayes::spread_draws(p[year]),aes(p),color="black",geom='line',size=1) +
#   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
#                  tidybayes::spread_draws(p[year]),aes(p),color="red",geom='line') +
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(p[site,year]),aes(p),color="blue",geom='line') +
#   geom_point(data=simData,
#              aes(x=fruitingPlantNumber/seedlingNumber,y=0),
#              position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
#   facet_wrap(~year,scales="free_y",nrow=2) +
#   theme_minimal()
# 
# 
# nohier <-apply(MCMCchains(samplesBinLinkBetaPriorComplete,params='p'),2,median)
# hier1 <- apply(MCMCchains(samplesBinLinkBetaPriorPartialMean,params='phi'),2,median)
# hier2 <-apply(MCMCchains(samplesBinLinkBetaPriorPartialMeanHier,params='phi'),2,median)
# 
# plot(nohier,hier1);abline(a=0,b=1)
# plot(nohier,hier2);abline(a=0,b=1)
# plot(hier1,hier2);abline(a=0,b=1)
# 
# samplesBinLinkBetaPriorPartialMean %>%
#   tidybayes::spread_draws(theta[year,plot])
# 
# ref<-simData %>% 
#   dplyr::select(year,plot) %>%
#   tidyr::unite(yearPlot,c(year,plot))
# 
# theta.med<-samplesBinLinkBetaPriorPartialMean %>%
#   tidybayes::spread_draws(theta[year,plot]) %>%
#   tidyr::unite(yearPlot,c(year,plot)) %>%
#   dplyr::filter(yearPlot %in% ref$yearPlot) %>%
#   tidyr::separate(yearPlot, c("year","plot")) %>%
#   dplyr::group_by(year,plot) %>%
#   dplyr::summarise(theta.med = median(theta))
# 
# theta.med$year <- as.numeric(theta.med$year)
# theta.med$plot <- as.numeric(theta.med$plot)
# 
# out<-simData %>%
#   dplyr::left_join(theta.med,by=c("year","plot"))
# 
# par(mfrow=c(1,1))
# plot(out$fruitingPlantNumber/out$seedlingNumber,
#      out$theta.med)
# 
# plot(out$seedlingNumber,out$theta.med-out$fruitingPlantNumber/out$seedlingNumber)
# 
# p.summary<-samplesBinLinkBetaPriorPartialMean %>%
#   tidybayes::spread_draws(p[year]) %>%
#   dplyr::summarise(p.med = median(p))
# p.summary$year <- as.numeric(p.summary$year)
# 
# phi.summary<-samplesBinLinkBetaPriorPartialMean %>%
#   tidybayes::spread_draws(phi[year]) %>%
#   dplyr::summarise(phi.med = median(phi))
# phi.summary$year <- as.numeric(phi.summary$year)
# 
# ggplot() +
#   geom_point(data=out,aes(x=fruitingPlantNumber/seedlingNumber,y=theta.med)) +
#   geom_abline(intercept=0,slope=1,lwd=0.5) +
#   geom_hline(data=p.summary,aes(yintercept=p.med),size=0.5) +
#   geom_hline(data=phi.summary,aes(yintercept=phi.med),size=0.5,color='red') +
#   
#   facet_wrap(~year) +
#   theme_minimal()
# 
# 
# theta.med<-samplesBinLinkBetaPriorPartialMean %>%
#   tidybayes::spread_draws(p[year]) %>%
#   dplyr::group_by(year) %>%
#   dplyr::summarise(p.med = median(p))
# 
# theta.med$year <- as.numeric(theta.med$year)
# theta.med$plot <- as.numeric(theta.med$plot)
# 
# out<-simData %>%
#   dplyr::left_join(theta.med,by=c("year","plot"))
# 
# par(mfrow=c(1,1))
# plot(out$fruitingPlantNumber/out$seedlingNumber,
#      out$theta.med)
# 
# 
# thetaBinLinkBetaPriorPartialMode <- MCMCchains(samplesBinLinkBetaPriorPartialMode,params="theta")
# 
# # plot(density(thetaBinLinkBetaPriorPartialMode[,1]), type='n', 
# #      xlim=c(-.2,1.2),ylim=c(0,30))
# # for(i in 211:241){
# #   lines(density(thetaBinLinkBetaPriorPartialMode[,i]))
# # }
# 
# omega<-MCMCchains(samplesBinLinkBetaPriorPartialMode,params="omega")
# 
# plot(density(omega[,1]), type='n',
#      xlim=c(-.2,1.2),ylim=c(0,80))
# for(i in 1:10){
#   lines(density(omega[,i]))
# }
# 
# kappa<-MCMCchains(samplesBinLinkBetaPriorPartialMode,params="kappa")
# 
# alpha = omega*(kappa-2)+1
# beta = (1-omega)*(kappa-2)+1
# plot(seq(0,1,by=0.001),
#      dbeta(seq(0,1,by=0.001),median(alpha),median(beta)),
#      type='n',
#      ylim=c(0,20))
# for(i in 1:20){
#   lines(seq(0,1,by=0.001),
#         dbeta(seq(0,1,by=0.001),median(alpha[i]),median(beta[i])),
#         type='l',lwd=.5,col='red')
# }
# 
# 
# thetaBinLinkBetaPriorPartialMean <- MCMCchains(samplesBinLinkBetaPriorPartialMean,params="theta")
# 
# plot(density(thetaBinLinkBetaPriorPartialMean[,1]), type='n', 
#      xlim=c(-.2,1.2),ylim=c(0,30))
# for(i in 211:241){
#   lines(density(thetaBinLinkBetaPriorPartialMean[,i]))
# }
# 
# phi<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params="phi")
# 
# 
# plot(density(phi[,1]), type='n',
#      xlim=c(-.2,1.2),ylim=c(0,80))
# for(i in 1:10){
#   lines(density(phi[,i]))
# }
# 
# 
# kappa<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params="kappa")
# alpha = kappa*phi
# beta = kappa*(1-phi)
# plot(seq(0,1,by=0.001),dbeta(seq(0,1,by=0.001),median(alpha),median(beta)),
#      type='n')
# for(i in 1:20){
#   lines(seq(0,1,by=0.001),
#         dbeta(seq(0,1,by=0.001),median(alpha[i]),median(beta[i])),
#         type='l',lwd=.5,col='red')
# }
# 
# p<-MCMCchains(samplesBinLinkBetaPriorComplete,params='p')
# 
# # compare posteriors for the hyperparameters
# # the parameterization via the mode seems to sample better in general
# phi<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params='phi')
# kappa.mean<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params='kappa')
# thin <- function(x) sample(x,1000)
# par(mfrow=c(2,5))
# for(i in 1:10){
#   plot(thin(phi[,i]),thin(kappa.mean[,i]),cex=.5,pch=16)
# }
# 
# 
# omega<-MCMCchains(samplesBinLinkBetaPriorPartialMode,params='omega')
# kappa.mode<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params='kappa')
# thin <- function(x) sample(x,1000)
# par(mfrow=c(2,5))
# for(i in 1:10){
#   plot(thin(omega[,i]),thin(kappa.mode[,i]),cex=.5,pch=16)
# }
# 
# plot(apply(p,2,median),apply(phi,2,median),xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)
# plot(apply(p,2,median),apply(omega,2,median),xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)
# plot(apply(phi,2,median),apply(omega,2,median),xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)
# 
# stat.summary <- function(x){
#   
#   x %<>% recover_types(simData)
#   
#   tmp <- x %>% 
#     tidybayes::spread_draws(p[year]) %>%
#     dplyr::summarise(ci.lo = quantile(theta, prob = c(.025)),
#                      med = quantile(theta, prob = c(.5)),
#                      ci.hi = quantile(theta, prob = c(.975)))
#   return(tmp)
# }
# 
# thetaBinLinkBetaPriorPartialMode <- stat.summary(samplesBinLinkBetaPriorPartialMode)
# thetaBinLinkBetaPriorPartialMean <- stat.summary(samplesBinLinkBetaPriorPartialMean)
# 
# d<-data.frame(thetaBinLinkBetaPriorPartialMode,thetaBinLinkBetaPriorPartialMean) %>% 
#   dplyr::mutate(abs(med.1-med)) %>%
#   dplyr::left_join(simData,by=c('year','plot'))
# 
# d<-d %>% dplyr::filter(!is.na(fruitingPlantNumber))
# 
# plot(d$med,d$med.1, 
#      pch=16,cex = 0.5,xlim=c(0,1),ylim=c(0,1))
# # segments(x0=thetaBinLinkBetaPriorPartialMode$med,y0=thetaBinLinkBetaPriorPartialMean$ci.lo,
# #          x1=thetaBinLinkBetaPriorPartialMode$med,y1=thetaBinLinkBetaPriorPartialMean$ci.hi)
# # segments(x0=thetaBinLinkBetaPriorPartialMode$ci.lo,y0=thetaBinLinkBetaPriorPartialMean$med,
# #          x1=thetaBinLinkBetaPriorPartialMode$ci.hi,y1=thetaBinLinkBetaPriorPartialMean$med)
# abline(a=0,b=1)
# 
# mle<-simData %>%
#   dplyr::group_by(year,plot) %>%
#   dplyr::summarise(y = sum(fruitingPlantNumber), n = sum(seedlingNumber)) %>%
#   dplyr::mutate(pHat = y/n) %>%
#   dplyr::left_join(thetaBinLinkBetaPriorPartialMode)
# 
# points(x=mle$med,y=mle$pHat,pch=2)
# 
# ## summary
# MCMCsummary(samplesBinLinkBetaPriorPartialMean,params=c("phi","kappa"))
# MCMCsummary(samplesBinLinkBetaPriorPartialMode,params=c("omega","kappa"))
# 
# samplesBinLinkBetaPriorPartialMean %<>% recover_types(simData)
# 
# d1 <- samplesBinLinkBetaPriorPartialMean %>%
#   tidybayes::spread_draws(theta[year,plot]) %>%
#   dplyr::filter(year==2014)
# 
# d1.summary <- d1 %>% dplyr::group_by(plot) %>%
#   dplyr::summarise(mean.med=median(theta))
# 
# samplesBinLinkBetaPriorPartialMode %<>% recover_types(simData)
# 
# d2 <- samplesBinLinkBetaPriorPartialMode %>%
#   tidybayes::spread_draws(theta[year,plot]) %>%
#   dplyr::filter(year==2014)
# 
# d2.summary <- d2 %>% dplyr::group_by(plot) %>%
#   dplyr::summarise(mode.med=median(theta))
# 
# d <- d1.summary %>%
#   dplyr::left_join(d2.summary, by= "plot")
# plot(d$mode.med,d$mean.med,xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)
# 
# simData %>%
#   dplyr::filter(year==2014) %>%
#   dplyr::left_join(d,by="plot") %>% View
# pdf(file=paste0(dirFigures,"appendix-x-mle_bayes.pdf"), width=8, height=4)
# 
