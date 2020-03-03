#################################################################################
################################################################################
################################################################################
# Code for figures to compare the following modeling approaches for the seedling survivorship data
# 1. binomial likelihood
# 2. binomial likelihood with beta prior, complete pooling
# 3. binomial likelihood with beta prior, partial pooling, parameterize mode
# 4. binomial likelihood with beta prior, partial pooling, parameterize mean
# 4. binomial likelihood with logit parameterization
#
# Scripts by Gregor Siegmund
# fausto.siegmund@gmail.com
# last updated 02-29-2020
#################################################################################
#################################################################################
#################################################################################
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)


directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/"
survivorshipFitFiles <- paste0(directory,list.files(directory))[-c(1:4)]
dirFigures = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/")

samplesBinLinkBetaPriorComplete <- readRDS(survivorshipFitFiles[[1]])
samplesBinLinkBetaPriorPartialMean <- readRDS(survivorshipFitFiles[[2]])
samplesBinLinkBetaPriorPartialMode <- readRDS(survivorshipFitFiles[[3]])
samplesBinLinkLogitLinkPartial<- readRDS(survivorshipFitFiles[[4]])
survDataAnalysisBayes <- readRDS(survivorshipFitFiles[[5]])
mleDataLong <- readRDS(survivorshipFitFiles[[6]])

nYears = 10
nSites = 20

################################################################################
# Load packages
################################################################################
library(janitor)
library(tidyverse)
library(tidybayes)
library(rjags)
library(magrittr) # for %<>% with recover_types
library(dplyr)
library(MCMCvis)

################################################################################
# Maximum likelihood estimate of survivorship, binomial likelihood
# by minimizing the negative log likelihood
################################################################################

# Figure
# Comparing the geographic pattern of MLE (red) to geographic pattern of all data
# ggplot() +
#   geom_point(data=survDataAnalysis%>% dplyr::left_join(siteEastingData,by="site"),aes(x=easting,y=fruitingPlantNumber/seedlingNumber),cex=0.25) +
#   geom_point(data=data.frame(mleDataLong)%>% dplyr::left_join(siteEastingData,by="site"),aes(x=easting,y=pHat),color="red") +
#   geom_smooth(data=survDataAnalysis%>% dplyr::left_join(siteEastingData,by="site"),aes(x=easting,y=fruitingPlantNumber/seedlingNumber),method="lm") +
#   geom_smooth(data=data.frame(mleDataLong)%>% dplyr::left_join(siteEastingData,by="site"),aes(x=easting,y=pHat),color="red",method="lm") +
#     facet_wrap(~year) +
#   theme_minimal()

################################################################################
# JAGS fit for model with binomial likelihood, beta prior, complete pooling
################################################################################

summaryBinLinkBetaPriorComplete<-samplesBinLinkBetaPriorComplete %>% 
  tidybayes::spread_draws(p[site,year]) %>%
  dplyr::summarise(bayesBBcomp.q50=median(p)) 

summaryBinLinkBetaPriorComplete$year <- as.character(summaryBinLinkBetaPriorComplete$year)
summaryBinLinkBetaPriorComplete$site <- as.factor(summaryBinLinkBetaPriorComplete$site)

### Figures for comparison

samplesizeData <- survDataAnalysisBayes %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n())

samplesizeData$site <- as.factor(samplesizeData$site)
samplesizeData$year <- as.numeric(as.character(samplesizeData$year))
summaryBinLinkBetaPriorComplete$year <- as.numeric(as.character(summaryBinLinkBetaPriorComplete$year))

estimatesComparisonData <- mleDataLong %>%
  dplyr::left_join(summaryBinLinkBetaPriorComplete,by=c("site","year")) %>%
  dplyr::left_join(samplesizeData,by=c("site","year"))

estimateComparisonDataArranged <- estimatesComparisonData %>%
  dplyr::arrange(year)

# calculate difference of posterior - max. likelihood estimate
post.chain.pars <-MCMCvis::MCMCchains(samplesBinLinkBetaPriorComplete,params="p")
estimateComparisonDataArranged<-estimatesComparisonData %>%
  dplyr::arrange(year)
tmp<-matrix(NA,ncol=200,nrow=length(samplesBinLinkBetaPriorComplete)*dim(samplesBinLinkBetaPriorComplete[[1]]))
for(i in 1:200){
  tmp[,i]<-post.chain.pars[,i]-estimateComparisonDataArranged$pHat[i]
}
delta.conf <- apply(tmp,2,quantile,probs=c(.025,.5,.975))

# plot comparing MLE fit and bayes fit 

pdf(file=paste0(dirFigures,"appendix-x-mle_bayes.pdf"), width=8, height=4)

par(mfrow=c(1,2))
par(mar=c(5, 6, 4, 2) + 0.1)
plot(x = estimatesComparisonData$pHat,
     y = estimatesComparisonData$bayesBBcomp.q50,
     xlab="Maximum likelihood estimate",
     ylab="Median of posterior for \n beta-binomial parameterization",
     xlim=c(0,1),ylim=c(0,1),
     pch=16,cex=.5)  
abline(a=0,b=1,lwd=0.5)

par(mar=c(5, 6, 4, 2) + 0.1)
plot(estimateComparisonDataArranged$`n()`,t(delta.conf)[,2],
     pch=16,cex=0.5,
     ylim=c(-1,1),
     xlab="Sample size (n)",
     ylab="Comparison of the maximum likelihood estimate\n and the beta-binomial parameterization\n (median and 95% CI)")
segments(x0=estimateComparisonDataArranged$`n()`,y0=t(delta.conf)[,1],
         x1=estimateComparisonDataArranged$`n()`,y1=t(delta.conf)[,3])
abline(h=0,lwd=1)

dev.off()

################################################################################
# JAGS fit for model with binomial likelihood, beta prior, partial pooling
################################################################################

summaryBinLinkBetaPriorPartialMode<-samplesBinLinkBetaPriorPartialMode %>% 
  tidybayes::spread_draws(theta[site,year]) %>%
  dplyr::summarise(bayesBBpartial.q50=median(theta)) 

summaryBinLinkBetaPriorPartialMode$year <- as.character(summaryBinLinkBetaPriorPartialMode$year)
summaryBinLinkBetaPriorPartialMode$site <- as.factor(summaryBinLinkBetaPriorPartialMode$site)
summaryBinLinkBetaPriorPartialMode$year <- as.numeric(as.character(summaryBinLinkBetaPriorPartialMode$year))

estimatesComparisonData<-  estimatesComparisonData %>%
  dplyr::left_join(summaryBinLinkBetaPriorPartialMode,by=c("site","year"))

# calculate difference of posterior - max. likelihood estimate
post.chain.pars <-MCMCvis::MCMCchains(samplesBinLinkBetaPriorPartialMode,params="theta")
estimateComparisonDataArranged<-estimatesComparisonData %>%
  dplyr::arrange(year)
tmp<-matrix(NA,ncol=200,nrow=length(samplesBinLinkBetaPriorPartialMode)*dim(samplesBinLinkBetaPriorPartialMode[[1]]))
for(i in 1:200){
  tmp[,i]<-post.chain.pars[,i]-estimateComparisonDataArranged$pHat[i]
}
delta.conf <- apply(tmp,2,quantile,probs=c(.025,.5,.975))

# plot comparing MLE fit and bayes fit 

pdf(paste0(dirFigures,"appendix-x-mle_bayeshier.pdf"), width=8, height=4)
par(mfrow=c(1,2))
par(mar=c(5, 6, 4, 2) + 0.1)
plot(estimatesComparisonData$pHat,estimatesComparisonData$bayesBBpartial.q50,
     xlab="Maximum likelihood estimate",
     ylab="Median of posterior for \n beta-binomial parameterization",
     xlim=c(0,1),ylim=c(0,1),
     pch=16,cex=.5)  
abline(a=0,b=1,lwd=0.5)

par(mar=c(5, 6, 4, 2) + 0.1)
plot(estimateComparisonDataArranged$`n()`,t(delta.conf)[,2],
     pch=16,cex=0.5,
     ylim=c(-1,1),
     xlab="Sample size (n)",
     ylab="Comparison of the maximum likelihood estimate\n and the beta-binomial parameterization\n (median and 95% CI)")
segments(x0=estimateComparisonDataArranged$`n()`,y0=t(delta.conf)[,1],
         x1=estimateComparisonDataArranged$`n()`,y1=t(delta.conf)[,3])
abline(h=0,lwd=1)

dev.off()

sum <- summaryBinLinkBetaPriorComplete %>%
  dplyr::left_join(summaryBinLinkBetaPriorPartialMode,by=c("site","year"))

plot(sum$bayesBBcomp.q50,sum$bayesBBpartial.q50)
abline(a=0,b=1)

tmp<-matrix(NA,ncol=200,nrow=length(samplesBinLinkBetaPriorPartialMode)*dim(samplesBinLinkBetaPriorPartialMode[[1]]))
for(i in 1:200){
  tmp[,i]<-MCMCvis::MCMCchains(samplesBinLinkBetaPriorComplete,params="p")[,i] - MCMCvis::MCMCchains(samplesBinLinkBetaPriorPartialMode)[,i]
}
delta.conf <- apply(tmp,2,quantile,probs=c(.025,.5,.975))

pdf(paste0(dirFigures,"appendix-x-bayescomp_bayeshier.pdf"), width=4, height=4)
par(mar=c(5, 6, 4, 2) + 0.1)
plot(estimateComparisonDataArranged$`n()`,t(delta.conf)[,2],
     pch=16,cex=0.5,
     ylim=c(-1,1),
     xlab="Sample size (n)",
     ylab="Comparison of the posteriors for \n complete and partial pooling parameterizations \n complete-partial (median and 95% CI)")
segments(x0=estimateComparisonDataArranged$`n()`,y0=t(delta.conf)[,1],
         x1=estimateComparisonDataArranged$`n()`,y1=t(delta.conf)[,3])
abline(h=0,lwd=1)
dev.off()


pdf(paste0(dirFigures,"appendix-x-mismatch.pdf"), width=8, height=6)
v<-(1:200)[estimateComparisonDataArranged$`n()`==1][!is.na((1:200)[estimateComparisonDataArranged$`n()`==1])]

omegaSite<-c(11,11,13,6,7,11,12,20)

par(mfrow=c(2,4))
for(i in 1:length(v)){
  plot(density(  MCMCvis::MCMCchains(samplesBinLinkBetaPriorComplete,params="p")[,v[i]] ),
       xlim=c(0,1),ylim=c(0,4),
       main="") 
  abline(v=median( MCMCvis::MCMCchains(samplesBinLinkBetaPriorComplete,params="p")[,v[i]] ))
  tmp <- sample(MCMCvis::MCMCchains(samplesBinLinkBetaPriorPartialMode,params="theta")[,v[i]],3000)
  lines( density(tmp),lty='dotted')
  abline(v=median( tmp),lty='dotted')
  
  lines( density(sample(MCMCvis::MCMCchains(samplesBinLinkBetaPriorPartialMode,params="omega")[,omegaSite[i]],3000)),col='red')
  
}
dev.off()

pdf(paste0(dirFigures,"appendix-x-match.pdf"), width=8, height=6)
v<-(1:200)[estimateComparisonDataArranged$`n()`==30][!is.na((1:200)[estimateComparisonDataArranged$`n()`==30])]

par(mfrow=c(2,4))
for(i in 1:8){
  plot(density(  MCMCvis::MCMCchains(samplesBinLinkBetaPriorComplete,params="p")[,v[i]] ),
       
       main="") 
  abline(v=median( MCMCvis::MCMCchains(samplesBinLinkBetaPriorComplete,params="p")[,v[i]] ))
  tmp <- sample(MCMCvis::MCMCchains(samplesBinLinkBetaPriorPartialMode,params="theta")[,v[i]],3000)
  lines( density(tmp),lty='dotted')
  abline(v=median( tmp),lty='dotted')
  
}
dev.off()

# summarize sites
summarySite = "LO"

summary.val <- samplesBinLinkBetaPriorPartialMode %>% 
  tidybayes::spread_draws(theta[site,year]) %>%
  dplyr::filter(site==summarySite) %>%
  tidyr::spread(year,theta) %>%
  dplyr::select(-c(site,.chain,.iteration,.draw))

omega.val <- samplesBinLinkBetaPriorPartialMode %>% 
  tidybayes::spread_draws(omega[site]) %>%
  dplyr::filter(site==summarySite) 
omega.val = omega.val$omega

summary.val = as.matrix(summary.val[,2:11])

d <- survDataAnalysisBayes %>%
  dplyr::filter(site==summarySite) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(y =sum(fruitingPlantNumber), n = sum(seedlingNumber)) %>%
  dplyr::mutate(prop = y/n)

pdf(file=paste0(dirFigures,"appendix-x-hierarchical",summarySite,".pdf"), width=8, height=6)

par(mfrow=c(2,5))
for(i in 1:nYears){
  hist(summary.val[,i],breaks=40,
       main=c(2006:2015)[i],
       xlab=expression("theta"),
       freq=FALSE)
  abline(v=d$prop[i],col='red',lty='dashed')
  lines(density(omega.val),col="red")
}

dev.off()


################################################################################
# JAGS fit for model with binomial likelihood, logit parameterization, partial pooling
################################################################################


summaryBinLinkLogitLinkPartial<-samplesBinLinkLogitLinkPartial %>% 
  tidybayes::spread_draws(alpha.i[site,year]) %>%
  dplyr::mutate(theta=boot::inv.logit(alpha.i)) %>%
  dplyr::summarise(bayesBBLpartial.q50=median(theta)) 

summaryBinLinkLogitLinkPartial$year <- as.character(summaryBinLinkLogitLinkPartial$year)
summaryBinLinkLogitLinkPartial$site <- as.factor(summaryBinLinkLogitLinkPartial$site)
summaryBinLinkLogitLinkPartial$year <- as.numeric(as.character(summaryBinLinkLogitLinkPartial$year))

estimatesComparisonData <-  estimatesComparisonData %>%
  dplyr::left_join(summaryBinLinkLogitLinkPartial,by=c("site","year"))

# plot comparing MLE fit and bayes fit 

pdf(paste0(dirFigures,"appendix-x-mle_bayeslogit.pdf"), width=8, height=4)
par(mfrow=c(1,2))
par(mar=c(5, 6, 4, 2) + 0.1)
plot(estimatesComparisonData$pHat,estimatesComparisonData$bayesBBLpartial.q50,
     xlab="Maximum likelihood estimate",
     ylab="Median of posterior for \n logit-link parameterization",
     xlim=c(0,1),ylim=c(0,1),
     pch=16,cex=.5)  
abline(a=0,b=1,lwd=0.5)

# calculate difference of posterior - max. likelihood estimate
post.chain.pars <-apply(MCMCvis::MCMCchains(samplesBinLinkLogitLinkPartial,params="alpha.i"),2,boot::inv.logit)
estimateComparisonDataArranged<-estimatesComparisonData %>%
  dplyr::arrange(year)
tmp<-matrix(NA,ncol=200,nrow=length(samplesBinLinkLogitLinkPartial)*dim(samplesBinLinkLogitLinkPartial[[1]]))
for(i in 1:200){
  tmp[,i]<-post.chain.pars[,i]-estimateComparisonDataArranged$pHat[i]
}
delta.conf <- apply(tmp,2,quantile,probs=c(.025,.5,.975))


par(mar=c(5, 6, 4, 2) + 0.1)
plot(estimateComparisonDataArranged$`n()`,t(delta.conf)[,2],
     pch=16,cex=0.5,
     ylim=c(-1,1),
     xlab="Sample size (n)",
     ylab="Comparison of the maximum likelihood estimate\n and the beta-binomial parameterization\n (median and 95% CI)")
segments(x0=estimateComparisonDataArranged$`n()`,y0=t(delta.conf)[,1],
         x1=estimateComparisonDataArranged$`n()`,y1=t(delta.conf)[,3])
abline(h=0,lwd=1)

dev.off()


sum <- summaryBinLinkBetaPriorComplete %>%
  dplyr::left_join(summaryBinLinkLogitLinkPartial,by=c("site","year"))

plot(sum$bayesBBcomp.q50,sum$bayesBBLpartial.q50)
abline(a=0,b=1)

tmp<-matrix(NA,ncol=200,nrow=length(samplesBinLinkLogitLinkPartial)*dim(samplesBinLinkLogitLinkPartial[[1]]))
for(i in 1:200){
  
  tmp[,i]<-MCMCvis::MCMCchains(samplesBinLinkBetaPriorComplete,params="p")[,i] - boot::inv.logit(MCMCvis::MCMCchains(samplesBinLinkLogitLinkPartial,params="alpha.i")[,i])
  
}

pdf(paste0(dirFigures,"appendix-x-bayescomp_bayeslogit.pdf"), width=4, height=4)
delta.conf <- apply(tmp,2,quantile,probs=c(.025,.5,.975))
par(mar=c(5, 6, 4, 2) + 0.1)
plot(estimateComparisonDataArranged$`n()`,t(delta.conf)[,2],
     pch=16,cex=0.5,
     ylim=c(-1,1),
     xlab="Sample size (n)",
     ylab="Comparison of the posteriors for \n complete and partial pooling parameterizations \n complete-partial (median and 95% CI)")
segments(x0=estimateComparisonDataArranged$`n()`,y0=t(delta.conf)[,1],
         x1=estimateComparisonDataArranged$`n()`,y1=t(delta.conf)[,3])
abline(h=0,lwd=1)
dev.off()


v<-(1:200)[estimateComparisonDataArranged$`n()`==1][!is.na((1:200)[estimateComparisonDataArranged$`n()`==1])]

omegaSite<-c(11,11,13,6,7,11,12,20)


par(mfrow=c(2,4))
for(i in 1:length(v)){
  plot(density(  MCMCvis::MCMCchains(samplesBinLinkBetaPriorComplete,params="p")[,v[i]] ),
       xlim=c(0,1),ylim=c(0,4),
       main="") 
  abline(v=median( MCMCvis::MCMCchains(samplesBinLinkBetaPriorComplete,params="p")[,v[i]] ))
  tmp <- boot::inv.logit(MCMCvis::MCMCchains(samplesBinLinkLogitLinkPartial,params="alpha.i"))[,v[i]]
  lines( density(tmp),lty='dotted')
  abline(v=median( tmp),lty='dotted')
  
} 
