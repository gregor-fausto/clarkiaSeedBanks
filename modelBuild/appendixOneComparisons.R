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
survivorshipFitFiles <- paste0(directory,list.files(directory))
dirFigures = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/")

samplesBinLinkBetaPriorComplete <- readRDS(survivorshipFitFiles[[1]])
samplesBinLinkBetaPriorPartialMean <- readRDS(survivorshipFitFiles[[2]])
samplesBinLinkBetaPriorPartialMode <- readRDS(survivorshipFitFiles[[3]])
samplesBinLinkLogitLinkPartial<- readRDS(survivorshipFitFiles[[4]])
survDataAnalysisBayes <- readRDS(survivorshipFitFiles[[5]])
mleDataLong <- readRDS(survivorshipFitFiles[[6]])

nYears = 10
nSites = 20

#MCMCsummary(samplesBinLinkBetaPriorPartialMean)
#MCMCsummary(samplesBinLinkBetaPriorPartialMode)

stat.summary <- function(x){
  
  x %<>% recover_types(survDataAnalysisBayes)
  
  tmp <- x %>% 
    tidybayes::spread_draws(theta[site,year]) %>%
    dplyr::summarise(ci.lo = quantile(theta, prob = c(.025)),
                     med = quantile(theta, prob = c(.5)),
                     ci.hi = quantile(theta, prob = c(.975)))
  return(tmp)
}

thetaBinLinkBetaPriorPartialMode <- stat.summary(samplesBinLinkBetaPriorPartialMode)
thetaBinLinkBetaPriorPartialMean <- stat.summary(samplesBinLinkBetaPriorPartialMean)

plot(thetaBinLinkBetaPriorPartialMode$med,thetaBinLinkBetaPriorPartialMean$med, 
     pch=16,cex = 0.5,
     xlim=c(0,1),ylim=c(0,1))
segments(x0=thetaBinLinkBetaPriorPartialMode$med,y0=thetaBinLinkBetaPriorPartialMean$ci.lo,
         x1=thetaBinLinkBetaPriorPartialMode$med,y1=thetaBinLinkBetaPriorPartialMean$ci.hi)
segments(x0=thetaBinLinkBetaPriorPartialMode$ci.lo,y0=thetaBinLinkBetaPriorPartialMean$med,
         x1=thetaBinLinkBetaPriorPartialMode$ci.hi,y1=thetaBinLinkBetaPriorPartialMean$med)
abline(a=0,b=1)

# calculate difference of posterior of mode - posterior of ean
post.chain.pars.mode <-MCMCvis::MCMCchains(samplesBinLinkBetaPriorPartialMode,params="theta")
post.chain.pars.mean <-MCMCvis::MCMCchains(samplesBinLinkBetaPriorPartialMean,params="theta")

# for plotting

samplesizeData <- survDataAnalysisBayes %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n())

samplesizeData$site <- as.factor(samplesizeData$site)
samplesizeData$year <- as.numeric(as.character(samplesizeData$year))

estimatesComparisonData <- mleDataLong %>%
  dplyr::left_join(samplesizeData,by=c("site","year"))

estimateComparisonDataArranged <- estimatesComparisonData %>%
  dplyr::arrange(year)

tmp<-matrix(NA,ncol=200,nrow=length(samplesBinLinkBetaPriorPartialMode)*dim(samplesBinLinkBetaPriorPartialMode[[1]]))
for(i in 1:200){
  tmp[,i]<-post.chain.pars.mode[,i]-post.chain.pars.mean[i]
}
delta.conf <- apply(tmp,2,quantile,probs=c(.025,.5,.975))

plot(estimateComparisonDataArranged$`n()`,t(delta.conf)[,2],
     pch=16,cex=0.5,
     ylim=c(-1,1),
     xlab="Sample size (n)",
     ylab="Comparison of the beta-binomial parameterized via mode\n and the beta-binomial parameterized via mean\n (median and 95% CI)")
segments(x0=estimateComparisonDataArranged$`n()`,y0=t(delta.conf)[,1],
         x1=estimateComparisonDataArranged$`n()`,y1=t(delta.conf)[,3])
abline(h=0)

estimateDistributionComparison<-data.frame(cbind(estimateComparisonDataArranged,t(delta.conf)))
estimateDistributionComparison <- estimateDistributionComparison %>%
  dplyr::select(-pHat)
names(estimateDistributionComparison) = c("site","year","n","ci.lo","med","ci.hi")

library(ggplot2)

# replicate graph above
ggplot(data=estimateDistributionComparison) +
  geom_point(aes(x=n,y=med)) +
  geom_errorbar(aes(x=n, ymin=ci.lo, ymax=ci.hi), width=.1) +
  theme_classic()

# facet by site
ggplot(data=estimateDistributionComparison) +
  geom_point(aes(x=n,y=med)) +
  geom_errorbar(aes(x=n, ymin=ci.lo, ymax=ci.hi), width=.1) +
  facet_wrap(~site) +
  geom_hline(yintercept=0) +
  theme_classic()

# facet by year
ggplot(data=estimateDistributionComparison) +
  geom_point(aes(x=n,y=med)) +
  geom_errorbar(aes(x=n, ymin=ci.lo, ymax=ci.hi), width=.1) +
  facet_wrap(~year) +
  geom_hline(yintercept=0) +
  theme_classic()
