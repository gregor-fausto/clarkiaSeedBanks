#################################################################################
################################################################################
################################################################################
# Code for comparing the following modeling approaches for the seedling survivorship data
# 1. binomial likelihood
# 2. binomial likelihood with beta prior, complete pooling
# 3. binomial likelihood with beta prior, partial pooling
# 4. binomial likelihood with logit parameterization
#
# Scripts by Gregor Siegmund
# fausto.siegmund@gmail.com
# last updated 02-29-2020
#################################################################################
#################################################################################
#################################################################################

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

#################################################################################
# get paths for JAGS script and figure directories
#################################################################################
dirJagsScripts = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/appendixOneScripts/")
dirFigures = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures")

################################################################################
# Load survivorship data
################################################################################
survivorshipData0615 <-
  readxl::read_excel(path = "~/Dropbox/dataLibrary/Clarkia-LTREB/20_demography_sites/Survivorship & Fecundity_06-15.xls")

survData <- survivorshipData0615 %>%
  janitor::clean_names(case = "lower_camel")

survVariables <- names(survData)
seedlingNumNames <-
  tidyselect::vars_select(survVariables, contains("seedlingNumber"))
fruitingPlantNumNames <-
  tidyselect::vars_select(survVariables, contains("fruitplNumber6"))

survDataSubset <- survData %>%
  dplyr::select(site,
                transect,
                position,
                seedlingNumNames,
                fruitingPlantNumNames)

seedlingData <- survDataSubset %>%
  dplyr::select(site, transect, position, seedlingNumNames) %>%
  tidyr::pivot_longer(cols = seedlingNumNames,
                      names_to = "year",
                      values_to = "seedlingNumber") %>%
  tidyr::separate(year, into = c("discard", "year")) %>%
  dplyr::mutate(year = as.numeric(paste0(20, year))) %>%
  dplyr::select(-discard)

fruitingplantData <- survDataSubset %>%
  dplyr::select(site, transect, position, fruitingPlantNumNames) %>%
  tidyr::pivot_longer(cols = fruitingPlantNumNames,
                      names_to = "year",
                      values_to = "fruitingPlantNumber") %>%
  tidyr::separate(year, into = c("discard", "year")) %>%
  dplyr::mutate(year = as.numeric(paste0(20, year))) %>%
  dplyr::select(-discard)

survDataReformatted <- seedlingData %>%
  dplyr::left_join(fruitingplantData, by = c("site", "transect", "position", "year")) %>%
  dplyr::rename(plot = position)

################################################################################
# Load site data
################################################################################
siteData <- read.csv(file="/Users/Gregor/Dropbox/dataLibrary/Datafile_Site_Environment_corr.csv",header=TRUE)

siteEastingData<-siteData %>%
  janitor::clean_names(case="lower_camel") %>%
  dplyr::select(site,easting) %>% unique %>%
  dplyr::filter(site %in% unique(survData$site))

################################################################################
# Prepare data for analysis
################################################################################
# filter dataset to rows in which
# the seedling number > 0 AND there are no NAs in the fruiting plant count
survDataAnalysis <- survDataReformatted %>%
    dplyr::filter(seedlingNumber>0 & !is.na(fruitingPlantNumber))

survDataAnalysis<-survDataAnalysis %>%
  dplyr::mutate(seedlingNumber = ifelse(seedlingNumber<fruitingPlantNumber,fruitingPlantNumber,seedlingNumber))

################################################################################
# Maximum likelihood estimate of survivorship, binomial likelihood
# by minimizing the negative log likelihood
################################################################################

# split the data frame into a list of lists by site and year
survDataAnalysisList <- split(survDataAnalysis,
                              list(survDataAnalysis$site,survDataAnalysis$year))


# create list (length 10, for 10 years) of lists (each length 20, for 20 sites)
survDataAnalysisListList <- list(survDataAnalysisList[1:20],
               survDataAnalysisList[21:40],
               survDataAnalysisList[41:60],
               survDataAnalysisList[61:80],
               survDataAnalysisList[81:100],
               survDataAnalysisList[101:120],
               survDataAnalysisList[121:140],
               survDataAnalysisList[141:160],
               survDataAnalysisList[161:180],
               survDataAnalysisList[181:200])

# empty matrix for maximum likelihood estimates
mles <- matrix(NA,20,10)
# empty matrix for number of trials
nTrials <- matrix(NA,20,10)

# obtain maximum likelihood estimates by minimizing the negative log likelihood for each year and site
for(i in 1:10){
  tmp.list <- survDataAnalysisListList[[i]]
  for(k in 1:20){
    tmp <- tmp.list[[k]]
    y = tmp$fruitingPlantNumber
    n = tmp$seedlingNumber
    negLogLik <- function(p,x=tmp$fruitingPlantNumber,size=tmp$seedlingNumber) -sum(dbinom(x=x, size=size, prob=p,log=TRUE))
    opt<-optimise(negLogLik,interval=c(0,1),maximum=FALSE)
    mles[k,i] <- opt$minimum
    nTrials[k,i] <- sum(tmp$seedlingNumber)
  }
}

# place maximum likelihood estimates in data frame
mleDataWide <- data.frame(site=unique(survDataAnalysis$site), mles )

names(mleDataWide)=c("site",2006:2015)

mleDataLong <- mleDataWide %>%
  tidyr::pivot_longer(cols=`2006`:`2015`,names_to="year",values_to="pHat")
mleDataLong$year<-as.numeric(mleDataLong$year)

# Figure
# Comparing the geographic pattern of MLE (red) to geographic pattern of all data
ggplot() +
  geom_point(data=survDataAnalysis%>% dplyr::left_join(siteEastingData,by="site"),aes(x=easting,y=fruitingPlantNumber/seedlingNumber),cex=0.25) +
  geom_point(data=data.frame(mleDataLong)%>% dplyr::left_join(siteEastingData,by="site"),aes(x=easting,y=pHat),color="red") +
  geom_smooth(data=survDataAnalysis%>% dplyr::left_join(siteEastingData,by="site"),aes(x=easting,y=fruitingPlantNumber/seedlingNumber),method="lm") +
  geom_smooth(data=data.frame(mleDataLong)%>% dplyr::left_join(siteEastingData,by="site"),aes(x=easting,y=pHat),color="red",method="lm") +
    facet_wrap(~year) +
  theme_minimal()

################################################################################
# Prepare data for fit with JAGS using tidybayes
################################################################################
survDataAnalysisBayes <- survDataAnalysis %>% 
  dplyr::select(-c(transect,plot)) %>%
  dplyr::mutate(year = as.factor(year))

survDataAnalysisBayesJags <- tidybayes::compose_data(survDataAnalysisBayes)

################################################################################
# JAGS fit for model with binomial likelihood, beta prior, complete pooling
################################################################################

################################################################################
# Prepare data for fit with JAGS using tidybayes
################################################################################
survDataAnalysisBayes <- survDataAnalysis %>% 
  dplyr::select(-c(transect,plot)) %>%
  dplyr::mutate(year = as.factor(year))

survDataAnalysisBayesJags <- tidybayes::compose_data(survDataAnalysisBayes)

################################################################################
# Initial values
################################################################################

nSites = survDataAnalysisBayesJags$n_site
nYears = survDataAnalysisBayesJags$n_year

inits = list(list(p = matrix(rep(.1,nSites*nYears),nrow=nSites,ncol=nYears)),
             list(p = matrix(rep(.5,nSites*nYears),nrow=nSites,ncol=nYears)),
             list(p = matrix(rep(.9,nSites*nYears),nrow=nSites,ncol=nYears))) 

################################################################################
# Set parameters to monitor
################################################################################

parsToMonitor = c("p")

################################################################################
# Set JAGS parameters
################################################################################

# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 500
n.update = 5000
n.iterations = 10000
n.thin = 1

################################################################################
# Fit model 
################################################################################

# set random seed
set.seed(2)

# tuning (n.adapt)
jm1 = jags.model(paste0(dirJagsScripts,"survivorshipModel-binLik-betaPrior-completePooling.R"), 
                 data = survDataAnalysisBayesJags, inits = inits,
                 n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm1, n.iterations = n.update)

# chain (n.iter)
samplesBinLinkBetaPriorComplete = coda.samples(jm1, variable.names = c(parsToMonitor), n.iter = n.iterations, thin = n.thin)

# MCMCvis::MCMCsummary(samples.rjags1)

samplesBinLinkBetaPriorComplete %<>% recover_types(survDataAnalysisBayes)

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
tmp<-matrix(NA,ncol=200,nrow=length(inits)*n.iterations)
for(i in 1:200){
  tmp[,i]<-post.chain.pars[,i]-estimateComparisonDataArranged$pHat[i]
}
delta.conf <- apply(tmp,2,quantile,probs=c(.025,.5,.975))

# plot comparing MLE fit and bayes fit 

setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
pdf(file="appendix-x-mle_bayes.pdf", width=8, height=4)

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

################################################################################
# Prepare data for fit with JAGS using tidybayes
################################################################################
survDataAnalysisBayes <- survDataAnalysis %>% 
  dplyr::select(-c(transect,plot)) %>%
  dplyr::mutate(year = as.factor(year))

survDataAnalysisBayesJags <- tidybayes::compose_data(survDataAnalysisBayes)

################################################################################
# Initial values
################################################################################

nSites = survDataAnalysisBayesJags$n_site
nYears = survDataAnalysisBayesJags$n_year

kappaInit = 100;

inits = list( list( theta = matrix(rep(.1,nSites*nYears),nrow=nSites,ncol=nYears),
                          omega=rep(.5,nSites) ,
                          kappaMinusTwo=rep(98,nSites) ),
              list( theta = matrix(rep(.5,nSites*nYears),nrow=nSites,ncol=nYears),
                    omega=rep(.5,nSites) ,
                    kappaMinusTwo=rep(98,nSites) ),
              list( theta = matrix(rep(.9,nSites*nYears),nrow=nSites,ncol=nYears),
                    omega=rep(.5,nSites) ,
                    kappaMinusTwo=rep(98,nSites) )
              )

################################################################################
# Set parameters to monitor
################################################################################

parsToMonitor = c("theta","omega","kappa")

################################################################################
# Set JAGS parameters
################################################################################

# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 500
n.update = 5000
n.iterations = 10000
n.thin = 1

################################################################################
# Fit model 
################################################################################

# set random seed
set.seed(2)

# tuning (n.adapt)
jm2 = jags.model(paste0(dirJagsScripts,"survivorshipModel-binLik-betaPrior-partialPooling.R"), 
                 data = survDataAnalysisBayesJags, inits = inits,
                 n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm2, n.iterations = n.update)

# chain (n.iter)
samplesBinLinkBetaPriorPartial = coda.samples(jm2, variable.names = c(parsToMonitor), 
                              n.iter = n.iterations, thin = n.thin)

#MCMCvis::MCMCsummary(samples.rjags3)

samplesBinLinkBetaPriorPartial %<>% recover_types(survDataAnalysisBayes)

summaryBinLinkBetaPriorPartial<-samplesBinLinkBetaPriorPartial %>% 
  tidybayes::spread_draws(theta[site,year]) %>%
  dplyr::summarise(bayesBBpartial.q50=median(theta)) 

summaryBinLinkBetaPriorPartial$year <- as.character(summaryBinLinkBetaPriorPartial$year)
summaryBinLinkBetaPriorPartial$site <- as.factor(summaryBinLinkBetaPriorPartial$site)
summaryBinLinkBetaPriorPartial$year <- as.numeric(as.character(summaryBinLinkBetaPriorPartial$year))

estimatesComparisonData<-  estimatesComparisonData %>%
  dplyr::left_join(summaryBinLinkBetaPriorPartial,by=c("site","year"))

# calculate difference of posterior - max. likelihood estimate
post.chain.pars <-MCMCvis::MCMCchains(samplesBinLinkBetaPriorPartial,params="theta")
estimateComparisonDataArranged<-estimatesComparisonData %>%
  dplyr::arrange(year)
tmp<-matrix(NA,ncol=200,nrow=length(inits)*n.iterations)
for(i in 1:200){
  tmp[,i]<-post.chain.pars[,i]-estimateComparisonDataArranged$pHat[i]
}
delta.conf <- apply(tmp,2,quantile,probs=c(.025,.5,.975))

# plot comparing MLE fit and bayes fit 

setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
pdf(file="appendix-x-mle_bayeshier.pdf", width=8, height=4)
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
  #dplyr::left_join(summaryBayes2,by=c("site","year")) %>%
  dplyr::left_join(summaryBinLinkBetaPriorPartial,by=c("site","year"))

plot(sum$bayesBBcomp.q50,sum$bayesBBpartial.q50)
abline(a=0,b=1)

tmp<-matrix(NA,ncol=200,nrow=30000)
for(i in 1:200){
tmp[,i]<-MCMCvis::MCMCchains(samplesBinLinkBetaPriorComplete,params="p")[,i] - MCMCvis::MCMCchains(samplesBinLinkBetaPriorPartial,params="theta")[,i]
}
delta.conf <- apply(tmp,2,quantile,probs=c(.025,.5,.975))

setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
pdf(file="appendix-x-bayescomp_bayeshier.pdf", width=4, height=4)
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


setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
pdf(file="appendix-x-mismatch.pdf", width=8, height=6)
v<-(1:200)[estimateComparisonDataArranged$`n()`==1][!is.na((1:200)[estimateComparisonDataArranged$`n()`==1])]

omegaSite<-c(11,11,13,6,7,11,12,20)

par(mfrow=c(2,4))
for(i in 1:length(v)){
plot(density(  MCMCvis::MCMCchains(samplesBinLinkBetaPriorComplete,params="p")[,v[i]] ),
     xlim=c(0,1),ylim=c(0,4),
     main="") 
  abline(v=median( MCMCvis::MCMCchains(samplesBinLinkBetaPriorComplete,params="p")[,v[i]] ))
tmp <- sample(MCMCvis::MCMCchains(samplesBinLinkBetaPriorPartial,params="theta")[,v[i]],3000)
lines( density(tmp),lty='dotted')
abline(v=median( tmp),lty='dotted')

lines( density(sample(MCMCvis::MCMCchains(samplesBinLinkBetaPriorPartial,params="omega")[,omegaSite[i]],3000)),col='red')

}
dev.off()

pdf(file="appendix-x-match.pdf", width=8, height=6)
v<-(1:200)[estimateComparisonDataArranged$`n()`==30][!is.na((1:200)[estimateComparisonDataArranged$`n()`==30])]

par(mfrow=c(2,4))
for(i in 1:8){
  plot(density(  MCMCvis::MCMCchains(samplesBinLinkBetaPriorComplete,params="p")[,v[i]] ),
       
       main="") 
  abline(v=median( MCMCvis::MCMCchains(samplesBinLinkBetaPriorComplete,params="p")[,v[i]] ))
  tmp <- sample(MCMCvis::MCMCchains(samplesBinLinkBetaPriorPartial,params="theta")[,v[i]],3000)
  lines( density(tmp),lty='dotted')
  abline(v=median( tmp),lty='dotted')

}
dev.off()

# summarize sites
summarySite = "LO"

summary.val <- samplesBinLinkBetaPriorPartial %>% 
  tidybayes::spread_draws(theta[site,year]) %>%
  dplyr::filter(site==summarySite) %>%
  tidyr::spread(year,theta) %>%
  dplyr::select(-c(site,.chain,.iteration,.draw))

omega.val <- samplesBinLinkBetaPriorPartial %>% 
  tidybayes::spread_draws(omega[site]) %>%
  dplyr::filter(site==summarySite) 
omega.val = omega.val$omega

summary.val = as.matrix(summary.val[,2:11])

d <- survDataAnalysisBayes %>%
  dplyr::filter(site==summarySite) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(y =sum(fruitingPlantNumber), n = sum(seedlingNumber)) %>%
  dplyr::mutate(prop = y/n)

setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
pdf(file=paste0("appendix-x-hierarchical",summarySite,".pdf"), width=8, height=6)

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
 
################################################################################
# Initial values
################################################################################
 
nSites = survDataAnalysisBayesJags$n_site
nYears = survDataAnalysisBayesJags$n_year
 
inits = list(list(mu.alpha = rep(rnorm(1),nSites), 
                  sigma.site = rep(rlnorm(1),nSites)),
             list(mu.alpha = rep(rnorm(1),nSites), 
                  sigma.site = rep(rlnorm(1),nSites)),
             list(mu.alpha = rep(rnorm(1),nSites),
                  sigma.site = rep(rlnorm(1),nSites))
) 


################################################################################
# Set parameters to monitor
################################################################################

parsToMonitor = c("mu.alpha","sigma.site","theta","alpha.i")

################################################################################
# Set JAGS parameters
################################################################################

# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 500
n.update = 5000
n.iterations = 10000
n.thin = 1

################################################################################
# Fit model 
################################################################################

# set random seed
set.seed(2)

# tuning (n.adapt)
jm3 = jags.model(paste0(dirJagsScripts,"survivorshipModel-binLik-logitLink-partialPooling.R"), 
                 data = survDataAnalysisBayesJags, inits = inits,
                 n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm3, n.iterations = n.update)

# chain (n.iter)
samplesBinLinkLogitLinkPartial = coda.samples(jm3, variable.names = c(parsToMonitor), 
                              n.iter = n.iterations, thin = n.thin)

#MCMCvis::MCMCsummary(samples.rjags4)

samplesBinLinkLogitLinkPartial %<>% recover_types(survDataAnalysisBayes)

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

setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
pdf(file="appendix-x-mle_bayeslogit.pdf", width=8, height=4)
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
tmp<-matrix(NA,ncol=200,nrow=length(inits)*n.iterations)
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

tmp<-matrix(NA,ncol=200,nrow=3000)
for(i in 1:200){
  
  tmp[,i]<-MCMCvis::MCMCchains(samplesBinLinkBetaPriorComplete,params="p")[,i] - boot::inv.logit(MCMCvis::MCMCchains(samplesBinLinkLogitLinkPartial,params="alpha.i")[,i])
  
}

setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
pdf(file="appendix-x-bayescomp_bayeslogit.pdf", width=4, height=4)
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
