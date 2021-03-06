
# log likelihood
# https://sourceforge.net/p/mcmc-jags/discussion/610036/thread/8211df61/

# 
# https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.43.903&rep=rep1&type=pdf

# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)
# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(rjags) # jags interface
library(MCMCvis)
library(dplyr)
#library(ggplot2)
library(reshape2)
library(tidyr)
library(HDInterval)
library(magrittr)
library(tidybayes)
library(parallel)
library(stringr)
library(loo)

#set.seed(10)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Read data
# -------------------------------------------------------------------
# -------------------------------------------------------------------

data = readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/seedBagMaster.RDS")

data = data %>%
  #dplyr::filter(site=="URS") %>%
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="-", remove=FALSE) %>%
  tidyr::unite(col='siteBag', c(site,bagNo), sep="-", remove=FALSE) %>%
  dplyr::mutate(siteBag = as.factor(siteBag))   %>%
  dplyr::mutate(year = as.factor(yearStart),
                age = as.factor(age)) %>%
  dplyr::select(siteBag,site,year,age,totalJan,seedStart,seedlingJan,intactOct) %>%
  dplyr::rename(siteBags = site,
                yearBags = year,
                ageBags = age)

refIndex <- data.frame(name=rep(c("totalJan","intactOct"),3),
                       ageBags=as.factor(c(1,1,2,2,3,3)),
                       betaIndex=c(1,2,3,4,5,6),
                       months=c(3,12,15,24,27,36),
                       g_1 = c(0,1,1,1,1,1),
                       g_2 = c(0,0,0,1,1,1),
                       g_3 = c(0,0,0,0,0,1)) %>%
  dplyr::mutate(months=months/36)

survivalData=data %>%
  dplyr::select(-seedlingJan) %>%
  tidyr::pivot_longer(cols=c(totalJan,intactOct)) %>%
  dplyr::left_join(refIndex,by=c("name","ageBags")) %>%
  dplyr::rename(response = name, y = value, siteSurvival = siteBags) %>%
  dplyr::mutate(gIndexSurvival = ageBags,
                yearSurvival = yearBags)

#survivalData=survivalData[!duplicated(survivalData$months),]
# survivalData=survivalData[!duplicated(interaction(survivalData$yearSurvival,
#                                                   survivalData$gIndexSurvival)),]

gCompIndex <- data.frame(yearSurvival = as.factor(c(rep(2006,6),rep(2007,4),rep(2008,2))),
                         ageBags=as.factor(c(1,1,2,2,3,3,1,1,2,2,1,1)),
                         compIndex=c(1:12),
                         response=rep(c("totalJan","intactOct"),6)) 

survivalData<-survivalData %>% 
  dplyr::left_join(gCompIndex,by=c('yearSurvival','ageBags',"response"))

#survivalData=survivalData[!duplicated(survivalData$compIndex),]


germinationData=data %>%
  dplyr::select(-c(seedStart,intactOct)) %>%
  dplyr::mutate(gIndex=ageBags) %>%
  dplyr::rename(siteGermination = siteBags,
                yearGermination = yearBags) %>%
  dplyr::mutate(seedlingJan2=seedlingJan)



# FECUNDITY ESTIMATES
censusSeedlingsFruitingPlants <- readRDS("~/Dropbox/dataLibrary/postProcessingData/censusSeedlingsFruitingPlants.RDS")

censusSeedlings <- censusSeedlingsFruitingPlants %>%
  dplyr::filter(year%in%c(2008,2009)) %>%
  dplyr::mutate(yearRef=year-1) %>%
  dplyr::rename(t1=year) %>%
  dplyr::select(-fruitplNumber)


countFruitsPerPlantTransects <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countFruitsPerPlantTransects.RDS")

countFruitsPerPlotTransects <- countFruitsPerPlantTransects %>%
  dplyr::filter(year%in%c(2007,2008)) %>%
  dplyr::group_by(site,year,transect,position,) %>%
  dplyr::summarise(totalFruitsPerPlot = sum(countFruitsPerPlant))

countSeedPerFruit <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countSeedPerFruit.RDS")

countSeedPerTotalFruitEquivalent <- countSeedPerFruit %>%
  dplyr::filter(demography==1) %>%
  dplyr::filter(year %in% c(2007,2008)) %>%
  dplyr::filter(damaged==0) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(meanSeedPerTotalFruitEquivalent = round(mean(sdno)))

dataUnknown<-countFruitsPerPlotTransects %>%
  dplyr::left_join(countSeedPerTotalFruitEquivalent,by=c("site","year")) %>%
  dplyr::mutate(fec.est = totalFruitsPerPlot*meanSeedPerTotalFruitEquivalent) %>%
  dplyr::mutate(t=year) %>% 
  dplyr::rename(yearRef=year) %>%
  dplyr::mutate(site=as.character(site),transect=as.character(transect),position=as.character(position)) %>%
  dplyr::left_join(censusSeedlings,by=c("site","transect","position","yearRef")) %>%
  dplyr::mutate(p = seedlingNumber/fec.est)

max(dataUnknown$p)

data.transect = dataUnknown %>%
  dplyr::group_by(site,yearRef,transect) %>%
  dplyr::summarise(fec.est = sum(fec.est), seedlingNumber = sum(seedlingNumber)) %>%
  dplyr::mutate(p = seedlingNumber/fec.est)
#plot(data.transect$fec.est,data.transect$seedlingNumber);abline(a=0,b=1)

data.transect = data.transect %>% 
  dplyr::rename(sitePlot = site) %>%
  dplyr::rename(yearPlot = yearRef) %>%
  dplyr::rename(fec = fec.est) %>%
  dplyr::rename(plotSeedlings = seedlingNumber) %>%
  dplyr::select(-p) %>%
  dplyr::mutate(yearPlot = as.factor(yearPlot)) #%>%
# dplyr::filter(sitePlot=="URS")





data <- tidybayes::compose_data(survivalData,germinationData,data.transect)

data$n1 = dim(survivalData)[1]
data$n2 = dim(germinationData)[1]
data$n3 = dim(data.transect)[1]

refDf<-data.frame(year=data$yearGermination,gIndex=data$gIndex) %>%
  unique %>%
  dplyr::mutate(germinationReference=1:6)

data$yearRefGerm = refDf$year
data$indexRefGerm = refDf$gIndex
data$refGerm = refDf$germinationReference

# d=survivalData %>%
#   dplyr::mutate(index=1:dim(survivalData)[1]) %>%
#   dplyr::group_by(months,yearBags) %>%
#   dplyr::top_n(n=1)  %>%
#   dplyr::ungroup() %>%
#   dplyr::rename(seedStart_pred=seedStart,
#                 siteSurvival_pred=siteSurvival,
#                 yearSurvival_pred=yearSurvival,
#                 months_pred=months,
#                 compIndex_pred=compIndex) %>%
#   dplyr::select(seedStart_pred,siteSurvival_pred,yearSurvival_pred,months_pred,compIndex_pred)
# 
# d
# d<-tidybayes::compose_data(d)
# names(d)[8]="n_pred"
# data=c(data,d)


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.chain = 3
n.adapt = 3000
n.update = 5000
n.iterations = 10000
n.thin = 1

dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/priorChecking/jagsScripts/")

# set inits for JAGS
inits <- list()

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Weibull residence time, age-dependent germination
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# # tuning (n.adapt)
#rjags::load.module("glm", quiet = TRUE)

jmWrAg = jags.model(paste0(dir,"binomialLikelihood-logLink-wrDecay-4.R"), data = data, n.adapt = n.adapt,
                    n.chains=3)

# burn-in (n.update)
update(jmWrAg, n.iter = n.update)

# samples.rjags = coda.samples(jmWrAg,
#                              variable.names =
#                                # GERMINATION
#                                c("mu0_g", "mu_g",
#                                  "sigma0_g","sigma_g",
#                                  #PRIOR PREDICTIVE
#                                  # "y_pred","y_prior.pred",
#                                  # SURVIVAL
#                                  "mu0_s","mu_s",
#                                 "sigma0_s","sigma_s",
#                                 "a", 
#                                 # PLOT DATA
#                                 "mu0_s0","mu_s0",
#                                 "sigma0_s0","sigma_s0"),
#                                 #  # LOG LIKELIHOODS
#                                 # "logLik_g","logLik_y","logLik_yplot",
#                                 # # POSTERIOR PREDICTIVE
#                                 # "y_sim","seedlingJan_sim","plotSeedlings_sim",
#                                 # # MODEL CHECKING 
#                                 # "chi2.plot.obs","chi2.plot.sim",
#                                 # "chi2.obs","chi2.sim",
#                                 # "chi2.yobs","chi2.ysim"
#                                 #),
#                              n.iter = n.iterations, thin = n.thin)

samples.rjags = coda.samples(jmWrAg,
                             variable.names =
                               c(# POSTERIOR PREDICTIVE
                                 "y_sim","seedlingJan_sim","plotSeedlings_sim",
                                 # MODEL CHECKING
                                 "chi2.plot.obs","chi2.plot.sim",
                                 "chi2.obs","chi2.sim",
                                 "chi2.yobs","chi2.ysim"),n.iter = n.iterations, thin = n.thin)


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Quick convergence checks
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# MCMCsummary(samples.rjags,params=c("mu0_g","mu0_s","mu0_s0"))
# MCMCsummary(samples.rjags,params=c("sigma0_g","sigma0_s","sigma0_s0"))
# options(max.print=999999)
# 
# MCMCsummary(samples.rjags,params=c("mu_g","mu_s","mu_s0"))
# MCMCsummary(samples.rjags,params=c("a"))
# # 
# hist(boot::inv.logit(MCMCchains(samples.rjags,params="mu0_s0")))
# hist(boot::inv.logit(MCMCchains(samples.rjags,params="mu_s0"))[,1])
# hist(boot::inv.logit(MCMCchains(samples.rjags,params="mu_s0"))[,2])
# 
# hist(boot::inv.logit(MCMCchains(samples.rjags,params="mu_g"))[,1])
# hist(boot::inv.logit(MCMCchains(samples.rjags,params="mu_g"))[,4])
# hist(boot::inv.logit(MCMCchains(samples.rjags,params="mu_g"))[,6])
# 
# mean(boot::inv.logit(MCMCchains(samples.rjags,params="mu0_g"))[,1])
# mean(boot::inv.logit(MCMCchains(samples.rjags,params="mu0_g"))[,2])
# mean(boot::inv.logit(MCMCchains(samples.rjags,params="mu_g"))[,6])

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Quick convergence checks
# -------------------------------------------------------------------
# -------------------------------------------------------------------


fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/decayModel/")
dir.create(file.path(fileDirectory), showWarnings = FALSE)
# 
samples.rjags=samples.rjags
saveRDS(samples.rjags,file=paste0(fileDirectory,"seedBurialSamplesChecks.rds"))
#saveRDS(data,file=paste0(fileDirectory,"data.rds"))
#saveRDS(seedBagsData,file=paste0(fileDirectory,"seedBagExperiment.rds"))
#saveRDS(viabilityRawData,file=paste0(fileDirectory,"viabilityExperiment.rds"))

# # -------------------------------------------------------------------
# # -------------------------------------------------------------------
# # How is germination doing?
# # -------------------------------------------------------------------
# # -------------------------------------------------------------------


# mcmcSamples = samples.rjags
# MCMCsummary(mcmcSamples,params="mu0_g")


# mu0_g=MCMCchains(mcmcSamples,params="mu0_g")
# mu_g=MCMCchains(mcmcSamples,params="mu_g")

# #pdf("~/Dropbox/clarkiaSeedBanks/products/figures/germination-population.pdf",width=8,height=6)

# par(mfrow=c(1,2))

# gamma1 = boot::inv.logit(mu0_g[,1:20])
# gamma1.sum=apply(gamma1,2,quantile,probs=c(0.025,.25,.5,.75,.975))

# plot(NA,NA,type='n',xlim=c(0,1),ylim=c(0,20),
     # axes=FALSE,frame=FALSE,
     # xlab="",ylab="")
# y.pt = 20:1
# for(i in 20:1){
  # tmp<-gamma1.sum[,i]
  # segments(x0=tmp[1],x1=tmp[5],y0=y.pt[i])
  # segments(x0=tmp[2],x1=tmp[4],y0=y.pt[i],lwd=3)
  
  # points(x=tmp[3],y=y.pt[i],pch=21,bg='white')
# }
# axis(1,  seq(0,1,by=.2), col.ticks = 1)
# axis(2, (1:20),
     # labels = rev(siteNames), las = 1, 
     # col = NA, col.ticks = 1, cex.axis = 1)
# mtext(expression(paste("Age 1 germination probability")),
      # side=1,line=2.5,adj=.5,col='black',cex=1)


# gamma1.sum<-data.frame(t(gamma1.sum),position)
# names(gamma1.sum)[1:5] = c("ci.lolo","ci.lo","ci.med","ci.hi","ci.hihi")

# plot(NA,NA,type='n',ylim=c(0,1),xlim=c(340,375),
     # axes=FALSE,frame=FALSE,
     # xlab="",ylab="")

# segments(x0=gamma1.sum$easting,y0=gamma1.sum$ci.lolo,y1=gamma1.sum$ci.hihi,lwd=1)
# segments(x0=gamma1.sum$easting,y0=gamma1.sum$ci.lo,y1=gamma1.sum$ci.hi,lwd=3)
# points(x=gamma1.sum$easting,y=gamma1.sum$ci.med,pch=21,bg='white')

# axis(2, seq(0,1,by=.2), col.ticks = 1)
# axis(1, seq(340,375,by=5),
     # labels = seq(340,375,by=5), las = 1, 
     # col.ticks = 1, cex.axis = 1)
# mtext("Age 1 germination probability",
      # side=2,line=2.5,adj=.5,col='black',cex=1)
# mtext("Easting (km)",
      # side=1,line=2.5,adj=.5,col='black',cex=1)

# par(mfrow=c(1,2))

# gamma2 = boot::inv.logit(mu0_g[,21:40])
# gamma2.sum=apply(gamma2,2,quantile,probs=c(0.025,.25,.5,.75,.975))

# gamma2.1 = boot::inv.logit(mu_g[,61:100])
# gamma2.sum.1=apply(gamma2.1,2,quantile,probs=c(0.025,.25,.5,.75,.975))


# plot(NA,NA,type='n',xlim=c(0,1),ylim=c(0,20),
     # axes=FALSE,frame=FALSE,
     # xlab="",ylab="")
# y.pt = 20:1
# for(i in 20:1){
  # tmp<-gamma2.sum[,i]
  # segments(x0=tmp[1],x1=tmp[5],y0=y.pt[i])
  # segments(x0=tmp[2],x1=tmp[4],y0=y.pt[i],lwd=3)
  
  # points(x=tmp[3],y=y.pt[i],pch=21,bg='white')
  
  
# }
# axis(1,  seq(0,1,by=.2), col.ticks = 1)
# axis(2, (1:20),
     # labels = rev(siteNames), las = 1, 
     # col = NA, col.ticks = 1, cex.axis = 1)
# mtext(expression(paste("Age 2 germination probability")),
      # side=1,line=2.5,adj=.5,col='black',cex=1)

# position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  # dplyr::select(site,easting) %>%
  # dplyr::mutate(easting=easting/1000)
# gamma2.sum<-data.frame(t(gamma2.sum),position)
# names(gamma2.sum)[1:5] = c("ci.lolo","ci.lo","ci.med","ci.hi","ci.hihi")

# plot(NA,NA,type='n',ylim=c(0,1),xlim=c(340,375),
     # axes=FALSE,frame=FALSE,
     # xlab="",ylab="")

# segments(x0=gamma2.sum$easting,y0=gamma2.sum$ci.lolo,y1=gamma2.sum$ci.hihi,lwd=1)
# segments(x0=gamma2.sum$easting,y0=gamma2.sum$ci.lo,y1=gamma2.sum$ci.hi,lwd=3)
# points(x=gamma2.sum$easting,y=gamma2.sum$ci.med,pch=21,bg='white')

# axis(2, seq(0,1,by=.2), col.ticks = 1)
# axis(1, seq(340,375,by=5),
     # labels = seq(340,375,by=5), las = 1, 
     # col.ticks = 1, cex.axis = 1)
# mtext("Age 2 germination probability",
      # side=2,line=2.5,adj=.5,col='black',cex=1)
# mtext("Easting (km)",
      # side=1,line=2.5,adj=.5,col='black',cex=1)

# par(mfrow=c(1,2))

# # note here this is the population-year level 
# # there is NO population-level parameter for age 3 germination
# gamma3 = boot::inv.logit(mu_g[,101:120])
# gamma3.sum=apply(gamma3,2,quantile,probs=c(0.025,.25,.5,.75,.975))

# plot(NA,NA,type='n',xlim=c(0,1),ylim=c(0,20),
     # axes=FALSE,frame=FALSE,
     # xlab="",ylab="")
# y.pt = 20:1
# for(i in 20:1){
  # tmp<-gamma3.sum[,i]
  # segments(x0=tmp[1],x1=tmp[5],y0=y.pt[i])
  # segments(x0=tmp[2],x1=tmp[4],y0=y.pt[i],lwd=3)
  
  # points(x=tmp[3],y=y.pt[i],pch=21,bg='white')
# }
# axis(1,  seq(0,1,by=.2), col.ticks = 1)
# axis(2, (1:20),
     # labels = rev(siteNames), las = 1, 
     # col = NA, col.ticks = 1, cex.axis = 1)
# mtext(expression(paste("Age 3 germination probability")),
      # side=1,line=2.5,adj=.5,col='black',cex=1)


# gamma3.sum<-data.frame(t(gamma3.sum),position)
# names(gamma3.sum)[1:5] = c("ci.lolo","ci.lo","ci.med","ci.hi","ci.hihi")

# plot(NA,NA,type='n',ylim=c(0,1),xlim=c(340,375),
     # axes=FALSE,frame=FALSE,
     # xlab="",ylab="")

# segments(x0=gamma3.sum$easting,y0=gamma3.sum$ci.lolo,y1=gamma3.sum$ci.hihi,lwd=1)
# segments(x0=gamma3.sum$easting,y0=gamma3.sum$ci.lo,y1=gamma3.sum$ci.hi,lwd=3)
# points(x=gamma3.sum$easting,y=gamma3.sum$ci.med,pch=21,bg='white')

# axis(2, seq(0,1,by=.2), col.ticks = 1)
# axis(1, seq(340,375,by=5),
     # labels = seq(340,375,by=5), las = 1, 
     # col.ticks = 1, cex.axis = 1)
# mtext("Age 3 germination probability",
      # side=2,line=2.5,adj=.5,col='black',cex=1)
# mtext("Easting (km)",
      # side=1,line=2.5,adj=.5,col='black',cex=1)