# -------------------------------------------------------------------
# Models for germination in seed bag burial experiment
# as a function of winter climate 
# -------------------------------------------------------------------
# Complete pooling model: 
# With low precipitation, increased temperature promotes germination
# At high precipitation, increased temperature reduces germination
# At low temperatures, increased precipitation promotes germination
# At high temperatures, increased precipitation reduces germination
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# No pooling model: 
# With low precipitation, increased temperature promotes germination
# At high precipitation, increased temperature reduces germination
# At low temperatures, increased precipitation promotes germination
# At high temperatures, increased precipitation reduces germination
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# log likelihood
# https://sourceforge.net/p/mcmc-jags/discussion/610036/thread/8211df61/
# https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.43.903&rep=rep1&type=pdf
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
# rm(list=setdiff(ls(all=TRUE),c("dataDirectory","modelDirectory","fileDirectory","n.adapt","n.update","n.iterations","n.thin"))) # if using in source(script)
options(stringsAsFactors = FALSE)
# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(rjags) 
library(tidybayes)
library(tidyverse)
library(parallel)

# -------------------------------------------------------------------
# Set directories
# -------------------------------------------------------------------
dataDirectory = "~/Dropbox/dataLibrary/postProcessingData/"
modelDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/"

# -------------------------------------------------------------------
# Import seed bag burial data
# -------------------------------------------------------------------
seedBagsData = readRDS(paste0(dataDirectory,"seedBagsData.rds"))

# -------------------------------------------------------------------
# Manually relabel a bag number and filter out 1 row of data
# -------------------------------------------------------------------
# relabel bag
seedBagsData[seedBagsData$site=="DLW" & seedBagsData$bagNo==43 & seedBagsData$age==1,]$bagNo = 42

# remove row of data with 108 seeds in January
seedBagsData <- seedBagsData %>% dplyr::filter(totalJan<=100)

# -------------------------------------------------------------------
# Prep data for JAGS
# -------------------------------------------------------------------
# Create a variable for the number of seeds at the start of the trial
seedBagsData$seedStart<-as.double(100)

# Rename data frame
data = seedBagsData

# rename variables and oragnize
data = data %>%
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="-", remove=FALSE) %>%
  tidyr::unite(col='siteBag', c(site,bagNo), sep="-", remove=FALSE) %>%
  dplyr::mutate(siteBag = as.factor(siteBag))   %>%
  dplyr::mutate(year = as.factor(yearStart),
                age = as.factor(age)) %>%
  dplyr::select(siteBag,site,year,age,totalJan,seedStart,seedlingJan,intactOct) %>%
  dplyr::rename(siteBags = site,
                yearBags = year,
                ageBags = age)

# create a reference dataframe
refIndex <- data.frame(name=rep(c("totalJan","intactOct"),3),
                       ageBags=as.factor(c(1,1,2,2,3,3)),
                       betaIndex=c(1,2,3,4,5,6),
                       months=c(3,12,15,24,27,36),
                       g_1 = c(0,1,1,1,1,1),
                       g_2 = c(0,0,0,1,1,1),
                       g_3 = c(0,0,0,0,0,1)) %>%
  dplyr::mutate(months=months/36)

# create a reference data frame for event histories
gCompIndex <- data.frame(yearSurvival = as.factor(c(rep(2006,6),rep(2007,4),rep(2008,2))),
                         ageBags=as.factor(c(1,1,2,2,3,3,1,1,2,2,1,1)),
                         compIndex=c(1:12),
                         response=rep(c("totalJan","intactOct"),6)) 

# create a data frame for germination data
germinationData=data %>%
  dplyr::select(-c(seedStart,intactOct)) %>%
  dplyr::mutate(gIndex=ageBags) %>%
  dplyr::rename(siteGermination = siteBags,
                yearGermination = yearBags) %>%
  dplyr::mutate(seedlingJan2=seedlingJan)

# indexing to avoid issue of dealing with ragged arrays in JAGS
germinationAge = c(1,2,3,1,2,1)
germinationYear = c(2006,2006,2006,2007,2007,2008)
germinationIndex = 1:6

df.index=data.frame(ageBags=as.factor(germinationAge),yearGermination=as.factor(germinationYear),germinationIndex=germinationIndex)

germinationData <- germinationData %>%
  dplyr::left_join(df.index,by=c("ageBags","yearGermination")) %>%
  dplyr::mutate(yearGermination=as.factor(yearGermination),ageBags=as.factor(ageBags))

# remove data with missing totalJan from germination dataset
germinationData <- germinationData %>%
  dplyr::filter(!is.na(totalJan))

germinationData = germinationData %>% dplyr::filter(ageBags==1)

# -------------------------------------------------------------------
# Read in climate data
# -------------------------------------------------------------------
climate=readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/climateData.RDS")
climate=climate %>% dplyr::ungroup()
climate <- climate %>% dplyr::filter(intenseDemography==1)
climate <- climate %>% dplyr::filter(year %in% 2005:2007&season=="winter") %>% dplyr::mutate(year = year+1)

climate = germinationData %>%
  dplyr::rename(site = siteGermination , year = yearGermination) %>%
  dplyr::mutate(site=as.character(site),year=as.numeric(as.character(year))) %>%
  dplyr::left_join(climate,by=c("site","year"))

muTemp = mean(climate$t,na.rm=TRUE)
sdTemp = sd(climate$t,na.rm=TRUE)
muPrecip = mean(climate$p,na.rm=TRUE)
sdPrecip = sd(climate$p,na.rm=TRUE)

climate$t.std = (climate$t - muTemp)/sdTemp
climate$p.std = (climate$p - muPrecip)/sdPrecip

# -------------------------------------------------------------------
# Visualize data
# -------------------------------------------------------------------
g1 = ggplot(climate) +
  geom_point(aes(x= t, y = seedlingJan/totalJan) , alpha = 3/10, shape = 21, colour = "black",  fill = "brown", size = 3) +
  theme_minimal()

g2 =ggplot(climate) +
  geom_point(aes(x= p, y = seedlingJan/totalJan) , alpha = 3/10, shape = 21, colour = "black",  fill = "brown", size = 3) +
  theme_minimal()

# -------------------------------------------------------------------
# Predictions
# -------------------------------------------------------------------
# for plotting probability of germination as function of temperature
temp_pred = seq(min(climate$t), max(climate$t), length.out=20) 
#Standardized for making predictions of probability of occupancy at new
temp_pred_std = (temp_pred - muTemp) / sdTemp

# for plotting probability of germination as function of temperature
precip_pred = seq(min(climate$p), max(climate$p), length.out=20) 
#Standardized for making predictions of probability of occupancy at new
precip_pred_std = (precip_pred - muPrecip) / sdPrecip

# -------------------------------------------------------------------
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------

data <- tidybayes::compose_data(germinationData)

data$n = dim(germinationData)[1]

# add cliamte data
data$clim.std = climate$t.std
clim_pred_std = precip_pred_std
muClim = muTemp
sdClim = sdTemp
gClim = g2

detach("package:tidyverse", unload=TRUE)

# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
# -------------------------------------------------------------------

n.adapt = 3000
n.update = 5000
n.iterations = 5000
n.thin = 1

# -------------------------------------------------------------------
# Describe functions
# -------------------------------------------------------------------

initsMu0 <- function(samps = 1){
  rnorm(n = samps, mean = 0, sd = 1)
}

initsSigma0Weak <- function(samps = 1){
  extraDistr::rhnorm(n = samps, sigma = 1)
}

initsSigmaWeakMat <- function(samps = data$n_siteGermination){
  extraDistr::rhnorm(n = samps, sigma = 1)
}

initsBeta <- function(samps = 1){
  rnorm(n = samps, mean = 0, sd = 1)
}

# -------------------------------------------------------------------
# Set inits for JAGS
# -------------------------------------------------------------------
inits <- list()
for(i in 1:3){
  inits[[i]] <- list(initsMu0(),initsBeta())
  
  names(inits[[i]]) = c("alpha", "beta")
  
}

# -------------------------------------------------------------------
# Call to JAGS
# -------------------------------------------------------------------
jm = jags.model("~/Dropbox/clarkiaSeedBanks/modelBuild/jags-germinationInteractions-completePool.R",
                data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

update(jm, n.iter = n.update)

parsToMonitor = c("alpha","beta")

zm.completePool = coda.samples(jm, variable.names = c(parsToMonitor),
                               n.iter = n.iterations, thin = 1)

# -------------------------------------------------------------------
# Check convergence
# -------------------------------------------------------------------
MCMCvis::MCMCsummary(zm.completePool,c("alpha","beta"))

library(MCMCvis)
MCMCplot(zm.completePool, params = c("alpha","beta"))

#----------------------------------------------------------------
# Plot parameter estimates
# -------------------------------------------------------------------
alpha = as.vector(MCMCvis::MCMCchains(zm.completePool,params="alpha"))
beta = MCMCvis::MCMCchains(zm.completePool,params="beta")

preds=alpha+outer(beta[,1],clim_pred_std)
preds=apply(preds,2,boot::inv.logit)
pred3 <- data.frame(t(apply(preds,2,hdi,.95)))
colnames(pred3) = c("lower","upper")
pred4 <- apply(preds,2,median)
pred.po.df <- cbind(clim_pred_std, pred3, median = pred4)

g4 <- gClim +
  geom_line(data = pred.po.df, aes(x = clim_pred_std*sdClim+muClim, y = median),col='red') +
  geom_ribbon(data = pred.po.df, aes(x = clim_pred_std*sdClim+muClim, ymin = lower, ymax = upper), alpha = 0.2, fill = "red")
g4

g4 + facet_wrap(~site)


# -------------------------------------------------------------------
# Describe functions
# -------------------------------------------------------------------

initsMu0 <- function(samps = 1){
  rnorm(n = samps, mean = 0, sd = 1)
}

initsBeta <- function(samps = 1){
  rnorm(n = samps, mean = 0, sd = 1)
}

# -------------------------------------------------------------------
# Set inits for JAGS
# -------------------------------------------------------------------
inits <- list()
for(i in 1:3){
  inits[[i]] <- list(initsMu0(samps=20),initsBeta())
  
  names(inits[[i]]) = c("alpha", "beta")
  
}

# -------------------------------------------------------------------
# Call to JAGS
# -------------------------------------------------------------------
jm = jags.model("~/Dropbox/clarkiaSeedBanks/modelBuild/jags-germinationInteractions-noPool.R",
                data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

update(jm, n.iter = n.update)

parsToMonitor = c("alpha","beta")

zm.noPool = coda.samples(jm, variable.names = c(parsToMonitor),
                         n.iter = n.iterations, thin = 1)

# -------------------------------------------------------------------
# Check convergence
# -------------------------------------------------------------------
MCMCvis::MCMCsummary(zm.noPool,c("alpha","beta"))

library(MCMCvis)
MCMCplot(zm.noPool, params = c("beta"))
MCMCplot(zm.noPool, params = c("alpha"),horiz=FALSE)

alpha = MCMCchains(zm.noPool,"alpha")
beta = MCMCchains(zm.noPool,"beta")

plot.list.noPool=list()

for(j in 1:20){
  preds=alpha[,j]+outer(beta[,1],clim_pred_std)
  preds=apply(preds,2,boot::inv.logit)
  pred3 <- apply(preds,2,hdi,.95)
  pred4 <- apply(preds,2,median)
  pred.po.df <- cbind(clim_pred_std, data.frame(t(pred3)), median = pred4)
  
  plot.list.noPool[[j]] = pred.po.df
}

siteNames = unique(climate$site)
par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  tmp = climate[climate$site==siteNames[i],]
  plot(tmp$t.std,tmp$seedlingJan/tmp$totalJan,xlim=c(min(clim_pred_std),max(clim_pred_std)),ylim=c(0,1))
  lines(plot.list.noPool[[i]]$clim_pred_std,plot.list.noPool[[i]]$median,type='l',lwd=2)
}


# -------------------------------------------------------------------
# Random intercept model
# -------------------------------------------------------------------

initsMu0 <- function(samps = 1){
  rnorm(n = samps, mean = 0, sd = 1)
}

initsSigma0Weak <- function(samps = 1){
  extraDistr::rhnorm(n = samps, sigma = 1)
}

initsSigmaWeakMat <- function(samps = data$n_siteGermination){
  extraDistr::rhnorm(n = samps, sigma = 1)
}

initsBeta <- function(samps = 1){
  rnorm(n = samps, mean = 0, sd = 1)
}

# -------------------------------------------------------------------
# Set inits for JAGS
# -------------------------------------------------------------------
inits <- list()
for(i in 1:3){
  inits[[i]] <- list(initsMu0(),initsSigma0Weak(),initsBeta())
  
  names(inits[[i]]) = c("mu.alpha","sigma.alpha", "beta")
  
}

# -------------------------------------------------------------------
# Call to JAGS
# -------------------------------------------------------------------
jm = jags.model("~/Dropbox/clarkiaSeedBanks/modelBuild/jags-germinationInteractions-randomIntercept.R",
                data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

update(jm, n.iter = n.update)

parsToMonitor = c("mu.alpha","sigma.alpha","alpha","beta")

zm.randomIntercept = coda.samples(jm, variable.names = c(parsToMonitor),
                                  n.iter = n.iterations, thin = 1)

# -------------------------------------------------------------------
# Check convergence
# -------------------------------------------------------------------
MCMCvis::MCMCsummary(zm.randomIntercept, c("mu.alpha","sigma.alpha","alpha","beta"))

dev.off()
library(MCMCvis)
MCMCplot(zm.randomIntercept, params = c("beta"))
MCMCplot(zm.randomIntercept, params = c("alpha"),horiz=FALSE)

mu.alpha = MCMCchains(zm.randomIntercept,"mu.alpha")
sigma.alpha = MCMCchains(zm.randomIntercept,"sigma.alpha")

beta = MCMCchains(zm.randomIntercept,"beta")

alpha = rnorm(dim(mu.alpha)[1],mu.alpha, sigma.alpha)

preds=alpha+outer(beta[,1],clim_pred_std)
preds=apply(preds,2,boot::inv.logit)
pred3 <- apply(preds,2,HDInterval::hdi,.95)
pred4 <- apply(preds,2,median)
pred.po.df <- cbind(clim_pred_std, data.frame(t(pred3)), median = pred4)

g4 <- g2 +
  geom_line(data = pred.po.df, aes(x = clim_pred_std*sdClim+muClim, y = median),col='red') +
  geom_ribbon(data = pred.po.df, aes(x = clim_pred_std*sdClim+muClim, ymin = lower, ymax = upper), alpha = 0.2, fill = "red")
g4

plot.listRanInt=list()

alpha = MCMCchains(zm.randomIntercept,"alpha")
beta = MCMCchains(zm.randomIntercept,"beta")

for(j in 1:20){
  preds=alpha[,j]+outer(beta[,1],clim_pred_std)
  preds=apply(preds,2,boot::inv.logit)
  pred3 <- apply(preds,2,hdi,.95)
  pred4 <- apply(preds,2,median)
  pred.po.df <- cbind(clim_pred_std, data.frame(t(pred3)), median = pred4)
  
  plot.listRanInt[[j]] = pred.po.df
}

siteNames = unique(climate$site)
par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  tmp = climate[climate$site==siteNames[i],]
  plot(tmp$t.std,tmp$seedlingJan/tmp$totalJan,xlim=c(min(clim_pred_std),max(clim_pred_std)),ylim=c(0,1))
  lines(plot.listRanInt[[i]]$clim_pred_std,plot.listRanInt[[i]]$median,type='l',lwd=2)
  lines(plot.list.noPool[[i]]$clim_pred_std,plot.list.noPool[[i]]$median,type='l',lwd=2,lty='dotted',col='orange')
}

siteNames = unique(climate$site)
par(mfrow=c(1,1),mar=c(2,2,2,2))
plot(NA,NA,xlim=c(min(clim_pred_std),max(clim_pred_std)),ylim=c(0,1))
for(i in 1:20){
  tmp = climate[climate$site==siteNames[i],]
  lines(plot.listRanInt[[i]]$clim_pred_std,plot.listRanInt[[i]]$median,type='l',lwd=2)
  lines(plot.list.noPool[[i]]$clim_pred_std,plot.list.noPool[[i]]$median,type='l',lwd=2,lty='dotted',col='orange')
}


alpha.noPool = boot::inv.logit(MCMCchains(zm.noPool,params="alpha"))
alpha.ranInt = boot::inv.logit(MCMCchains(zm.randomIntercept,params="alpha"))
alpha.gm = boot::inv.logit(MCMCchains(zm.randomIntercept,params="mu.alpha"))

a.np=apply(alpha.noPool,2,quantile,c(.025,.5,.975))
a.ri=apply(alpha.ranInt,2,quantile,c(.025,.5,.975))
a.gm=apply(alpha.gm,2,quantile,c(.025,.5,.975))

par(mfrow=c(1,1))
plot(a.np[2,],a.ri[2,],xlim=c(.075,.3),ylim=c(.075,.3));abline(a=0,b=1)
segments(x0=a.np[2,],y0=a.ri[1,],y1=a.ri[3,])
segments(y0=a.ri[2,],x0=a.np[1,],x1=a.np[3,])
abline(h=a.gm[2,],lty='dotted')


position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)

par(mfrow=c(1,1))
plot(position$easting,a.ri[2,],ylim=c(.075,.3),pch=19);abline(a=0,b=1)
segments(x0=position$easting,y0=a.ri[1,],y1=a.ri[3,])


# -------------------------------------------------------------------
# Random intercept model + elevation
# -------------------------------------------------------------------

elevation<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,elevation)
elev = elevation$elevation
sdElevation = sd(elev)
muElevation = mean(elev)
elev.std = (elev-mean(elev))/sd(elev)

# or  climate
climate2=readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/climateData.RDS")
climate2=climate2 %>% dplyr::ungroup()
climate2 <- climate2 %>% dplyr::filter(intenseDemography==1)
climate2 <- climate2 %>% dplyr::filter(season=="winter") 
winter_temp=climate2 %>% dplyr::group_by(site) %>% dplyr::summarise(mu.p = mean(t))
wt = winter_temp$mu.p
elev = wt
sdElevation = sd(wt)
muElevation = mean(wt)
elev.std = (wt-mean(wt))/sd(wt)


data$elev.std = elev.std

initsMu0 <- function(samps = 1){
  rnorm(n = samps, mean = 0, sd = 1)
}

initsSigma0Weak <- function(samps = 1){
  extraDistr::rhnorm(n = samps, sigma = 1)
}

initsSigmaWeakMat <- function(samps = data$n_siteGermination){
  extraDistr::rhnorm(n = samps, sigma = 1)
}

initsBeta <- function(samps = 1){
  rnorm(n = samps, mean = 0, sd = 1)
}

# -------------------------------------------------------------------
# Set inits for JAGS
# -------------------------------------------------------------------
inits <- list()
for(i in 1:3){
  inits[[i]] <- list(initsMu0(),initsMu0(),initsSigma0Weak(),initsBeta())
  
  names(inits[[i]]) = c("kappa","eta","sigma.alpha", "beta")
  
}

# -------------------------------------------------------------------
# Call to JAGS
# -------------------------------------------------------------------
jm = jags.model("~/Dropbox/clarkiaSeedBanks/modelBuild/jags-germinationInteractions-randomInterceptGroup.R",
                data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

update(jm, n.iter = n.update)

parsToMonitor = c("kappa","eta","sigma.alpha","alpha","beta")

zm.randomInterceptGroup = coda.samples(jm, variable.names = c(parsToMonitor),
                                  n.iter = n.iterations, thin = 1)

# -------------------------------------------------------------------
# Check convergence
# -------------------------------------------------------------------
MCMCvis::MCMCsummary(zm.randomInterceptGroup, c("kappa","eta","sigma.alpha","alpha","beta"))

dev.off()
library(MCMCvis)
MCMCplot(zm.randomInterceptGroup, params = c("beta"))
MCMCplot(zm.randomInterceptGroup, params = c("alpha"),horiz=FALSE)

sigma.alpha = MCMCchains(zm.randomInterceptGroup,"sigma.alpha")

kappa = MCMCchains(zm.randomInterceptGroup,"kappa")
eta = MCMCchains(zm.randomInterceptGroup,"eta")
mu.alpha = kappa
beta = MCMCchains(zm.randomInterceptGroup,"beta")

alpha = rnorm(dim(mu.alpha)[1],mu.alpha, sigma.alpha)

preds=alpha+outer(beta[,1],clim_pred_std)
preds=apply(preds,2,boot::inv.logit)
pred3 <- apply(preds,2,HDInterval::hdi,.95)
pred4 <- apply(preds,2,median)
pred.po.df <- cbind(clim_pred_std, data.frame(t(pred3)), median = pred4)

g4 <- g2 +
  geom_line(data = pred.po.df, aes(x = clim_pred_std*sdClim+muClim, y = median),col='red') +
  geom_ribbon(data = pred.po.df, aes(x = clim_pred_std*sdClim+muClim, ymin = lower, ymax = upper), alpha = 0.2, fill = "red")
g4

plot.listRanIntGroup=list()

alpha = MCMCchains(zm.randomInterceptGroup,"alpha")
beta = MCMCchains(zm.randomInterceptGroup,"beta")

for(j in 1:20){
  preds=alpha[,j]+outer(beta[,1],clim_pred_std)
  preds=apply(preds,2,boot::inv.logit)
  pred3 <- apply(preds,2,hdi,.95)
  pred4 <- apply(preds,2,median)
  pred.po.df <- cbind(clim_pred_std, data.frame(t(pred3)), median = pred4)
  
  plot.listRanIntGroup[[j]] = pred.po.df
}

siteNames = unique(climate$site)
par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  tmp = climate[climate$site==siteNames[i],]
  plot(tmp$t.std,tmp$seedlingJan/tmp$totalJan,xlim=c(min(clim_pred_std),max(clim_pred_std)),ylim=c(0,1))
  lines(plot.listRanIntGroup[[i]]$clim_pred_std,plot.listRanIntGroup[[i]]$median,type='l',lwd=2)
  lines(plot.listRanInt[[i]]$clim_pred_std,plot.listRanInt[[i]]$median,type='l',lwd=2)
  lines(plot.list.noPool[[i]]$clim_pred_std,plot.list.noPool[[i]]$median,type='l',lwd=2,lty='dotted',col='orange')
}

siteNames = unique(climate$site)
par(mfrow=c(1,1),mar=c(2,2,2,2))
plot(NA,NA,xlim=c(min(clim_pred_std),max(clim_pred_std)),ylim=c(0,1))
for(i in 1:20){
  tmp = climate[climate$site==siteNames[i],]
  lines(plot.listRanInt[[i]]$clim_pred_std,plot.listRanInt[[i]]$median,type='l',lwd=2)
  lines(plot.list.noPool[[i]]$clim_pred_std,plot.list.noPool[[i]]$median,type='l',lwd=2,lty='dotted',col='orange')
}

alpha.ranIntGroup = boot::inv.logit(MCMCchains(zm.randomInterceptGroup,params="alpha"))

alpha.noPool = boot::inv.logit(MCMCchains(zm.noPool,params="alpha"))
alpha.ranInt = boot::inv.logit(MCMCchains(zm.randomIntercept,params="alpha"))
alpha.gm = boot::inv.logit(MCMCchains(zm.randomIntercept,params="mu.alpha"))

a.rig=apply(alpha.ranIntGroup,2,quantile,c(.025,.5,.975))
a.np=apply(alpha.noPool,2,quantile,c(.025,.5,.975))
a.ri=apply(alpha.ranInt,2,quantile,c(.025,.5,.975))
a.gm=apply(alpha.gm,2,quantile,c(.025,.5,.975))

par(mfrow=c(1,1))
plot(a.np[2,],a.ri[2,],xlim=c(.075,.3),ylim=c(.075,.3));abline(a=0,b=1)
segments(x0=a.np[2,],y0=a.ri[1,],y1=a.ri[3,])
segments(y0=a.ri[2,],x0=a.np[1,],x1=a.np[3,])
abline(h=a.gm[2,],lty='dotted')

par(mfrow=c(1,1))
plot(a.rig[2,],a.ri[2,],xlim=c(.075,.3),ylim=c(.075,.3));abline(a=0,b=1)
segments(x0=a.rig[2,],y0=a.ri[1,],y1=a.ri[3,])
segments(y0=a.ri[2,],x0=a.rig[1,],x1=a.rig[3,])
abline(h=a.gm[2,],lty='dotted')

position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)

par(mfrow=c(1,1))
plot(position$easting,a.rig[2,],ylim=c(.075,.3),pch=19);abline(a=0,b=1)
segments(x0=position$easting,y0=a.rig[1,],y1=a.rig[3,])


elev_pred=seq(min(elev.std),max(elev.std),length.out=50)
mu.alpha=as.vector(kappa) + eta %*% elev_pred
tmp=apply(boot::inv.logit(mu.alpha),2,quantile,c(.025,.5,.975))
elev_pred = elev_pred*sdElevation+muElevation
plot(elev_pred,tmp[2,],type='l',ylim=c(0,.5))
polygon(c(elev_pred,rev(elev_pred)),c(tmp[1,],rev(tmp[3,])),border=0,col='gray95')
lines(elev_pred,tmp[2,],lwd=2)

a.ri=apply(alpha.ranInt,2,quantile,c(.025,.5,.975))
points(winter_temp$mu.p,a.ri[2,],pch=19)
segments(x0=winter_temp$mu.p,y0=a.ri[1,],y1=a.ri[3,])



# -------------------------------------------------------------------
# Random intercept+random slope model
# -------------------------------------------------------------------

initsMu0 <- function(samps = 1){
  rnorm(n = samps, mean = 0, sd = 1)
}

initsSigma0Weak <- function(samps = 1){
  extraDistr::rhnorm(n = samps, sigma = 1)
}

initsSigmaWeakMat <- function(samps = data$n_siteGermination){
  extraDistr::rhnorm(n = samps, sigma = 1)
}

initsBeta <- function(samps = 1){
  rnorm(n = samps, mean = 0, sd = 1)
}

# -------------------------------------------------------------------
# Set inits for JAGS
# -------------------------------------------------------------------
# inits <- list()
# for(i in 1:3){
#   inits[[i]] <- list(initsMu0(),initsSigma0Weak(),initsBeta())
#   
#   names(inits[[i]]) = c("mu.alpha","sigma.alpha", "beta")
#   
# }


B = matrix(nrow = data$n_siteGermination, ncol = 2)
B[,1] = .1
B[,2] = 1.5

inits = list(
  list(B = B, sigma = 50, mu.alpha = 0, mu.beta = 1.5, sigma.alpha = 10, sigma.beta = 10, rho = -.5),
  list(B = B*.5, sigma = 20, mu.alpha = -.2, mu.beta = .8, sigma.alpha = 50, sigma.beta = 50, rho = .5),
  list(B = B*.3, sigma = 10, mu.alpha = -.4, mu.beta = 1.8, sigma.alpha = 30, sigma.beta = 20, rho = 0))

# -------------------------------------------------------------------
# Call to JAGS
# -------------------------------------------------------------------
jm = jags.model("~/Dropbox/clarkiaSeedBanks/modelBuild/jags-germinationInteractions-randomInterceptSlope.R",
                data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

update(jm, n.iter = n.update)

parsToMonitor = c('mu.alpha', "mu.beta", "rho", "beta", "alpha",
                  "sigma.alpha", "sigma.beta")

zm.randomInterceptSlope = coda.samples(jm, variable.names = c(parsToMonitor),
                                  n.iter = n.iterations, thin = 1)

# -------------------------------------------------------------------
# Check convergence
# -------------------------------------------------------------------
MCMCvis::MCMCsummary(zm.randomInterceptSlope, c("mu.alpha","sigma.alpha","mu.beta","mu.alpha","rho"),n.eff=TRUE,round=3)

dev.off()
library(MCMCvis)
par(mfrow=c(2,1))
MCMCplot(zm.randomInterceptSlope, params = c("beta"),horiz=FALSE)
MCMCplot(zm.randomInterceptSlope, params = c("alpha"),horiz=FALSE)

alpha=MCMCchains(zm.randomInterceptSlope, params = c("alpha"))
beta=MCMCchains(zm.randomInterceptSlope, params = c("beta"))
n.iter=dim(alpha)[1]
index=sample(1:n.iter,1000)
par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
 plot(alpha[index,i],beta[index,i]) 
}

alpha = MCMCchains(zm.randomInterceptSlope,"alpha")
beta = MCMCchains(zm.randomInterceptSlope,"beta")

plot.listRanIntSlope = list()
for(j in 1:20){
  preds=alpha[,j]+outer(beta[,j],clim_pred_std)
  preds=apply(preds,2,boot::inv.logit)
  pred3 <- apply(preds,2,hdi,.95)
  pred4 <- apply(preds,2,median)
  pred.po.df <- cbind(clim_pred_std, data.frame(t(pred3)), median = pred4)
  
  plot.listRanIntSlope[[j]] = pred.po.df
}

siteNames = unique(climate$site)
par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  tmp = climate[climate$site==siteNames[i],]
  plot(tmp$p.std,tmp$seedlingJan/tmp$totalJan,xlim=c(min(clim_pred_std),max(clim_pred_std)),ylim=c(0,1))
  lines(plot.list.noPool[[i]]$clim_pred_std,plot.list.noPool[[i]]$median,type='l',lwd=2,lty='dotdash',col='green')
  lines(plot.listRanInt[[i]]$clim_pred_std,plot.listRanInt[[i]]$median,type='l',lwd=2,lty='dotted',col='orange')
  lines(plot.listRanIntSlope[[i]]$clim_pred_std,plot.listRanIntSlope[[i]]$median,type='l',lwd=2)
  text(.8*max(clim_pred_std),.95*1,siteNames[i])
}


alpha.ranIntSlope = (MCMCchains(zm.randomInterceptSlope,params="alpha"))
beta.ranIntSlope = (MCMCchains(zm.randomInterceptSlope,params="beta"))

b.ris=apply(beta.ranIntSlope,2,quantile,c(.025,.5,.975))
a.ris=apply(alpha.ranIntSlope,2,quantile,c(.025,.5,.975))

par(mfrow=c(1,1),mar=c(2,2,2,2))
plot(a.ris[2,],b.ris[2,],xlim=c(-4.5,-1),ylim=c(-3,2),pch=16)
segments(x0=a.ris[2,],y0=b.ris[1,],y1=b.ris[3,])
segments(y0=b.ris[2,],x0=a.ris[1,],x1=a.ris[3,])


position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)

siteNames = position$site 

par(mfrow=c(1,1))
a.ris=apply(boot::inv.logit(alpha.ranIntSlope),2,quantile,c(.025,.5,.975))
plot(position$easting,a.ris[2,],pch=19,ylim=c(0,.25))
segments(x0=position$easting,y0=a.ris[1,],y1=a.ris[3,])

par(mfrow=c(1,1))
b.ris=apply((beta.ranIntSlope),2,quantile,c(.025,.5,.975))
plot(position$easting,b.ris[2,],pch=19,ylim=c(-3,2))
segments(x0=position$easting,y0=b.ris[1,],y1=b.ris[3,])

par(mfrow=c(1,1),mar=c(2,2,2,2))
a.ris=apply((alpha.ranIntSlope),2,quantile,c(.025,.5,.975))
plot(a.ris[2,],b.ris[2,],xlim=c(-4.5,-1),ylim=c(-3,2),pch=16)
segments(x0=a.ris[2,],y0=b.ris[1,],y1=b.ris[3,])
segments(y0=b.ris[2,],x0=a.ris[1,],x1=a.ris[3,])

cor.vec = c()
for(i in 1:20){
  cor.vec[i]  =cor(alpha.ranIntSlope[,i],beta.ranIntSlope[,i])
  
}

cor.vec
par(mfrow=c(1,1),mar=c(2,2,2,2))
sa=read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE)
plot(sa$easting,cor.vec,ylim=c(-1,1))
text(sa$easting,cor.vec,siteNames)



position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)
library(lme4)
gm1 = glmer(cbind(seedlingJan,totalJan-seedlingJan) ~ (1 + t.std + p.std + t.std:p.std|site), data=climate, family = "binomial")
summary(gm1)


plot(position$easting,boot::inv.logit(fixef(gm1)+ranef(gm1)$site$`(Intercept)`))

gm1.ran=ranef(gm1)$site

# sensitivity to temperature reverses from positive to negative
plot(position$easting,gm1.ran$t.std)
# sensitivity to precipitation switches from negative to positive/0
plot(position$easting,gm1.ran$p.std)
# sensitivity to temperature:precipitation switches from positive to negative
plot(position$easting,gm1.ran$`t.std:p.std`)


groupMeans=climate %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(mu.t = mean(t),sd.t = sd(t), mu.p = mean(p),sd.p = sd(p))

climate=climate %>% 
  dplyr::left_join(groupMeans,by=c('site')) %>%
  dplyr::mutate(t.std.group = (t-mu.t)/sd.t, 
                p.std.group = (p-mu.p)/sd.p)


library(lme4)
gm1 = glmer(cbind(seedlingJan,totalJan-seedlingJan) ~ (1 + t.std.group + p.std.group + t.std.group:p.std.group|site), data=climate, family = "binomial")
summary(gm1)


plot(position$easting,boot::inv.logit(fixef(gm1)+ranef(gm1)$site$`(Intercept)`))

gm1.ran=ranef(gm1)$site

# sensitivity to temperature reverses from positive to negative
plot(position$easting,gm1.ran$t.std.group)
# sensitivity to precipitation switches from negative to positive/0
plot(position$easting,gm1.ran$p.std.group)
# sensitivity to temperature:precipitation switches from positive to negative
plot(position$easting,gm1.ran$`t.std.group:p.std.group`)

