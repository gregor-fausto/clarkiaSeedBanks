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
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)
# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(rjags) 
library(tidybayes)
library(tidyverse)
library(parallel)
library(lme4)
# -------------------------------------------------------------------
# Set directories
# -------------------------------------------------------------------
dataDirectory = "~/Dropbox/dataLibrary/postProcessingData-2021/"
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
# Add site centered climate variables
# -------------------------------------------------------------------

groupMeans=climate %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(mu.t = mean(t),sd.t = sd(t), mu.p = mean(p),sd.p = sd(p))

climate=climate %>% 
  dplyr::left_join(groupMeans,by=c('site')) %>%
  dplyr::mutate(t.std.group = (t-mu.t)/sd.t, 
                p.std.group = (p-mu.p)/sd.p)
# -------------------------------------------------------------------
# Position data
# -------------------------------------------------------------------

position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)


par(mfrow = c(1, 3),mar=c(3,3,3,3))

# Temperature and precipitation are positively correlated
plot(p ~ t, climate)
lines(lowess(climate$t, climate$p), col = "red")

# Temperature is uncorrelated with germination
plot(seedlingJan/totalJan ~ t, climate)
lines(lowess(climate$t, climate$seedlingJan/climate$totalJan), col = "red")

# Precipitation is negatively then positively correlated with germination
plot(seedlingJan/totalJan ~ p, climate)
lines(lowess(climate$p, climate$seedlingJan/climate$totalJan), col = "red")

# conditioning plots
# effect of temp is variable across levels of precip
par(mfrow=c(1,1))
library(lattice)
coplot(seedlingJan/totalJan ~ t | p, data = climate, 
       number = 4, rows = 1,
       panel = panel.smooth)

# conditioning plots
# effect of precip is positive at low temp
# but weakens at higher temps
par(mfrow=c(1,1))
coplot(seedlingJan/totalJan ~ p | t, data = climate, 
       number = 4, rows = 1,
       panel = panel.smooth)

g1 = ggplot(climate) +
  geom_point(aes(x= t, y = seedlingJan/totalJan) , alpha = 3/10, shape = 21, colour = "black",  fill = "brown", size = 3) +
  theme_minimal()

g2 =ggplot(climate) +
  geom_point(aes(x= p, y = seedlingJan/totalJan) , alpha = 3/10, shape = 21, colour = "black",  fill = "brown", size = 3) +
  theme_minimal()

gridExtra::grid.arrange(g1,g2)
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

# add climate data
data$t.std = climate$t.std
data$p.std = climate$p.std
data$t.mu.group = climate$mu.t
data$p.mu.group = climate$mu.p
data$t.std.group = climate$t.std.group
data$p.std.group = climate$p.std.group

# clim_pred_std = precip_pred_std
# muClim = muTemp
# sdClim = sdTemp
# gClim = g2

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
# Plot climate 
# -------------------------------------------------------------------

df.clim = climate %>%
  dplyr::select(site,year,t,p,t.std,p.std,t.std.group,p.std.group) %>%
  unique()

ggplot(df.clim) +
  geom_point(aes(x=t,y=p,color=as.factor(year)))

ggplot(df.clim) +
  geom_point(aes(x=t.std,y=p.std))

ggplot(df.clim) +
  geom_point(aes(x=year,y=p.std.group))

hist(df.clim$p);
hist(df.clim$p.std);
hist(df.clim$p.std.group)

# -------------------------------------------------------------------
# GLMMs with random intercept and slope: no centering
# -------------------------------------------------------------------

gm1 = glmer(cbind(seedlingJan,totalJan-seedlingJan) ~ 1 + t + (1 + t |site), data=climate, family = "binomial")
summary(gm1)
coef(gm1)
fixef(gm1)
ranef(gm1)

t_pred=seq(min(climate$t),max(climate$t),.1)
out=boot::inv.logit(fixef(gm1)[1]+fixef(gm1)[2]*t_pred)
plot(t_pred,out,type='l',ylim=c(0,1),xlab= "Temperature",ylab= "P(germination)")
points(climate$t,climate$seedlingJan/climate$totalJan,pch=1,col='gray80')
lines(t_pred,out,lwd=2,col='orange')

extract_pred=function(model,pred_seq){
  kappa = fixef(model)[1]
  eta = ranef(model)$site$`(Intercept)`
  
  beta = fixef(model)[2]
  nu = ranef(model)$site$t
  
  alpha.site=kappa+eta[i]
  beta.site=beta+nu[i]
  
  pred = alpha.site + beta.site*pred_seq
  return(pred)
}

pred.mat=matrix(NA,ncol=20,nrow=length(t_pred))
for(i in 1:20){
  pred.mat[,i] = extract_pred(gm1,t_pred)
}

out.mat.gm1=apply(pred.mat,2,boot::inv.logit)
plot(t_pred,out,type='l',ylim=c(0,1),xlab= "Temperature",ylab= "P(germination)")
points(climate$t,climate$seedlingJan/climate$totalJan,pch=1,col='gray80')
lines(t_pred,out,lwd=2,col='orange')
for(i in 1:20){
  lines(t_pred,out.mat.gm1[,i],lwd=1)
}

siteNames = unique(climate$site)
par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  plot(t_pred,out.mat.gm1[,i],type='l',ylim=c(0,1),xlab= "Temperature",ylab= "P(germination)")
  tmp = climate[climate$site==siteNames[i],]
  points(tmp$t,tmp$seedlingJan/tmp$totalJan,pch=1,col='gray80')
  lines(t_pred,out,lwd=2,col='orange')
  text(max(t_pred)*.9,.95,siteNames[i])
}

par(mfrow=c(1,1),mar=c(3,3,3,3))
plot(position$easting,boot::inv.logit(fixef(gm1)[1]+ranef(gm1)$site$`(Intercept)`))


# -------------------------------------------------------------------
# GLMMs with random intercept and slope: grand mean centering
# -------------------------------------------------------------------

gm2 = glmer(cbind(seedlingJan,totalJan-seedlingJan) ~ 1 + t.std + (1 + t.std |site), data=climate, family = "binomial")
summary(gm2)
coef(gm2)
fixef(gm2)
ranef(gm2)

t_pred=seq(min(climate$t),max(climate$t),.1)
t_pred=(t_pred-muTemp)/sdTemp
out=boot::inv.logit(fixef(gm2)[1]+fixef(gm2)[2]*t_pred)
plot(t_pred,out,type='l',ylim=c(0,1),xlab= "Temperature",ylab= "P(germination)")
points(climate$t.std,climate$seedlingJan/climate$totalJan,pch=1,col='gray80')
lines(t_pred,out,lwd=2,col='orange')

extract_pred=function(model,pred_seq){
  kappa = fixef(model)[1]
  eta = ranef(model)$site$`(Intercept)`
  
  beta = fixef(model)[2]
  nu = ranef(model)$site$t
  
  alpha.site=kappa+eta[i]
  beta.site=beta+nu[i]
  
  pred = alpha.site + beta.site*pred_seq
  return(pred)
}

pred.mat=matrix(NA,ncol=20,nrow=length(t_pred))
for(i in 1:20){
  pred.mat[,i] = extract_pred(gm2,t_pred)
}

out.mat.gm2=apply(pred.mat,2,boot::inv.logit)
plot(t_pred,out,type='l',ylim=c(0,1),xlab= "Temperature",ylab= "P(germination)")
points(climate$t.std,climate$seedlingJan/climate$totalJan,pch=1,col='gray80')
lines(t_pred,out,lwd=2,col='orange')
for(i in 1:20){
  lines(t_pred,out.mat.gm2[,i],lwd=1)
}

siteNames = unique(climate$site)
par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  plot(t_pred,out.mat.gm2[,i],type='l',ylim=c(0,1),xlab= "Temperature",ylab= "P(germination)")
  tmp = climate[climate$site==siteNames[i],]
  points(tmp$t.std,tmp$seedlingJan/tmp$totalJan,pch=1,col='gray80')
  lines(t_pred,out,lwd=2,col='orange')
  text(max(t_pred)*.9,.95,siteNames[i])
}

par(mfrow=c(1,1),mar=c(3,3,3,3))
plot(position$easting,boot::inv.logit(fixef(gm2)[1]+ranef(gm2)$site$`(Intercept)`))


# -------------------------------------------------------------------
# GLMMs with random intercept and slope: group mean centering
# -------------------------------------------------------------------

gm3 = glmer(cbind(seedlingJan,totalJan-seedlingJan) ~ 1 + t.std.group + (1 + t.std.group |site), data=climate, family = "binomial")
summary(gm3)
coef(gm3)
fixef(gm3)
ranef(gm3)

# t_pred=seq(min(climate$t),max(climate$t),.1)
# t_pred=(t_pred-muTemp)/sdTemp
# out=boot::inv.logit(fixef(gm3)[1]+fixef(gm3)[2]*t_pred)
# plot(t_pred,out,type='l',ylim=c(0,1),xlab= "Temperature",ylab= "P(germination)")
# points(climate$t.std.group,climate$seedlingJan/climate$totalJan,pch=1,col='gray80')
# lines(t_pred,out,lwd=2,col='orange')

extract_pred=function(model,pred_seq){
  kappa = fixef(model)[1]
  eta = ranef(model)$site$`(Intercept)`
  
  beta = fixef(model)[2]
  nu = ranef(model)$site$t
  
  alpha.site=kappa+eta[i]
  beta.site=beta+nu[i]
  
  pred = alpha.site + beta.site*pred_seq
  return(pred)
}

pred.mat=matrix(NA,ncol=20,nrow=length(t_pred))
for(i in 1:20){
  t_pred=seq(min(climate$t),max(climate$t),.1)
  t_pred=(t_pred-groupMeans$mu.t[i])/groupMeans$sd.t[i]
  pred.mat[,i] = extract_pred(gm3,t_pred)
}

out.mat.gm3=apply(pred.mat,2,boot::inv.logit)
# plot(t_pred,out,type='l',ylim=c(0,1),xlab= "Temperature",ylab= "P(germination)")
# points(climate$t.std.group,climate$seedlingJan/climate$totalJan,pch=1,col='gray80')
# lines(t_pred,out,lwd=2,col='orange')
# for(i in 1:20){
#   lines(t_pred,out.mat[,i],lwd=1)
# }

siteNames = unique(climate$site)
par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  t_pred=seq(min(climate$t),max(climate$t),.1)
  t_pred=(t_pred-groupMeans$mu.t[i])/groupMeans$sd.t[i]
  plot(t_pred,out.mat.gm3[,i],type='l',ylim=c(0,1),xlab= "Temperature",ylab= "P(germination)")
  tmp = climate[climate$site==siteNames[i],]
  points(tmp$t.std.group,tmp$seedlingJan/tmp$totalJan,pch=1,col='gray80')
  lines(t_pred,out,lwd=2,col='orange')
  text(max(t_pred)*.9,.95,siteNames[i])
}

# -------------------------------------------------------------------
# Compare estimates
# -------------------------------------------------------------------
fixef(gm1)[1] + ranef(gm1)$site[,1]
fixef(gm2)[1] + ranef(gm2)$site[,1]

siteNames = unique(climate$site)
par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  t_pred=seq(min(climate$t),max(climate$t),.1)
  plot(t_pred,out.mat.gm1[,i],type='l',ylim=c(0,1),xlim=c(min(t_pred)-1,max(t_pred)+1),
       xlab= "Temperature",ylab= "P(germination)")
  lines(t_pred,out.mat.gm2[,i])
  t_pred=seq(min(climate$t),max(climate$t),.1)
  lines(t_pred,out.mat.gm3[,i])
  tmp = climate[climate$site==siteNames[i],]
  points(tmp$t,tmp$seedlingJan/tmp$totalJan,pch=1,col='gray80')
  lines(t_pred,out,lwd=2,col='orange')
  text(max(t_pred)*.9,.95,siteNames[i])
  text(unique(tmp$t),c(.4,.5,.6),round(unique(tmp$p),0))
}

siteNames = unique(climate$site)
par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  tmp.clim = climate %>% dplyr::select(site,year,t,p) %>% unique
  plot(tmp.clim$t,tmp.clim$p,type='p', xlim=c(5,11), ylim=c(50,200),
       xlab= "Temperature",ylab= "P(germination)",pch=16,col='gray90')
  tmp = tmp.clim[tmp.clim$site==siteNames[i],]
  points(tmp$t,tmp$p,col='orange',pch=16)
  text(max(tmp.clim$t)*.9,max(tmp.clim$p)*.9,siteNames[i])
}

siteNames = unique(climate$site)
par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  tmp.clim = climate %>% dplyr::select(site,year,t,p) %>% unique
  plot(tmp.clim$t,tmp.clim$p,type='p', #xlim=c(5,11), ylim=c(50,200),
       xlab= "Temperature",ylab= "P(germination)",pch=16,col='gray90')
  abline(v= muTemp, h=muPrecip, lty = 'dotted')
  tmp = tmp.clim[tmp.clim$site==siteNames[i],]
  points(tmp$t,tmp$p,col='orange',pch=16)
  text(max(tmp.clim$t)*.9,max(tmp.clim$p)*.9,siteNames[i])
}

siteNames = unique(climate$site)
par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  tmp.clim = climate %>% dplyr::select(site,year,t.std,p.std) %>% unique
  plot(tmp.clim$t.std,tmp.clim$p.std,type='p', 
       xlab= "Temperature",ylab= "P(germination)",pch=16,col='gray90')
  abline(v= 0, h=0, lty = 'dotted')
  tmp = tmp.clim[tmp.clim$site==siteNames[i],]
  points(tmp$t.std,tmp$p.std,col='orange',pch=16)
  text(max(tmp.clim$t.std)*.9,max(tmp.clim$p.std)*.9,siteNames[i])
}


siteNames = unique(climate$site)
par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  tmp.clim = climate %>% dplyr::select(site,year,t.std.group,p.std.group) %>% unique
  plot(tmp.clim$t.std.group,tmp.clim$p.std.group,type='p', #xlim=c(-2,3),ylim=c(-4,4),
       xlab= "Temperature",ylab= "P(germination)",pch=16,col='gray90')
  abline(v= 0, h=0, lty = 'dotted')
  tmp = tmp.clim[tmp.clim$site==siteNames[i],]
  points(tmp$t.std.group,tmp$p.std.group,col='orange',pch=16)
  text(max(tmp.clim$t.std.group)*.9,max(tmp.clim$p.std.group)*.9,siteNames[i])
}


# -------------------------------------------------------------------
# Random intercept+random slope model in JAGS
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
jm = jags.model("~/Dropbox/clarkiaSeedBanks/modelBuild/jags-germinationInteractions-randomInterceptSlope.R",
                data = data,  n.adapt = n.adapt)

update(jm, n.iter = n.update)

parsToMonitor = c('mu.alpha', "mu.beta", "rho", "beta", "alpha",
                  "sigma.alpha", "sigma.beta")

zm.randomInterceptSlope = coda.samples(jm, variable.names = c(parsToMonitor),
                                       n.iter = n.iterations, thin = 1)

# -------------------------------------------------------------------
# Check convergence
# -------------------------------------------------------------------
MCMCvis::MCMCsummary(zm.randomInterceptSlope, c("mu.alpha","sigma.alpha","mu.beta","sigma.beta","rho"),n.eff=TRUE,round=3)

# compare to gm2
summary(gm2)
library(MCMCvis)
dev.off()
par(mfrow=c(2,2))
hist(MCMCchains(zm.randomInterceptSlope,"mu.alpha"),breaks=50,col='gray80',border=0,
     main="Intercept, mean",xlab='')
abline(v=fixef(gm2)[1],lwd=2,col='orange')

hist(MCMCchains(zm.randomInterceptSlope,"mu.beta"),breaks=50,col='gray80',border=0,
     main="Slope, mean",xlab='')
abline(v=fixef(gm2)[2],lwd=2,col='orange')

hist(MCMCchains(zm.randomInterceptSlope,"sigma.alpha"),breaks=50,col='gray80',border=0,
     main="Intercept, standard deviation",xlab='')
alpha.sd=attr(summary(gm2)$varcor$site, "stddev")[1]
abline(v=alpha.sd,lwd=2,col='orange')

hist(MCMCchains(zm.randomInterceptSlope,"sigma.beta"),breaks=50,col='gray80',border=0,
     main="Slope, standard deviation",xlab='')
beta.sd=attr(summary(gm2)$varcor$site, "stddev")[2]
abline(v=beta.sd,lwd=2,col='orange')

par(mfrow=c(1,1))
hist(MCMCchains(zm.randomInterceptSlope,"rho"),breaks=50,col='gray80',border=0,
     main="Correlation",xlab='')
cor.t=attr(summary(gm2)$varcor$site, "correlation")[1,2]
abline(v=cor.t,lwd=2,col='orange')

# fixed effects correlation are less interpretable
par(mfrow=c(1,1))
cor(MCMCchains(zm.randomInterceptSlope,"mu.alpha"),MCMCchains(zm.randomInterceptSlope,"mu.beta"))
summary(gm2)$vcov@factors$correlation
abline(v=cor.t,lwd=2,col='orange')


dev.off()
library(MCMCvis)
par(mfrow=c(2,1))
beta = MCMCchains(zm.randomInterceptSlope,"beta")
alpha = MCMCchains(zm.randomInterceptSlope,"alpha")

lm.alpha=ranef(gm2)$site[,1]+fixef(gm2)[1]

par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  hist(alpha[,i],breaks=25,col='gray80',border=0,
       main="",xlab='')
  abline(v=lm.alpha[i],lwd=2,col='orange')
}

lm.beta=ranef(gm2)$site[,2]+fixef(gm2)[2]

par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  hist(beta[,i],breaks=25,col='gray80',border=0,
       main="",xlab='')
  abline(v=lm.beta[i],lwd=2,col='orange')
}


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

temp_pred=seq(min(climate$t),max(climate$t),.1)
temp_pred.std = (temp_pred-muTemp)/sdTemp

plot.listRanIntSlope = list()
for(j in 1:20){
  preds=alpha[,j]+outer(beta[,j],temp_pred.std)
  preds=apply(preds,2,boot::inv.logit)
  pred3 <- data.frame(t(apply(preds,2,hdi,.95)))
  colnames(pred3) = c('lower','upper')
  pred4 <- apply(preds,2,median)
  pred.po.df <- cbind(temp_pred.std, pred3, median = pred4)
  plot.listRanIntSlope[[j]] = pred.po.df
}

siteNames = unique(climate$site)
par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  tmp = climate[climate$site==siteNames[i],]
  plot(tmp$t.std,tmp$seedlingJan/tmp$totalJan,xlim=c(min(temp_pred.std),max(temp_pred.std)),ylim=c(0,1))
  polygon(c(plot.listRanIntSlope[[i]]$temp_pred.std,
            rev(plot.listRanIntSlope[[i]]$temp_pred.std)),
          c(plot.listRanIntSlope[[i]]$lower,
            rev(plot.listRanIntSlope[[i]]$upper)),col='gray90',border=0)
  lines(plot.listRanIntSlope[[i]]$temp_pred.std,plot.listRanIntSlope[[i]]$median,type='l',lwd=2)
  text(.8*max(clim_pred_std),.95*1,siteNames[i])
}

n.iter = dim(alpha)[1]
plot.listRanIntSlope = list()
for(i in 1:20){
  index=sample(1:n.iter,100)
  mat = matrix(NA,nrow=100,ncol=length(temp_pred.std))
  for(j in 1:100){
    preds=alpha[index[j],i]+outer(beta[index[j],i],temp_pred.std)
    preds=apply(preds,2,boot::inv.logit)
    mat[j,] = preds
  }
  plot.listRanIntSlope[[i]] = mat
}

siteNames = unique(climate$site)
par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  tmp = climate[climate$site==siteNames[i],]
  plot(tmp$t.std,tmp$seedlingJan/tmp$totalJan,xlim=c(min(temp_pred.std),max(temp_pred.std)),ylim=c(0,.5))
  for(j in 1:100){
    lines(temp_pred.std,plot.listRanIntSlope[[i]][j,],type='l',lwd=.5,col=sample(1:10,1))
  }
  text(.8*max(clim_pred_std),.95*1,siteNames[i])
}

position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)

siteNames = position$site 

par(mfrow=c(1,1),mar=c(3,3,3,3))
a.ris=apply(boot::inv.logit(alpha),2,quantile,c(.025,.5,.975))
plot(position$easting,a.ris[2,],pch=19,ylim=c(0,.25))
segments(x0=position$easting,y0=a.ris[1,],y1=a.ris[3,])

par(mfrow=c(1,1))
b.ris=apply((beta),2,quantile,c(.025,.5,.975))
plot(position$easting,b.ris[2,],pch=19,ylim=c(-3,2))
segments(x0=position$easting,y0=b.ris[1,],y1=b.ris[3,])


# -------------------------------------------------------------------
# Random intercept+random slope model in JAGS (Wishart)
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Set inits for JAGS
# -------------------------------------------------------------------
# Inits function
inits <- function() {
  list(mu.alpha=rnorm(1),mu.beta=rnorm(1), Tau_B_raw=bayesm::rwishart(3, diag(2)*var_vec)$W, sigma_res=rlnorm(1), xi_prior=rlnorm(3), mu_raw_prior=rnorm(3), Tau_B_raw_prior=rwishart(4, diag(3)*var_vec)$W)
}
# -------------------------------------------------------------------
# Call to JAGS
# -------------------------------------------------------------------
# http://mmeredith.net/blog/2020/Correlated_priors.htm
# using a scaled wishart prior
data$s = c(3,3)
#data$df = 2
load.module("glm")

jm = jags.model("~/Dropbox/clarkiaSeedBanks/modelBuild/jags-germinationInteractions-randomInterceptSlopeWishart.R",
                data = data,  n.adapt = n.adapt)

update(jm, n.iter = n.update)

parsToMonitor = c('mu.alpha', "mu.beta", "rho.a.b1", "beta", "alpha",
                  "sigma.B")

zm.randomInterceptSlopeWishart = coda.samples(jm, variable.names = c(parsToMonitor),
                                       n.iter = n.iterations, thin = 1)

MCMCvis::MCMCsummary(zm.randomInterceptSlopeWishart,params=c("mu.alpha","mu.beta","rho.a.b1","sigma.B"))
MCMCvis::MCMCsummary(zm.randomInterceptSlope, c("mu.alpha","sigma.alpha","mu.beta","sigma.beta","rho"),n.eff=TRUE,round=3)


# compare to gm2
summary(gm2)
library(MCMCvis)
dev.off()
par(mfrow=c(2,2))
hist(MCMCchains(zm.randomInterceptSlopeWishart,"mu.alpha"),breaks=50,col='gray80',border=0,
     main="Intercept, mean",xlab='')
abline(v=fixef(gm2)[1],lwd=2,col='orange')

hist(MCMCchains(zm.randomInterceptSlopeWishart,"mu.beta"),breaks=50,col='gray80',border=0,
     main="Slope, mean",xlab='')
abline(v=fixef(gm2)[2],lwd=2,col='orange')

hist(MCMCchains(zm.randomInterceptSlopeWishart,"sigma.alpha"),breaks=50,col='gray80',border=0,
     main="Intercept, standard deviation",xlab='')
alpha.sd=attr(summary(gm2)$varcor$site, "stddev")[1]
abline(v=alpha.sd,lwd=2,col='orange')

hist(MCMCchains(zm.randomInterceptSlope,"sigma.beta"),breaks=50,col='gray80',border=0,
     main="Slope, standard deviation",xlab='')
beta.sd=attr(summary(gm2)$varcor$site, "stddev")[2]
abline(v=beta.sd,lwd=2,col='orange')

par(mfrow=c(1,1))
hist(MCMCchains(zm.randomInterceptSlope,"rho"),breaks=50,col='gray80',border=0,
     main="Correlation",xlab='')
cor.t=attr(summary(gm2)$varcor$site, "correlation")[1,2]
abline(v=cor.t,lwd=2,col='orange')


dev.off()
library(MCMCvis)
par(mfrow=c(2,1))
beta = MCMCchains(zm.randomInterceptSlopeWishart,"beta")
alpha = MCMCchains(zm.randomInterceptSlopeWishart,"alpha")

lm.alpha=ranef(gm2)$site[,1]+fixef(gm2)[1]

par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  hist(alpha[,i],breaks=25,col='gray80',border=0,
       main="",xlab='')
  abline(v=lm.alpha[i],lwd=2,col='orange')
}

lm.beta=ranef(gm2)$site[,2]+fixef(gm2)[2]

par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  hist(beta[,i],breaks=25,col='gray80',border=0,
       main="",xlab='')
  abline(v=lm.beta[i],lwd=2,col='orange')
}


alpha=MCMCchains(zm.randomInterceptSlopeWishart, params = c("alpha"))
beta=MCMCchains(zm.randomInterceptSlopeWishart, params = c("beta"))
n.iter=dim(alpha)[1]
index=sample(1:n.iter,1000)
par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  plot(alpha[index,i],beta[index,i]) 
}

alpha = MCMCchains(zm.randomInterceptSlopeWishart,"alpha")
beta = MCMCchains(zm.randomInterceptSlopeWishart,"beta")

temp_pred=seq(min(climate$t),max(climate$t),.1)
temp_pred.std = (temp_pred-muTemp)/sdTemp

zm.randomInterceptSlopeWishart = list()
for(j in 1:20){
  preds=alpha[,j]+outer(beta[,j],temp_pred.std)
  preds=apply(preds,2,boot::inv.logit)
  pred3 <- data.frame(t(apply(preds,2,hdi,.95)))
  colnames(pred3) = c('lower','upper')
  pred4 <- apply(preds,2,median)
  pred.po.df <- cbind(temp_pred.std, pred3, median = pred4)
  zm.randomInterceptSlopeWishart[[j]] = pred.po.df
}

siteNames = unique(climate$site)
par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  tmp = climate[climate$site==siteNames[i],]
  plot(tmp$t.std,tmp$seedlingJan/tmp$totalJan,xlim=c(min(temp_pred.std),max(temp_pred.std)),ylim=c(0,1))
  polygon(c(zm.randomInterceptSlopeWishart[[i]]$temp_pred.std,
            rev(zm.randomInterceptSlopeWishart[[i]]$temp_pred.std)),
          c(zm.randomInterceptSlopeWishart[[i]]$lower,
            rev(zm.randomInterceptSlopeWishart[[i]]$upper)),col='gray90',border=0)
  lines(zm.randomInterceptSlopeWishart[[i]]$temp_pred.std,zm.randomInterceptSlopeWishart[[i]]$median,type='l',lwd=2)
  text(.8*max(temp_pred.std),.95*1,siteNames[i])
}

n.iter = dim(alpha)[1]
plot.listRanIntSlope = list()
for(i in 1:20){
  index=sample(1:n.iter,100)
  mat = matrix(NA,nrow=100,ncol=length(temp_pred.std))
  for(j in 1:100){
    preds=alpha[index[j],i]+outer(beta[index[j],i],temp_pred.std)
    preds=apply(preds,2,boot::inv.logit)
    mat[j,] = preds
  }
  plot.listRanIntSlope[[i]] = mat
}

siteNames = unique(climate$site)
par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  tmp = climate[climate$site==siteNames[i],]
  plot(tmp$t.std,tmp$seedlingJan/tmp$totalJan,xlim=c(min(temp_pred.std),max(temp_pred.std)),ylim=c(0,.5))
  for(j in 1:100){
    lines(temp_pred.std,plot.listRanIntSlope[[i]][j,],type='l',lwd=.5,col=sample(1:10,1))
  }
  text(.8*max(temp_pred.std),.95*1,siteNames[i])
}

position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)

siteNames = position$site 

par(mfrow=c(1,1),mar=c(3,3,3,3))
a.ris=apply(boot::inv.logit(alpha),2,quantile,c(.025,.5,.975))
plot(position$easting,a.ris[2,],pch=19,ylim=c(0,.25))
segments(x0=position$easting,y0=a.ris[1,],y1=a.ris[3,])

par(mfrow=c(1,1))
b.ris=apply((beta),2,quantile,c(.025,.5,.975))
plot(position$easting,b.ris[2,],pch=19,ylim=c(-3,2))
segments(x0=position$easting,y0=b.ris[1,],y1=b.ris[3,])

# -------------------------------------------------------------------
# Random intercept+random slope model with 2 slopes in JAGS (Wishart)
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Call to JAGS
# -------------------------------------------------------------------
# http://mmeredith.net/blog/2020/Correlated_priors.htm
# using a scaled wishart prior
data$s = c(3,3,3)
#data$df = 2
load.module("glm")

jm = jags.model("~/Dropbox/clarkiaSeedBanks/modelBuild/jags-germinationInteractions-randomInterceptSlopeWishart2.R",
                data = data,  n.adapt = n.adapt)

update(jm, n.iter = n.update)

parsToMonitor = c('mu.alpha', "mu.beta_1", "mu.beta_2", 
                  "alpha","beta_1", "beta_2",
                  "sigma.B",
                  "rho.a.b1", "rho.a.b2", "rho.b1.b2")

zm.randomInterceptSlopeWishart2 = coda.samples(jm, variable.names = c(parsToMonitor),
                                              n.iter = n.iterations, thin = 1)

MCMCvis::MCMCsummary(zm.randomInterceptSlopeWishart2,params=c("mu.alpha","mu.beta_1","mu.beta_2"))
MCMCvis::MCMCsummary(zm.randomInterceptSlopeWishart2,params=c("sigma.B"))
MCMCvis::MCMCsummary(zm.randomInterceptSlopeWishart2,params=c("rho.a.b1","rho.a.b2","rho.b1.b2"))

gm4 = glmer(cbind(seedlingJan,totalJan-seedlingJan) ~ 1 + t.std + p.std + (1 + t.std + p.std |site), data=climate, family = "binomial")

MCMCvis::MCMCsummary(zm.randomInterceptSlopeWishart2,params=c("mu.alpha","mu.beta_1","mu.beta_2"))
par(mfrow=c(1,3))
hist(MCMCvis::MCMCchains(zm.randomInterceptSlopeWishart2,"mu.alpha"),col='gray90',border="white",breaks=20)
abline(v=fixef(gm4)[1],col='orange',lwd=2)
hist(MCMCvis::MCMCchains(zm.randomInterceptSlopeWishart2,"mu.beta_1"),col='gray90',border="white",breaks=20)
abline(v=fixef(gm4)[2],col='orange',lwd=2)
hist(MCMCvis::MCMCchains(zm.randomInterceptSlopeWishart2,"mu.beta_2"),col='gray90',border="white",breaks=20)
abline(v=fixef(gm4)[3],col='orange',lwd=2)

alpha = MCMCvis::MCMCchains(zm.randomInterceptSlopeWishart2,"alpha")
alpha_lme4=fixef(gm4)[1]+ranef(gm4)$site[,1]
par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  hist(alpha[,i],col='gray90',border="white",breaks=20,main='')
  abline(v=alpha_lme4[i],col='orange',lwd=2)
}



# -------------------------------------------------------------------
# Random intercept+random slope model with 3 slopes in JAGS (Wishart)
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Call to JAGS
# -------------------------------------------------------------------
# http://mmeredith.net/blog/2020/Correlated_priors.htm
# using a scaled wishart prior
data$s = c(3,3,3,3)
#data$df = 2
load.module("glm")

jm = jags.model("~/Dropbox/clarkiaSeedBanks/modelBuild/jags-germinationInteractions-randomInterceptSlopeWishart3.R",
                data = data,  n.adapt = n.adapt)

update(jm, n.iter = n.update)

parsToMonitor = c('mu.alpha', "mu.beta_1", "mu.beta_2", "mu.beta_3",
                  "alpha","beta_1", "beta_2", "beta_3",
                  "sigma.B",
                  "rho.a.b1", "rho.a.b2", "rho.a.b3",
                              "rho.b1.b2","rho.b1.b3",
                                          "rho.b2.b3")

zm.randomInterceptSlopeWishart3 = coda.samples(jm, variable.names = c(parsToMonitor),
                                               n.iter = n.iterations, thin = 1)

MCMCvis::MCMCsummary(zm.randomInterceptSlopeWishart3,params=c("mu.alpha","mu.beta_1","mu.beta_2","mu.beta_3"))
MCMCvis::MCMCsummary(zm.randomInterceptSlopeWishart3,params=c("sigma.B"))
MCMCvis::MCMCsummary(zm.randomInterceptSlopeWishart3,params=c("rho.a.b1", "rho.a.b2", "rho.a.b3",
                                                              "rho.b1.b2","rho.b1.b3",
                                                              "rho.b2.b3"))

gm5 = glmer(cbind(seedlingJan,totalJan-seedlingJan) ~ 1 + t.std + p.std + t.std*p.std + (1 + t.std + p.std + t.std*p.std |site), data=climate, family = "binomial")

MCMCvis::MCMCsummary(zm.randomInterceptSlopeWishart3,params=c("mu.alpha","mu.beta_1","mu.beta_2","mu.beta_3"))
par(mfrow=c(1,4))
hist(MCMCvis::MCMCchains(zm.randomInterceptSlopeWishart3,"mu.alpha"),col='gray90',border="white",breaks=20)
abline(v=fixef(gm5)[1],col='orange',lwd=2)
hist(MCMCvis::MCMCchains(zm.randomInterceptSlopeWishart3,"mu.beta_1"),col='gray90',border="white",breaks=20)
abline(v=fixef(gm5)[2],col='orange',lwd=2)
hist(MCMCvis::MCMCchains(zm.randomInterceptSlopeWishart3,"mu.beta_2"),col='gray90',border="white",breaks=20)
abline(v=fixef(gm5)[3],col='orange',lwd=2)
hist(MCMCvis::MCMCchains(zm.randomInterceptSlopeWishart3,"mu.beta_3"),col='gray90',border="white",breaks=20)
abline(v=fixef(gm5)[4],col='orange',lwd=2)

alpha = MCMCvis::MCMCchains(zm.randomInterceptSlopeWishart3,"alpha")
alpha_lme4=fixef(gm5)[1]+ranef(gm5)$site[,1]
par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  hist(alpha[,i],col='gray90',border="white",breaks=20,main='')
  abline(v=alpha_lme4[i],col='orange',lwd=2)
}

mu.alpha = as.vector( MCMCvis::MCMCchains(zm.randomInterceptSlopeWishart3,"mu.alpha"))
mu.beta_1 = as.vector( MCMCvis::MCMCchains(zm.randomInterceptSlopeWishart3,"mu.beta_1"))
mu.beta_2 = as.vector(MCMCvis::MCMCchains(zm.randomInterceptSlopeWishart3,"mu.beta_2"))
mu.beta_3 = as.vector(MCMCvis::MCMCchains(zm.randomInterceptSlopeWishart3,"mu.beta_3"))

# varying temperature at average site
temp_pred=seq(min(climate$t),max(climate$t),.1)
temp_pred.std = (temp_pred-muTemp)/sdTemp

p.minus1 = -1
p.zero = 0
p.plus1 = 1
precip_category = c(-1,0,1)

preds = matrix(NA,nrow=dim(alpha)[1],ncol=length(temp_pred.std))
plot.list = list()
for(i in 1:3){
 tmp <- as.vector(mu.alpha)+outer(mu.beta_1,temp_pred.std) + as.vector(mu.beta_2*precip_category[i]) + outer(mu.beta_3,temp_pred.std)*precip_category[i]
 plot.list[[i]] = apply( boot::inv.logit(tmp),2,quantile,c(.025,.5,.975))
 }

par(mfrow=c(1,3),mar=c(3,3,3,3))
for(i in 1:3){
plot(temp_pred.std,rep(NA,length(temp_pred.std)),type='l',ylim=c(0,.75),
     main = paste0("Precipitation: ",precip_category[i]," SD"),
     xlab = "Temp")
polygon(c(temp_pred.std,rev(temp_pred.std)),
        c(plot.list[[i]][1,],rev(plot.list[[i]][3,])),
          col='gray90')
lines(temp_pred.std,plot.list[[i]][2,],lwd=2,col='orange')
}

# varying precipitation at average site
precip_pred=seq(min(climate$p),max(climate$p),4)
precip_pred.std = (precip_pred-muPrecip)/sdPrecip

temp_category = c(-1,0,1)

preds = matrix(NA,nrow=dim(alpha)[1],ncol=length(precip_pred.std))
plot.list = list()
for(i in 1:3){
  tmp <- as.vector(mu.alpha)+ as.vector(mu.beta_1*temp_category[i]) + outer(mu.beta_2,precip_pred.std) + outer(mu.beta_3,precip_pred.std)*temp_category[i]
  plot.list[[i]] = apply( boot::inv.logit(tmp),2,quantile,c(.025,.5,.975))
}

par(mfrow=c(1,3),mar=c(3,3,3,3))
for(i in 1:3){
  plot(precip_pred.std,rep(NA,length(precip_pred.std)),type='l',ylim=c(0,1),
       main = paste0("Temperature: ",temp_category[i]," SD"),
       xlab="Precip")
  polygon(c(precip_pred.std,rev(precip_pred.std)),
          c(plot.list[[i]][1,],rev(plot.list[[i]][3,])),
          col='gray90')
  lines(precip_pred.std,plot.list[[i]][2,],lwd=2,col='orange')
}

# site specific responses
siteNames = unique(climate$site)

alpha = ( MCMCvis::MCMCchains(zm.randomInterceptSlopeWishart3,"alpha"))
beta_1 = ( MCMCvis::MCMCchains(zm.randomInterceptSlopeWishart3,"beta_1"))
beta_2 = (MCMCvis::MCMCchains(zm.randomInterceptSlopeWishart3,"beta_2"))
beta_3 = (MCMCvis::MCMCchains(zm.randomInterceptSlopeWishart3,"beta_3"))

precip_pred=seq(min(climate$p),max(climate$p),4)
precip_pred.std = (precip_pred-muPrecip)/sdPrecip
temp_category = c(-1,0,1)

preds = matrix(NA,nrow=dim(alpha)[1],ncol=length(precip_pred.std))
plot.list.min1 = list()
plot.list.zero = list()
plot.list.plus1 = list()
for(i in 1:20){
  tmp = climate[climate$site==siteNames[i],]
  #temp_category = mean(unique(tmp$t.std))
  tmp <- as.vector(alpha[,i])+ as.vector(beta_1[,i]*temp_category[1]) + outer(beta_2[,i],precip_pred.std) + outer(beta_3[,i],precip_pred.std)*temp_category[1]
  plot.list.min1[[i]] = apply( boot::inv.logit(tmp),2,quantile,c(.025,.5,.975))
  
  tmp <- as.vector(alpha[,i])+ as.vector(beta_1[,i]*temp_category[2]) + outer(beta_2[,i],precip_pred.std) + outer(beta_3[,i],precip_pred.std)*temp_category[2]
  plot.list.zero[[i]] = apply( boot::inv.logit(tmp),2,quantile,c(.025,.5,.975))
  
  tmp <- as.vector(alpha[,i])+ as.vector(beta_1[,i]*temp_category[3]) + outer(beta_2[,i],precip_pred.std) + outer(beta_3[,i],precip_pred.std)*temp_category[3]
  plot.list.plus1[[i]] = apply( boot::inv.logit(tmp),2,quantile,c(.025,.5,.975))
}

par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  plot(precip_pred.std,rep(NA,length(precip_pred.std)),type='l',ylim=c(0,1),
     #  main = paste0("Temperature: ",temp_category[i]," SD"),
       xlab="Precip")
  preds1 = plot.list.min1[[i]]
  preds2 = plot.list.zero[[i]]
  preds3 = plot.list.plus1[[i]]
  
  # polygon(c(precip_pred.std,rev(precip_pred.std)),
  #         c(preds1[1,],rev(preds1[3,])),
  #         col='gray90')
  lines(precip_pred.std,preds1[2,],lwd=2,col='lightblue')
  
  # polygon(c(precip_pred.std,rev(precip_pred.std)),
  #         c(preds2[1,],rev(preds2[3,])),
  #         col='gray90')
  lines(precip_pred.std,preds2[2,],lwd=2,col='orange')
  
  # polygon(c(precip_pred.std,rev(precip_pred.std)),
  #         c(preds3[1,],rev(preds3[3,])),
  #         col='gray90')
  lines(precip_pred.std,preds3[2,],lwd=2,col='red')
  
  tmp = climate[climate$site==siteNames[i],]
  points(tmp$p.std,tmp$seedlingJan/tmp$totalJan,pch=1,col='black')
  text(2,.9,siteNames[i])
}



par(mfrow=c(1,3),mar=c(3,3,3,3))
plot(precip_pred.std,rep(NA,length(precip_pred.std)),type='l',ylim=c(0,1),
     #  main = paste0("Temperature: ",temp_category[i]," SD"),
     xlab="Precip")
for(i in 1:20){
  preds1 = plot.list.min1[[i]]
  preds2 = plot.list.zero[[i]]
  preds3 = plot.list.plus1[[i]]
  
  # polygon(c(precip_pred.std,rev(precip_pred.std)),
  #         c(preds1[1,],rev(preds1[3,])),
  #         col='gray90')
  lines(precip_pred.std,preds1[2,],lwd=2,col='lightblue')
  
  # polygon(c(precip_pred.std,rev(precip_pred.std)),
  #         c(preds2[1,],rev(preds2[3,])),
  #         col='gray90')
  #lines(precip_pred.std,preds2[2,],lwd=2,col='orange')
  
  # polygon(c(precip_pred.std,rev(precip_pred.std)),
  #         c(preds3[1,],rev(preds3[3,])),
  #         col='gray90')
 # lines(precip_pred.std,preds3[2,],lwd=2,col='red')
  
  # tmp = climate[climate$site==siteNames[i],]
  # points(tmp$p.std,tmp$seedlingJan/tmp$totalJan,pch=1,col='black')
  # text(2,.9,siteNames[i])
}

plot(precip_pred.std,rep(NA,length(precip_pred.std)),type='l',ylim=c(0,1),
     #  main = paste0("Temperature: ",temp_category[i]," SD"),
     xlab="Precip")
for(i in 1:20){
  preds1 = plot.list.min1[[i]]
  preds2 = plot.list.zero[[i]]
  preds3 = plot.list.plus1[[i]]
  
  # polygon(c(precip_pred.std,rev(precip_pred.std)),
  #         c(preds1[1,],rev(preds1[3,])),
  #         col='gray90')
  #lines(precip_pred.std,preds1[2,],lwd=2,col='lightblue')
  
  # polygon(c(precip_pred.std,rev(precip_pred.std)),
  #         c(preds2[1,],rev(preds2[3,])),
  #         col='gray90')
  lines(precip_pred.std,preds2[2,],lwd=2,col='orange')
  
  # polygon(c(precip_pred.std,rev(precip_pred.std)),
  #         c(preds3[1,],rev(preds3[3,])),
  #         col='gray90')
  # lines(precip_pred.std,preds3[2,],lwd=2,col='red')
  
  # tmp = climate[climate$site==siteNames[i],]
  # points(tmp$p.std,tmp$seedlingJan/tmp$totalJan,pch=1,col='black')
  # text(2,.9,siteNames[i])
}

plot(precip_pred.std,rep(NA,length(precip_pred.std)),type='l',ylim=c(0,1),
     #  main = paste0("Temperature: ",temp_category[i]," SD"),
     xlab="Precip")
for(i in 1:20){
  preds1 = plot.list.min1[[i]]
  preds2 = plot.list.zero[[i]]
  preds3 = plot.list.plus1[[i]]
  
  # polygon(c(precip_pred.std,rev(precip_pred.std)),
  #         c(preds1[1,],rev(preds1[3,])),
  #         col='gray90')
  #lines(precip_pred.std,preds1[2,],lwd=2,col='lightblue')
  
  # polygon(c(precip_pred.std,rev(precip_pred.std)),
  #         c(preds2[1,],rev(preds2[3,])),
  #         col='gray90')
  #lines(precip_pred.std,preds2[2,],lwd=2,col='orange')
  
  # polygon(c(precip_pred.std,rev(precip_pred.std)),
  #         c(preds3[1,],rev(preds3[3,])),
  #         col='gray90')
   lines(precip_pred.std,preds3[2,],lwd=2,col='red')
  
  # tmp = climate[climate$site==siteNames[i],]
  # points(tmp$p.std,tmp$seedlingJan/tmp$totalJan,pch=1,col='black')
  # text(2,.9,siteNames[i])
}

df=climate %>%
  dplyr::group_by(site,year,t,p) %>%
  dplyr::summarise(y=sum(seedlingJan),n=sum(totalJan)) %>%
  dplyr::mutate(g1=y/n)
siteNames = unique(df$site)

df.sum = df %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(p.mu = mean(t),g1.mu = mean(g1))

par(mfrow=c(1,1))
m1=lm(g1.mu~p.mu,data=df.sum)
plot(df.sum$p.mu,df.sum$g1.mu,xlim=c(0,12),ylim=c(0,.4))
abline(a=coef(m1)[1],b=coef(m1)[2],lwd=2)
for(i in 1:20){
  tmp=df[df$site==siteNames[i],]
  m.tmp=lm(g1~t,data=tmp)
  y.0=coef(m.tmp)[1]+min(tmp$t)*coef(m.tmp)[2]
  y.1=coef(m.tmp)[1]+max(tmp$t)*coef(m.tmp)[2]
  segments(x0=min(tmp$t),x1=max(tmp$t),y0=y.0,y1=y.1,lwd=1.5)
  #segments(x0=mean(tmp$p),y0=-1,y1=coef(m.tmp)[1]+mean(tmp$p)*coef(m.tmp)[2],lty='dotted')
  #segments(x0=mean(tmp$p),x1=-1,y0=coef(m.tmp)[1]+mean(tmp$p)*coef(m.tmp)[2],lty='dotted')
}

# -------------------------------------------------------------------
# Within subjects and between subjects
# -------------------------------------------------------------------

climate$t.group=climate$t-climate$mu.t
rrm1 = glmer(cbind(seedlingJan,totalJan-seedlingJan) ~ 1 + t.group + mu.t + (1 | site), data=climate, family = "binomial")
summary(rrm1)

summary(rrm1)
coef(rrm1)
fixef(rrm1)
ranef(rrm1)

t_pred=seq(min(climate$t),max(climate$t),.1)
out=boot::inv.logit(fixef(rrm1)["(Intercept)"]+fixef(rrm1)["mu.t"]*t_pred)
plot(t_pred,out,type='l',ylim=c(0,1),xlab= "Temperature",ylab= "P(germination)")
points(climate$t,climate$seedlingJan/climate$totalJan,col='gray80',pch=16,cex=.75)
lines(t_pred,out,lwd=2,col='black')

siteNames = unique(climate$site)


extract_pred=function(model,pred_seq){
  mu_alpha = fixef(model)["(Intercept)"]
  beta_w = fixef(model)["t.group"]
  beta_b = fixef(model)["mu.t"]
  alpha_j = ranef(model)$site$`(Intercept)`
  
  # beta = fixef(model)[2]
  # nu = ranef(model)$site$t
  
  alpha.site=mu_alpha+alpha_j[i]
 # beta.site=beta+nu[i]
  
  pred = alpha.site + beta_w*pred_seq + (beta_b-beta_w)*site_mean
  return(pred)
}

length.preds = 10
runList <- list()
pred.mat=matrix(NA,ncol=2,nrow=length.preds)
for(i in 1:20){
  site_temp=unique(climate$t[climate$site==siteNames[i]])
  t_pred = seq(min(site_temp),max(site_temp),length.out=length.preds)
  site_mean = unique(climate$mu.t[climate$site==siteNames[i]])
  pred.mat[,1] = t_pred
  pred.mat[,2] = extract_pred(rrm1,t_pred)
  runList[[i]] = pred.mat
}


t_pred=seq(min(climate$t),max(climate$t),.1)
out=boot::inv.logit(fixef(rrm1)["(Intercept)"]+fixef(rrm1)["mu.t"]*t_pred)
plot(t_pred,out,type='l',ylim=c(0,.25),xlab= "Temperature",ylab= "P(germination)")
points(climate$t,climate$seedlingJan/climate$totalJan,col='gray80',pch=16,cex=.75)
lines(t_pred,out,lwd=2,col='black')
for(i in 1:20){
  tmp = runList[[i]]
  lines(tmp[,1],boot::inv.logit(tmp[,2]),lwd=1)
}


# -------------------------------------------------------------------
# Within subjects and between subjects: precip
# -------------------------------------------------------------------

climate$p.group=climate$p-climate$mu.p
rrm2 = glmer(cbind(seedlingJan,totalJan-seedlingJan) ~ 1 + p.group + mu.p + (1 | site), data=climate, family = "binomial")
summary(rrm2)

summary(rrm2)
coef(rrm2)
fixef(rrm2)
ranef(rrm2)

p_pred=seq(min(climate$p),max(climate$p),5)
out=boot::inv.logit(fixef(rrm2)["(Intercept)"]+fixef(rrm2)["mu.p"]*p_pred)
plot(p_pred,out,type='l',ylim=c(0,1),xlab= "Precipitation",ylab= "P(germination)")
points(climate$p,climate$seedlingJan/climate$totalJan,col='gray80',pch=16,cex=.75)
lines(p_pred,out,lwd=2,col='black')

siteNames = unique(climate$site)


extract_pred=function(model,pred_seq){
  mu_alpha = fixef(model)["(Intercept)"]
  beta_w = fixef(model)["p.group"]
  beta_b = fixef(model)["mu.p"]
  alpha_j = ranef(model)$site$`(Intercept)`
  
  # beta = fixef(model)[2]
  # nu = ranef(model)$site$t
  
  alpha.site=mu_alpha+alpha_j[i]
  # beta.site=beta+nu[i]
  
  pred = alpha.site + beta_w*pred_seq + (beta_b-beta_w)*site_mean
  return(pred)
}

length.preds = 10
runList <- list()
pred.mat=matrix(NA,ncol=2,nrow=length.preds)
for(i in 1:20){
  site_precip=unique(climate$p[climate$site==siteNames[i]])
  p_pred = seq(min(site_precip),max(site_precip),length.out=length.preds)
  site_mean = unique(climate$mu.p[climate$site==siteNames[i]])
  pred.mat[,1] = p_pred
  pred.mat[,2] = extract_pred(rrm2,p_pred)
  runList[[i]] = pred.mat
}


p_pred=seq(min(climate$p),max(climate$p),.1)
out=boot::inv.logit(fixef(rrm2)["(Intercept)"]+fixef(rrm2)["mu.p"]*p_pred)
plot(p_pred,out,type='l',ylim=c(0,.25),xlab= "Precip",ylab= "P(germination)")
points(climate$p,climate$seedlingJan/climate$totalJan,col='gray80',pch=16,cex=.75)
lines(p_pred,out,lwd=2,col='black')
for(i in 1:20){
  tmp = runList[[i]]
  lines(tmp[,1],boot::inv.logit(tmp[,2]),lwd=1)
}


# -------------------------------------------------------------------
# Within subjects and between subjects: precip
# -------------------------------------------------------------------

climate$p.group=climate$p-climate$mu.p
rrm3 = glmer(cbind(seedlingJan,totalJan-seedlingJan) ~ 1 + p.group + mu.p + (1 + p.group | site), data=climate, family = "binomial")
summary(rrm3)

summary(rrm3)
coef(rrm3)
fixef(rrm3)
ranef(rrm3)

p_pred=seq(min(climate$p),max(climate$p),5)
out=boot::inv.logit(fixef(rrm3)["(Intercept)"]+fixef(rrm3)["mu.p"]*p_pred)
plot(p_pred,out,type='l',ylim=c(0,1),xlab= "Precipitation",ylab= "P(germination)")
points(climate$p,climate$seedlingJan/climate$totalJan,col='gray80',pch=16,cex=.75)
lines(p_pred,out,lwd=2,col='black')

siteNames = unique(climate$site)


extract_pred=function(model,pred_seq){
  mu_alpha = fixef(model)["(Intercept)"]
  beta_w = fixef(model)["p.group"]
  beta_b = fixef(model)["mu.p"]
  
  alpha_j = ranef(model)$site$`(Intercept)`
  beta_wj = ranef(model)$site$`p.group`
  
  alpha.site = mu_alpha+alpha_j[i]
  beta.site = beta_w + beta_wj[i]
  
  pred = alpha.site + beta.site*pred_seq + (beta_b-beta_w-beta_wj[i])*site_mean
  return(pred)
}

length.preds = 10
runList <- list()
pred.mat=matrix(NA,ncol=2,nrow=length.preds)
for(i in 1:20){
  site_precip=unique(climate$p[climate$site==siteNames[i]])
  p_pred = seq(min(site_precip),max(site_precip),length.out=length.preds)
  site_mean = unique(climate$mu.p[climate$site==siteNames[i]])
  pred.mat[,1] = p_pred
  pred.mat[,2] = extract_pred(rrm3,p_pred)
  runList[[i]] = pred.mat
}


p_pred=seq(min(climate$p),max(climate$p),.1)
out=boot::inv.logit(fixef(rrm3)["(Intercept)"]+fixef(rrm3)["mu.p"]*p_pred)
plot(p_pred,out,type='l',ylim=c(0,.25),xlab= "Precip",ylab= "P(germination)")
points(climate$p,climate$seedlingJan/climate$totalJan,col='gray80',pch=16,cex=.75)
lines(p_pred,out,lwd=2,col='black')
for(i in 1:20){
  tmp = runList[[i]]
  lines(tmp[,1],boot::inv.logit(tmp[,2]),lwd=1)
}

par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  plot(p_pred,out,type='l',ylim=c(0,.25),xlab= "Precip",ylab= "P(germination)")
  climate.tmp = climate[climate$site==siteNames[i],]
  points(climate.tmp$p,climate.tmp$seedlingJan/climate.tmp$totalJan,col='gray80',pch=16,cex=.75)
  lines(p_pred,out,lwd=2,col='black')
    tmp = runList[[i]]
    lines(tmp[,1],boot::inv.logit(tmp[,2]),lwd=1,col='red')
    text(x=160,y=0.02,siteNames[i])
}

plot(position$easting,ranef(rrm3)$site[,1])
plot(position$easting,ranef(rrm3)$site[,2])

fixef(rrm3)[2]-fixef(rrm3)[3]

# -------------------------------------------------------------------
# Within subjects and between subjects: temp
# -------------------------------------------------------------------

climate$t.group=climate$t-climate$mu.t
rrm4 = glmer(cbind(seedlingJan,totalJan-seedlingJan) ~ 1 + t.group + mu.t + (1 + t.group | site), data=climate, family = "binomial")
summary(rrm4)

summary(rrm4)
coef(rrm4)
fixef(rrm4)
ranef(rrm4)

t_pred=seq(min(climate$t),max(climate$t),.1)
out=boot::inv.logit(fixef(rrm4)["(Intercept)"]+fixef(rrm4)["mu.t"]*t_pred)
plot(t_pred,out,type='l',ylim=c(0,1),xlab= "Precipitation",ylab= "P(germination)")
points(climate$t,climate$seedlingJan/climate$totalJan,col='gray80',pch=16,cex=.75)
lines(t_pred,out,lwd=2,col='black')

siteNames = unique(climate$site)


extract_pred=function(model,pred_seq){
  mu_alpha = fixef(model)["(Intercept)"]
  beta_w = fixef(model)["t.group"]
  beta_b = fixef(model)["mu.t"]
  
  alpha_j = ranef(model)$site$`(Intercept)`
  beta_wj = ranef(model)$site$`t.group`
  
  alpha.site = mu_alpha+alpha_j[i]
  beta.site = beta_w + beta_wj[i]
  
  pred = alpha.site + beta.site*pred_seq + (beta_b-beta_w-beta_wj[i])*site_mean
  return(pred)
}

length.preds = 10
runList <- list()
pred.mat=matrix(NA,ncol=2,nrow=length.preds)
for(i in 1:20){
  site_temp=unique(climate$t[climate$site==siteNames[i]])
  t_pred = seq(min(site_temp),max(site_temp),length.out=length.preds)
  site_mean = unique(climate$mu.t[climate$site==siteNames[i]])
  pred.mat[,1] = t_pred
  pred.mat[,2] = extract_pred(rrm4,t_pred)
  runList[[i]] = pred.mat
}


t_pred=seq(min(climate$t),max(climate$t),.1)
out=boot::inv.logit(fixef(rrm4)["(Intercept)"]+fixef(rrm4)["mu.t"]*t_pred)
plot(t_pred,out,type='l',ylim=c(0,.25),xlab= "Precip",ylab= "P(germination)")
points(climate$t,climate$seedlingJan/climate$totalJan,col='gray80',pch=16,cex=.75)
lines(t_pred,out,lwd=2,col='black')
for(i in 1:20){
  tmp = runList[[i]]
  lines(tmp[,1],boot::inv.logit(tmp[,2]),lwd=1)
}

par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:20){
  plot(t_pred,out,type='l',ylim=c(0,.25),xlab= "Temp",ylab= "P(germination)")
  climate.tmp = climate[climate$site==siteNames[i],]
  points(climate.tmp$t,climate.tmp$seedlingJan/climate.tmp$totalJan,col='gray80',pch=16,cex=.75)
  lines(t_pred,out,lwd=2,col='black')
  tmp = runList[[i]]
  lines(tmp[,1],boot::inv.logit(tmp[,2]),lwd=1,col='red')
  text(x=160,y=0.02,siteNames[i])
}

plot(position$easting,ranef(rrm4)$site[,1])
plot(position$easting,ranef(rrm4)$site[,2])

fixef(rrm4)[2]-fixef(rrm4)[3]
