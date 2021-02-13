# to fit distributions with shape parameter add offset
# https://sourceforge.net/p/mcmc-jags/discussion/610037/thread/c21ef62a/

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
library(ggplot2)
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
  dplyr::filter(site=="EC") %>%
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
                       g_3 = c(0,0,0,0,0,1))

survivalData=data %>%
  dplyr::select(-seedlingJan) %>%
  tidyr::pivot_longer(cols=c(totalJan,intactOct)) %>%
  dplyr::left_join(refIndex,by=c("name","ageBags")) %>%
  dplyr::rename(response = name, y = value, siteSurvival = siteBags) %>%
  dplyr::mutate(gIndexSurvival = ageBags,
                yearSurvival = yearBags)

germinationData=data %>%
  dplyr::select(-c(seedStart,intactOct)) %>%
  dplyr::mutate(gIndex=ageBags) %>%
  dplyr::rename(siteGermination = siteBags,
                yearGermination = yearBags)

data <- tidybayes::compose_data(survivalData,germinationData)

data$n1 = dim(survivalData)[1]
data$n2 = dim(germinationData)[1]

refDf<-data.frame(year=data$yearGermination,gIndex=data$gIndex) %>%
  unique %>%
  dplyr::mutate(germinationReference=1:6)

data$yearRefGerm = refDf$year
data$indexRefGerm = refDf$gIndex
data$refGerm = refDf$germinationReference

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

dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/seedModelComparison/fullDataset/")

# set inits for JAGS
inits <- list()

# # -------------------------------------------------------------------
# # -------------------------------------------------------------------
# # Weibull residence time, constant germination
# # -------------------------------------------------------------------
# # -------------------------------------------------------------------
# # # tuning (n.adapt)
# jmWrCg = jags.model(paste0(dir,"jagsWrCg.R"), data = data, n.adapt = n.adapt,
#                     n.chains=2)
# #  inits = inits,
# #  n.chains = length(inits), n.adapt = n.adapt)
# 
# # burn-in (n.update)
# update(jmWrCg, n.iter = n.update)
# 
# samples.rjagsWrCg = coda.samples(jmWrCg,
#                                  variable.names = c("beta","a","b","g","logLik_y","logLik_g"),
#                                  n.iter = n.iterations, thin = n.thin)
# 
# MCMCsummary(samples.rjagsWrCg)
# 
# 
# LLmat <- MCMCchains(samples.rjagsWrCg,params="logLik_y")
# rel_n_eff <- relative_eff(exp(LLmat), chain_id = rep(1, each = 2*n.iterations))
# looWrCg_y <- loo(LLmat, r_eff = rel_n_eff, cores = 2)
# 
# print(looWrCg_y)
# plot(looWrCg_y)
# 
# LLmat <- MCMCchains(samples.rjagsWrCg,params="logLik_g")
# rel_n_eff <- relative_eff(exp(LLmat), chain_id = rep(1, each = 2*n.iterations))
# looWrCg_g <- loo(LLmat, r_eff = rel_n_eff, cores = 2)
# 
# print(looWrCg_g)
# plot(looWrCg_g)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Weibull residence time, age-dependent germination
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# # tuning (n.adapt)
rjags::load.module("glm", quiet = TRUE)

jmWrAg = jags.model(paste0(dir,"jagsWrAgRounds.R"), data = data, n.adapt = n.adapt,
                    n.chains=2)
#  inits = inits,
#  n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jmWrAg, n.iter = n.update)


# samples.rjagsWrAg = coda.samples(jmWrAg,
#                                  variable.names = c("beta","a","b","g","logLik_y","logLik_g"),
#                                  n.iter = n.iterations, thin = n.thin)

samples.rjagsWrAg = coda.samples(jmWrAg,
                                 variable.names = 
                                   c("a", "b", "mu0_s", "tau0_s", "mu_s", "tau_s", "beta", 
                                     "eta_surv","theta_g1","theta_s",
                                     "mu0_g","tau0_g","mu_g","tau_g",
                                     "logLik_y","logLik_g"),
                                 n.iter = n.iterations, thin = n.thin)
# 
 MCMCsummary(samples.rjagsWrAg,c("eta_surv"))
 MCMCchains(samples.rjagsWrAg,params="eta_surv")[,1] %>% hist
 MCMCchains(samples.rjagsWrAg,params="a")[,1] %>% hist
# df<-cbind(data$yearSurvival,apply(MCMCchains(samples.rjagsWrAg,"eta_surv"),2,mean))
# MCMCtrace(samples.rjagsWrAg,"a")
# MCMCchains(samples.rjagsWrAg,"a")
# plot(df[,1],-df[,2])
# MCMCsummary(samples.rjagsWrAg,params=c("mu0_g","tau0_g"))
# #MCMCtrace(samples.rjagsWrAg,params=c("tau0_g"))
# MCMCsummary(samples.rjagsWrAg,params=c("mu_g"))
# 
# MCMCsummary(samples.rjagsWrAg,params=c("mu0_s","tau0_s"))
# #MCMCtrace(samples.rjagsWrAg,params=c("tau0_g"))
# MCMCsummary(samples.rjagsWrAg,params=c("mu_s","tau_s"))
# 
# exp(-MCMCsummary(samples.rjagsWrAg,params="mu_s")$mean/MCMCsummary(samples.rjagsWrAg,params="a")$mean)
# 
# # plot 2
# library(RColorBrewer)
# colors <- colorRampPalette(brewer.pal(3, "Set1"))(3)
# 
# plot(rep(1,length=data$n2),data$seedlingJan/data$totalJan,
#      type='n',
#      xlim=c(0,3),ylim=c(0,1),
#      xlab="Years",ylab="P(germination)",
#      main="Germination model")
# 
# points(c(1,2,3),exp(MCMCsummary(samples.rjagsWrAg,params="mu0_g")$mean),pch=16)
# 
# points(c(1,1,1,2,2,3),exp(MCMCsummary(samples.rjagsWrAg,params="mu_g")$mean),
#        col=colors[c(1,2,3,1,2,1)],pch=1,cex=1)
# 
# for(i in 1:3){
# 
# points(data$gIndex[data$yearGermination==i],(data$seedlingJan/data$totalJan)[data$yearGermination==i],
#        xlim=c(0,3),ylim=c(0,1),pch=16,col=colors[i],cex=.5,
#        xlab="Years",ylab="P(germination)",
#        main="Germination model")
# }
# 
# exp(MCMCsummary(samples.rjagsWrAg,params="mu_g")$mean)
# germinationData %>%
#   dplyr::group_by(yearGermination,gIndex) %>%
#   dplyr::summarise(y = sum(seedlingJan),n=sum(totalJan)) %>%
#   dplyr::mutate(mu.p=y/n)


LLmat <- MCMCchains(samples.rjagsWrAg,params="logLik_y")
rel_n_eff <- relative_eff(exp(LLmat), chain_id = rep(1, each = 2*n.iterations))
looWrAg_y <- loo(LLmat, r_eff = rel_n_eff, cores = 2)

print(looWrAg_y)
plot(looWrAg_y)

LLmat <- MCMCchains(samples.rjagsWrAg,params="logLik_g")
rel_n_eff <- relative_eff(exp(LLmat), chain_id = rep(1, each = 2*n.iterations))
looWrAg_g <- loo(LLmat, r_eff = rel_n_eff, cores = 2)

print(looWrAg_g)
plot(looWrAg_g)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Model comparison
# -------------------------------------------------------------------
# -------------------------------------------------------------------

loo_compare(looWrCg_y, looWrAg_y)
loo_compare(looWrCg_g, looWrAg_g)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Extract parameters
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# 
# parmsWrCg = c()
# parmsWrCg.a = apply(MCMCchains(samples.rjagsWrCg,"a"),2,mean)
# parmsWrCg.b = apply(MCMCchains(samples.rjagsWrCg,"b"),2,mean)
# 
# funWrCg = function(t,par.a=parmsWrCg.a[1],par.b=parmsWrCg.b[1]){
#   a = par.a;
#   b = par.b;
#   theta = exp(-(t/b)^a);
#   return(theta)
# }

parmsWrAg = c()
parmsWrAg.a = apply(MCMCchains(samples.rjagsWrAg,"a"),2,mean)
parmsWrAg.b0 = exp(-boot::inv.logit(apply(MCMCchains(samples.rjagsWrAg,"mu0_s"),2,mean))/parmsWrAg.a)
parmsWrAg.b = exp(-boot::inv.logit(apply(MCMCchains(samples.rjagsWrAg,"mu_s"),2,mean))/parmsWrAg.a)

par(mfrow=c(1,2))
MCMCsummary(samples.rjagsWrAg,"a")
MCMCplot(samples.rjagsWrAg, 
         params = 'a',labels=unique(survivalData$siteSurvival))
MCMCplot(samples.rjagsWrAg, 
         params = 'b',labels=unique(survivalData$siteSurvival))

longevity = function(par.a=parmsWrAg.a[1],par.b=parmsWrAg.b[1]){
  a = par.a;
  b = par.b;
  longev = b*gamma(1+(1/a))
  return(longev)
}
tmp<-c()
for(i in 1:20){
  tmp[i]=longevity(par.a=parmsWrAg.a[i],par.b=parmsWrAg.b[i])
}
tmp=ifelse(tmp>1000,NA,tmp)
plot(site,tmp)
data.frame(site=unique(survivalData$siteSurvival),mean.longev=tmp/12)

funWrAg = function(t,par.a=parmsWrAg.a[1],par.b=parmsWrAg.b[1]){
  a = par.a;
  b = par.b;
  theta = exp(-(t/b)^a);
  return(theta)
}

gWrCg <- apply(MCMCchains(samples.rjagsWrCg,"g"),2,mean)
gWrAg<-matrix(NA,nrow=3,ncol=3)
for(i in 1:3){
  #index=grep(paste0("\\[",i,","),colnames(MCMCchains(samples.rjagsWrAg,"g")))
  index=grep(paste0("\\[1,",i),colnames(MCMCchains(samples.rjagsWrAg,"mu0_g")))
  
   # gWrAg[i,] <- apply(MCMCchains(samples.rjagsWrAg,"g")[,index],2,mean)
  gWrAg[i,] <- apply(MCMCchains(samples.rjagsWrAg,"mu0_g")[,index],2,mean)
  
}

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Plotting
# -------------------------------------------------------------------
# -------------------------------------------------------------------
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(4, "PuOr"))(20)

# -------------------------------------------------------------------
# Models
# -------------------------------------------------------------------
par(mfrow=c(1,2))
# plot 1
plot(rep(3,length=data$n1),data$y/data$seedStart,
     type='n',
     xlim=c(0,40),ylim=c(0,1),
     xlab="Months",ylab="P(survival)",
     main="Survival model")

for(i in 1:20){
  #lines(seq(0,36,by=0.1),funWrCg(seq(0,36,by=0.1),par.a=parmsWrCg.a[i],par.b=parmsWrCg.b[i]),lty='solid',col=colors[i])
  text(x=39,y=funWrAg(39,par.a=parmsWrAg.a[i],par.b=parmsWrAg.b[i]),
       unique(survivalData$siteSurvival)[i],cex=.5)
  lines(seq(0,36,by=0.1),funWrAg(seq(0,36,by=0.1),par.a=parmsWrAg.a[i],par.b=parmsWrAg.b[i]),lty='solid',col=colors[i])
}

# plot 2
plot(rep(1,length=data$n2),data$seedlingJan/data$totalJan,
     type='n',
     xlim=c(0,3),ylim=c(0,1),
     xlab="Years",ylab="P(germination)",
     main="Germination model")

points(data$yearGermination,data$seedlingJan/data$totalJan,
     xlim=c(0,3),ylim=c(0,1),pch=16,col='gray',cex=.5,
     xlab="Years",ylab="P(germination)",
     main="Germination model")

for(i in 1:20){
  #abline(h=gWrCg[i],lty=1,col=colors[i])
  # text(x=39,y=funWrCg(39,par.a=parmsWrCg.a[i],par.b=parmsWrCg.b[i]),
  #      unique(survivalData$siteSurvival)[i],cex=.5)
  #lines(seq(0,36,by=0.1),funWrAg(seq(0,36,by=0.1)),lty='dotted',col='orange')
  points(x=c(1,2,3),y=gWrAg[i,],col=colors[i],pch=16,cex=.5)
}


# -------------------------------------------------------------------
# Model vs. data
# -------------------------------------------------------------------
# mar default c(5, 4, 4, 2) + 0.1
par(mfrow=c(4,5),mar = c(2,2,1,1))

# for(i in 1:20){
# # plot 5
# plot(rep(3,length=data$n1),data$y/data$seedStart,
#      type='n',
#      xlim=c(0,36),ylim=c(0,1),
#      xlab="",ylab="",
#      main="",cex.main=.75)
# 
#   lines(seq(0,36,by=0.1),funWrCg(seq(0,36,by=0.1),par.a=parmsWrCg.a[i],par.b=parmsWrCg.b[i]),lty='solid',col='black',lwd=1)
#   
#   
# points(survivalData[data$siteSurvival==i,]$months,
#        survivalData[data$siteSurvival==i,]$y/survivalData[data$siteSurvival==i,]$seedStart,
#        pch=16,cex=.5,col='gray')
# 
# # df.summary = data.frame(cbind(months=data$months,
# #                               p=data$y/data$seedStart)) %>%
# #   dplyr::group_by(months) %>%
# #   dplyr::summarise(mu=mean(p))
# # points(df.summary$months,df.summary$mu,col=colors,pch=16)
# 
# 
# text(x=34.5,y=.98,
#      unique(survivalData$siteSurvival)[i])
# 
# # points(3,funWrCg(3),pch=16,col='black')
# # points(12,(1-gWrCg)*funWrCg(12),pch=16,col='black')
# # points(15,(1-gWrCg)*funWrCg(15),pch=16,col='black')
# # points(24,(1-gWrCg)*(1-gWrCg)*funWrCg(24),pch=16,col='black')
# # points(27,(1-gWrCg)*(1-gWrCg)*funWrCg(27),pch=16,col='black')
# # points(36,(1-gWrCg)*(1-gWrCg)*(1-gWrCg)*funWrCg(36),pch=16,col='black')
# }

# mar default c(5, 4, 4, 2) + 0.1
par(mfrow=c(4,5),mar = c(2,2,1,1))

for(i in 1:3){
# Weibull-residence time\n Age-dependent germination
plot(rep(3,length=data$n1),data$y/data$seedStart,
     type='n',
     xlim=c(0,36),ylim=c(0,1),
     xlab="",ylab="",
     main="",cex.main=.75)

  lines(seq(0,36,by=0.1),funWrAg(seq(0,36,by=0.1),par.a=parmsWrAg.a[1],par.b=parmsWrAg.b[i]),lty='solid',col='black',lwd=1)
  
  points(survivalData$months,
         survivalData$y/survivalData$seedStart,
         pch=16,cex=.5,col='gray')
  
  text(x=33,y=.98,
       unique(survivalData$siteSurvival)[i])
  
#   
df.summary = data.frame(cbind(months=survivalData[data$siteSurvival==i,]$months,
                               p=survivalData[data$siteSurvival==i,]$y/survivalData[data$siteSurvival==i,]$seedStart)) %>%
   dplyr::group_by(months) %>%
   dplyr::summarise(mu=mean(p))
 points(df.summary$months,df.summary$mu,col=colors,pch=16)

 surv = funWrAg(c(3,12,15,24,27,36),par.a=parmsWrAg.a[i],par.b=parmsWrAg.b[i])

points(3,surv[1],pch=16,col='black')
points(12,(1-gWrAg[i,1])*surv[2],pch=16,col='black')
points(15,(1-gWrAg[i,1])*surv[3],pch=16,col='black')
points(24,(1-gWrAg[i,1])*(1-gWrAg[i,2])*surv[4],pch=16,col='black')
points(27,(1-gWrAg[i,1])*(1-gWrAg[i,2])*surv[5],pch=16,col='black')
points(36,(1-gWrAg[i,1])*(1-gWrAg[i,2])*(1-gWrAg[i,3])*surv[6],pch=16,col='black')

}


# -------------------------------------------------------------------
# Predicted vs. observed
# -------------------------------------------------------------------
# mar default c(5, 4, 4, 2) + 0.1
par(mfrow=c(4,5),mar = c(2,2,1,1))
t=c(3,12,15,24,27,36)

mu.obs= survivalData %>%
  dplyr::group_by(siteSurvival,months) %>%
  dplyr::summarise(y.mu = mean(y)/100) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(siteIndex = group_indices_(., .dots = "siteSurvival"))


# plot 5
for(i in 1:20){
mu.pred.WrCg = c(funWrCg(3,par.a=parmsWrCg.a[i],par.b=parmsWrCg.b[i]),
                 (1-gWrCg[i])*funWrCg(12,par.a=parmsWrCg.a[i],par.b=parmsWrCg.b[i]),
                 (1-gWrCg[i])*funWrCg(15,par.a=parmsWrCg.a[i],par.b=parmsWrCg.b[i]),
                 (1-gWrCg[i])*(1-gWrCg[i])*funWrCg(24,par.a=parmsWrCg.a[i],par.b=parmsWrCg.b[i]),
                 (1-gWrCg[i])*(1-gWrCg[i])*funWrCg(27,par.a=parmsWrCg.a[i],par.b=parmsWrCg.b[i]),
                 (1-gWrCg[i])*(1-gWrCg[i])*(1-gWrCg[i])*funWrCg(36,par.a=parmsWrCg.a[i],par.b=parmsWrCg.b[i]))

plot(mu.pred.WrCg,mu.obs[mu.obs$siteIndex==i,]$y.mu,pch=16,
     xlim=c(0,1),ylim=c(0,1),
     xlab="",ylab="",
     main="",cex.main=.75)
abline(a=0,b=1)

text(x=.1,y=.98,
     unique(survivalData$siteSurvival)[i])
}
# plot 6
# Weibull residence time\n Age-dependent germination
par(mfrow=c(4,5),mar = c(2,2,1,1))

mu.obs= survivalData %>%
  dplyr::group_by(siteSurvival,months) %>%
  dplyr::summarise(y.mu = mean(y)/100) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(siteIndex = group_indices_(., .dots = "siteSurvival"))

for(i in 1:20){
  
surv = funWrAg(c(3,12,15,24,27,36),par.a=parmsWrAg.a[i],par.b=parmsWrAg.b[i])

mu.pred.WrAg = c(surv[1],
                 (1-gWrAg[i,1])*surv[2],
                 (1-gWrAg[i,1])*surv[3],
                 (1-gWrAg[i,1])*(1-gWrAg[i,2])*surv[4],
                 (1-gWrAg[i,1])*(1-gWrAg[i,2])*surv[5],
                 (1-gWrAg[i,1])*(1-gWrAg[i,2])*(1-gWrAg[i,3])*surv[6])

plot(mu.pred.WrAg,mu.obs[mu.obs$siteIndex==i,]$y.mu,pch=16,
     xlim=c(0,1),ylim=c(0,1),
     xlab="",ylab="",
     main="",cex.main=.75)
abline(a=0,b=1)

text(x=.1,y=.98,
     unique(survivalData$siteSurvival)[i])
}
# -------------------------------------------------------------------
# Predicted vs. observed
# -------------------------------------------------------------------
par(mfrow=c(1,2))

t=c(3,12,15,24,27,36)


# plot 5
mu.pred.WrCg = c(funWrCg(3),
                 (1-gWrCg)*funWrCg(12),
                 (1-gWrCg)*funWrCg(15),
                 (1-gWrCg)*(1-gWrCg)*funWrCg(24),
                 (1-gWrCg)*(1-gWrCg)*funWrCg(27),
                 (1-gWrCg)*(1-gWrCg)*(1-gWrCg)*funWrCg(36))

pred.obs = data.frame(y.pred=mu.pred.WrCg,months=t) %>%
  dplyr::left_join(data.frame(y=data$y/100,months=data$months) ,by="months")

plot(pred.obs$y.pred,pred.obs$y,pch=16,
     xlim=c(0,1),ylim=c(0,1),cex=.5,col='gray',
     xlab="Predicted intact",ylab="Observed intact",
     main="Weibull residence time\n Constant germination",cex.main=.75)
abline(a=0,b=1)

# plot 6

mu.pred.WrAg = c(funWrAg(3),
                 (1-gWrAg[1])*funWrAg(12),
                 (1-gWrAg[1])*funWrAg(15),
                 (1-gWrAg[1])*(1-gWrAg[2])*funWrAg(24),
                 (1-gWrAg[1])*(1-gWrAg[2])*funWrAg(27),
                 (1-gWrAg[1])*(1-gWrAg[2])*(1-gWrAg[3])*funWrAg(36))

pred.obs = data.frame(y.pred=mu.pred.WrAg,months=t) %>%
  dplyr::left_join(data.frame(y=data$y/100,months=data$months) ,by="months")

plot(pred.obs$y.pred,pred.obs$y,pch=16,
     xlim=c(0,1),ylim=c(0,1),cex=.5,col='gray',
     xlab="Predicted intact",ylab="Observed intact",
     main="Weibull residence time\n Age-dependent germination",cex.main=.75)
abline(a=0,b=1)

