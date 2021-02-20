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

set.seed(10)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Read data
# -------------------------------------------------------------------
# -------------------------------------------------------------------

data = readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/seedBagMaster.RDS")

data = data %>%
  dplyr::filter(site=="CF") %>%
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


data <- tidybayes::compose_data(survivalData,germinationData)

data$n1 = dim(survivalData)[1]
data$n2 = dim(germinationData)[1]

refDf<-data.frame(year=data$yearGermination,gIndex=data$gIndex) %>%
  unique %>%
  dplyr::mutate(germinationReference=1:6)

data$yearRefGerm = refDf$year
data$indexRefGerm = refDf$gIndex
data$refGerm = refDf$germinationReference

d=survivalData %>%
  dplyr::mutate(index=1:dim(survivalData)[1]) %>%
  dplyr::group_by(months,yearBags) %>%
  dplyr::top_n(n=1)  %>%
  dplyr::ungroup() %>%
  dplyr::rename(seedStart_pred=seedStart,
                siteSurvival_pred=siteSurvival,
                yearSurvival_pred=yearSurvival,
                months_pred=months,
                compIndex_pred=compIndex) %>%
  dplyr::select(seedStart_pred,siteSurvival_pred,yearSurvival_pred,months_pred,compIndex_pred)

d
d<-tidybayes::compose_data(d)
names(d)[8]="n_pred"
data=c(data,d)

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
rjags::load.module("glm", quiet = TRUE)

jmWrAg = jags.model(paste0(dir,"binomialLikelihood-logLink-wrDecay-2.R"), data = data, n.adapt = n.adapt,
                    n.chains=2)

# burn-in (n.update)
update(jmWrAg, n.iter = n.update)

samples.rjags = coda.samples(jmWrAg,
                                 variable.names = 
                                   # GERMINATION
                                   c("mu0_g","tau0_g","mu_g","tau_g",
                                     "sigma0_g","sigma_g","g",
                                     "mu0_pred","tau0_pred","mu_pred","tau_pred",
                                      "sigma0_pred","sigma_pred","g_pred",
                                     #PRIOR PREDICTIVE
                                     "y_pred","y_prior.pred",
                                     # SURVIVAL
                                     "mu0_s","mu_s","inv.b","theta_c",
                                     "sigma0_s","sigma_s",
                                     "beta","beta_mod","alpha_s","theta_s",
                                     "mu","a",
                                     # POSTERIOR PREDICTIVE
                                     "y_sim","y_surv","mu_survival",
                                     "logLik_y"),
                                 n.iter = n.iterations, thin = n.thin)


# exponential
jmExAg = jags.model(paste0(dir,"binomialLikelihood-logLink-expDecay-3.R"), data = data, n.adapt = n.adapt,
                    n.chains=2)

# burn-in (n.update)
update(jmExAg, n.iter = n.update)

samples.rjags.jmExAg = coda.samples(jmExAg,
                             variable.names = 
                               # GERMINATION
                               c("mu0_g","tau0_g","mu_g","tau_g",
                                 "sigma0_g","sigma_g","g",
                                 "mu0_pred","tau0_pred","mu_pred","tau_pred",
                                 "sigma0_pred","sigma_pred","g_pred",
                                 #PRIOR PREDICTIVE
                                 "y_pred","y_prior.pred",
                                 # SURVIVAL
                                 "mu0_s","mu_s","inv.b","theta_c",
                                 "sigma0_s","sigma_s",
                                 "beta","beta_mod","alpha_s","theta_s",
                                 "mu","a",
                                 # POSTERIOR PREDICTIVE
                                 "y_sim","y_surv","mu_survival",
                                 "logLik_y"),
                             n.iter = n.iterations, thin = n.thin)

 
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Model comparison
# -------------------------------------------------------------------
# -------------------------------------------------------------------
LLmat.wb <- MCMCchains(samples.rjags,params="logLik_y")
rel_n_eff.wb <- relative_eff(exp(LLmat.wb), chain_id = rep(1, each = 2*n.iterations))
looWb <- loo(LLmat.wb, r_eff = rel_n_eff.wb, cores = 2)

print(looWb)
plot(looWb)

LLmat.ne <- MCMCchains(samples.rjags.jmExAg,params="logLik_y")
rel_n_eff.ne <- relative_eff(exp(LLmat.ne), chain_id = rep(1, each = 2*n.iterations))
looNe <- loo(LLmat.ne, r_eff = rel_n_eff.ne, cores = 2)

print(looNe)
plot(looNe)

loo_compare(looNe,looWb)


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Evaluate priors
# -------------------------------------------------------------------
# -------------------------------------------------------------------

hist.chain=function(x){
  name<-colnames(x)
  tmp = hist(x,breaks=100,plot=FALSE)
  tmp$xname = name
  plot(tmp,freq=FALSE)
}

MCMCsummary(samples.rjags,params=c("mu0_s","sigma0_s"))
MCMCsummary(samples.rjags,params=c("mu_s","sigma_s"))

MCMCsummary(samples.rjags,params=c("mu0_g","sigma0_g"))
MCMCsummary(samples.rjags,params=c("mu_g","sigma_g"))

MCMCsummary(samples.rjags.jmExAg,params=c("mu0_s","sigma0_s"))
MCMCsummary(samples.rjags.jmExAg,params=c("mu_s","sigma_s"))

MCMCsummary(samples.rjags.jmExAg,params=c("mu0_g","sigma0_g"))
MCMCsummary(samples.rjags.jmExAg,params=c("mu_g","sigma_g"))


MCMCsummary(samples.rjags,params=c("a"))

## Priors on hyperparameters
par(mfrow=c(2,2))
hist.chain(MCMCchains(samples.rjags,params="mu0_s")[,1,drop=FALSE])
hist.chain(MCMCchains(samples.rjags,params="sigma0_s")[,1,drop=FALSE])

hist.chain(MCMCchains(samples.rjags,params="mu_s")[,1,drop=FALSE])
hist.chain(MCMCchains(samples.rjags,params="sigma_s")[,1,drop=FALSE])

par(mfrow=c(2,2))
hist.chain(MCMCchains(samples.rjags,params="mu0_g")[,1,drop=FALSE])
hist.chain(MCMCchains(samples.rjags,params="sigma0_g")[,1,drop=FALSE])

hist.chain(MCMCchains(samples.rjags,params="mu_g")[,1,drop=FALSE])
hist.chain(MCMCchains(samples.rjags,params="sigma_g")[,1,drop=FALSE])

## Priors on germination history
theta_c=MCMCchains(samples.rjags,params="theta_c")
par(mfrow=c(2,3))
for(i in 1:6){
  hist.chain(theta_c[,i])
}

## Priors on lambda
inv.b=MCMCchains(samples.rjags,params="inv.b")
par(mfrow=c(2,3))
for(i in 1:6){
  hist.chain(inv.b[,i])
}

a=MCMCchains(samples.rjags,params="a")
par(mfrow=c(1,1))
for(i in 1:1){
  hist.chain(a[,i])
}
abline(v=1,lty='dotted')
# eta_surv[i] ~ dnorm(mu_s[siteSurvival[i],yearSurvival[i]],tau_s[siteSurvival[i],yearSurvival[i]])
# inv.b[i] <- exp(-(eta_surv[i])/(a[siteSurvival[i]]))  
# 


## Priors on latent survival
mu_s=MCMCchains(samples.rjags,params="mu_survival")
par(mfrow=c(2,3))
for(i in 1:6){
  hist.chain(mu_s[,i])
}

d=survivalData %>%
  dplyr::mutate(index=1:dim(survivalData)[1]) %>%
  dplyr::group_by(months,yearBags) %>%
  dplyr::top_n(n=1) 
index=d$index


## Priors on latent survival*germination history
mu=MCMCchains(samples.rjags,params="mu")
par(mfrow=c(2,3))
for(i in 1:6){
  hist.chain(mu[,i])
}


## Priors on latent prob
theta_s=MCMCchains(samples.rjags,params="theta_s")
par(mfrow=c(2,3))
for(i in 1:6){
  hist.chain(theta_s[,i])
}

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Evaluate prior predictive
# -------------------------------------------------------------------
# -------------------------------------------------------------------
y_prior.pred=MCMCchains(samples.rjags,params="y_prior.pred")
par(mfrow=c(2,3))
for(i in 1:6){
  hist.chain(y_prior.pred[,i])
}

# data$months_pred
# 
# d=survivalData %>%
#   dplyr::mutate(index=1:dim(survivalData)[1]) %>%
#   dplyr::group_by(months,yearBags) %>%
#   dplyr::top_n(n=1) 
# index=d$index

#index=c(1,2,21,22,41,42,59,60,77,78,99,100)

par(mfrow=c(2,3))
samples=sample(1:dim(y_prior.pred)[1],1000)

df.plot<-survivalData%>%
  dplyr::group_by(months,yearBags) %>%
  dplyr::mutate(p.mu=mean(y))

#year1
plot(data$months_pred,y_prior.pred[1,],type='n',ylim=c(0,100))
for(i in 1:length(samples)){
  points(data$months_pred[1:6]+rnorm(1,0,sd=.01),
         y_prior.pred[samples[i],1:6],
         pch=16,cex=.5)
}
points(df.plot[df.plot$yearBags==2006,]$months,
       df.plot[df.plot$yearBags==2006,]$p.mu,ylim=c(0,100),
       pch=16,cex=1.5,
       col='gray')

#year2
plot(data$months_pred,y_prior.pred[1,],type='n',ylim=c(0,100))
for(i in 1:length(samples)){
  points(data$months_pred[7:10]+rnorm(1,0,sd=.01),
         y_prior.pred[samples[i],7:10],
         pch=16,cex=.5)
}
points(df.plot[df.plot$yearBags==2007,]$months,
       df.plot[df.plot$yearBags==2007,]$p.mu,ylim=c(0,100),
       pch=16,cex=1.5,
       col='gray')
#year3
plot(data$months_pred,y_prior.pred[1,],type='n',ylim=c(0,100))
for(i in 1:length(samples)){
  points(data$months_pred[11:12]+rnorm(1,0,sd=.01),
         y_prior.pred[samples[i],11:12],
         pch=16,cex=.5)
}
points(df.plot[df.plot$yearBags==2008,]$months,
       df.plot[df.plot$yearBags==2008,]$p.mu,ylim=c(0,100),
       pch=16,cex=1.5,
       col='gray')

## COMPARE TWO MODELS

y_prior.pred.jmExAg=MCMCchains(samples.rjags.jmExAg,params="y_prior.pred")

samples=sample(1:dim(y_prior.pred.jmExAg)[1],1000)

#year1
plot(data$months_pred,y_prior.pred[1,],type='n',ylim=c(0,100))
for(i in 1:length(samples)){
  points(data$months_pred[1:6]+rnorm(1,0,sd=.01),
         y_prior.pred.jmExAg[samples[i],1:6],
         pch=16,cex=.5)
}
points(df.plot[df.plot$yearBags==2006,]$months,
       df.plot[df.plot$yearBags==2006,]$p.mu,ylim=c(0,100),
       pch=16,cex=1.5,
       col='gray')

#year2
plot(data$months_pred,y_prior.pred[1,],type='n',ylim=c(0,100))
for(i in 1:length(samples)){
  points(data$months_pred[7:10]+rnorm(1,0,sd=.01),
         y_prior.pred.jmExAg[samples[i],7:10],
         pch=16,cex=.5)
}
points(df.plot[df.plot$yearBags==2007,]$months,
       df.plot[df.plot$yearBags==2007,]$p.mu,ylim=c(0,100),
       pch=16,cex=1.5,
       col='gray')
#year3
plot(data$months_pred,y_prior.pred[1,],type='n',ylim=c(0,100))
for(i in 1:length(samples)){
  points(data$months_pred[11:12]+rnorm(1,0,sd=.01),
         y_prior.pred.jmExAg[samples[i],11:12],
         pch=16,cex=.5)
}
points(df.plot[df.plot$yearBags==2008,]$months,
       df.plot[df.plot$yearBags==2008,]$p.mu,ylim=c(0,100),
       pch=16,cex=1.5,
       col='gray')

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Evaluate prior predictive of survival
# -------------------------------------------------------------------
# -------------------------------------------------------------------

mu0_s.exp=MCMCchains(samples.rjags.jmExAg,params="mu0_s")
mu_s.exp=MCMCchains(samples.rjags.jmExAg,params="mu_s")

b0=exp(mu0_s.exp)
b0.parm=apply(b0,2,quantile,probs=c(0.025,.5,.975))

b=exp(mu_s.exp)
b.parm=apply(b,2,quantile,probs=c(0.025,.5,.975))

# par(mfrow=(c(1,3)))
# mu_s=MCMCchains(samples.rjags.jmExAg,params="mu_s")
# plot(density(mu_s[,1]),xlim=c(-2,2),ylim=c(0,3));
# lines(density(mu_s[,2]))
# lines(density(mu_s[,3]))
# 
# par(mfrow=c(1,1))
# plot(density(b[,1]),xlim=c(0,6),ylim=c(0,3))
# lines(density(b[,2]),xlim=c(0,3))
# lines(density(b[,3]),xlim=c(0,3))

t=seq(0,max(data$months),.01)

# Plot survival curve
par(mfrow=c(1,1))
plot(t,exp(-b0.parm[2,1]*t),
     xlim=c(0,max(data$months)),ylim=c(0,1),type='l')
lines(t,exp(-b0.parm[1,1]*t),lty='dotted')
lines(t,exp(-b0.parm[3,1]*t),lty='dotted')

# colors=c('black',"purple","orange")
# for(i in 1:3){
#   lines(t,exp(-b0.parm[2,i]*t),
#         xlim=c(0,max(data$months)),ylim=c(0,1),type='l',
#         col=colors[i])
# }
#lines(t,exp(-median(b0)*t),lty='dotted')

par(mfrow=c(2,3))
## AGE 1
plot(t,exp(-b.parm[2,1]*t),
     xlim=c(0,max(data$months)),ylim=c(0,1),type='l')
for(i in 1:100){
  lines(t,exp(-sample(b[,1],1)*t),
        lwd=.5,col='lightgray')
}
lines(t,exp(-b.parm[2,1]*t),lwd=2)
lines(t,exp(-median(b0)*t),lty='dotted')

# df<-survivalData[survivalData$yearSurvival==2006,]
# points(df$months,df$y/df$seedStart,pch=16,cex=.75)

## AGE 2
plot(t,exp(-b.parm[2,2]*t),
     xlim=c(0,max(data$months)),ylim=c(0,1),type='l')
for(i in 1:100){
  lines(t,exp(-sample(b[,2],1)*t),
        lwd=.5,col='lightgray')
}
lines(t,exp(-b.parm[2,2]*t),lwd=2)
lines(t,exp(-median(b0)*t),lty='dotted')

# df<-survivalData[survivalData$yearSurvival==2007,]
# points(df$months,df$y/df$seedStart,pch=16,cex=.75)

## AGE 3
plot(t,exp(-b.parm[2,3]*t),
     xlim=c(0,max(data$months)),ylim=c(0,1),type='l')
for(i in 1:100){
  lines(t,exp(-sample(b[,3],1)*t),
        lwd=.5,col='lightgray')
}
lines(t,exp(-b.parm[2,3]*t),lwd=2)
lines(t,exp(-median(b0)*t),lty='dotted')

# df<-survivalData[survivalData$yearSurvival==2008,]
# points(df$months,df$y/df$seedStart,pch=16,cex=.75)

# WEIBULL

a.wb=MCMCchains(samples.rjags,params="a")
b.wb=MCMCchains(samples.rjags,params="mu_s")
b0.wb=MCMCchains(samples.rjags,params="mu0_s")

inv.b0=cbind(exp(-b0.wb/a.wb))
inv.b0.parm=apply(inv.b0,2,quantile,probs=c(0.025,.5,.975))

inv.b=cbind(exp(-b.wb[,1]/a.wb),exp(-b.wb[,2]/a.wb),exp(-b.wb[,3]/a.wb))
inv.b.parm=apply(inv.b,2,quantile,probs=c(0.025,.5,.975))

a.parm=apply(a.wb,2,quantile,probs=c(0.025,.5,.975))

t=seq(0,max(data$months),.01)

par(mfrow=c(1,1))
plot(t,exp(-(t/inv.b0.parm[2,1])^a.parm[2,1]),
     xlim=c(0,max(data$months)),ylim=c(0,1),type='l')
lines(t,exp(-(t/inv.b0.parm[1,1])^a.parm[2,1]),lty='dotted')
lines(t,exp(-(t/inv.b0.parm[3,1])^a.parm[2,1]),lty='dotted')


# sample from the same rows of the posterior to keep correlations between parameters
# par(mfrow=c(1,3))
# plot(a.wb,b.wb[,1])
# plot(a.wb,b.wb[,2])
# plot(a.wb,b.wb[,3])

par(mfrow=c(1,3))
## AGE 1
plot(t,exp(-(t/inv.b.parm[2,1])^a.parm[2,1]),
     xlim=c(0,max(data$months)),ylim=c(0,1),type='l')
for(i in 1:100){
  tmp=sample(1:dim(a.wb)[1],1)
  lines(t,exp(-(t/inv.b[tmp,1])^a.wb[tmp,1]),
        lwd=.5,col=ifelse(a.wb[tmp,1]>1,"purple","orange"))
}
lines(t,exp(-(t/inv.b.parm[2,1])^a.parm[2,1]),lwd=2)

# df<-survivalData[survivalData$yearSurvival==2006,]
# points(df$months,df$y/df$seedStart,pch=16,cex=.75)

## AGE 2
plot(t,exp(-(t/inv.b.parm[2,2])^a.parm[2,1]),
     xlim=c(0,max(data$months)),ylim=c(0,1),type='l')
for(i in 1:100){
  tmp=sample(1:dim(a.wb)[1],1)
  lines(t,exp(-(t/inv.b[tmp,2])^a.wb[tmp,1]),
        lwd=.5,col=ifelse(a.wb[tmp,1]>1,"purple","orange"))
}
lines(t,exp(-(t/inv.b.parm[2,2])^a.parm[2,1]),lwd=2)

# 
# df<-survivalData[survivalData$yearSurvival==2007,]
# points(df$months,df$y/df$seedStart,pch=16,cex=.75)


## AGE 3
plot(t,exp(-(t/inv.b.parm[2,3])^a.parm[2,1]),
     xlim=c(0,max(data$months)),ylim=c(0,1),type='l')
for(i in 1:100){
  tmp=sample(1:dim(a.wb)[1],1)
  lines(t,exp(-(t/inv.b[tmp,3])^a.wb[tmp,1]),
        lwd=.5,col=ifelse(a.wb[tmp,1]>1,"purple","orange"))
}
lines(t,exp(-(t/inv.b.parm[2,3])^a.parm[2,1]),lwd=2)
# 
# df<-survivalData[survivalData$yearSurvival==2008,]
# points(df$months,df$y/df$seedStart,pch=16,cex=.75)




# 
# samples=sample(1:dim(y_prior.pred.jmExAg)[1],100)
# 
# par(mfrow=c(1,2))
# plot(data$y,apply(y_prior.pred.jmExAg,2,mean),type='n');abline(a=0,b=1)
# for(i in 1:length(samples)){
#   points(data$y,y_prior.pred.jmExAg[samples[i],],pch=16,cex=.5)
# }
# 
# plot(data$y,apply(y_prior.pred,2,mean),type='n');abline(a=0,b=1)
# for(i in 1:length(samples)){
#   points(data$y,y_prior.pred[samples[i],],pch=16,cex=.5)
# }
# 
# plot(data$y,apply(y_prior.pred,2,mean),type='n');abline(a=0,b=1)
# for(i in 1:length(samples)){
#   points(y_prior.pred.jmExAg[samples[i],],y_prior.pred[samples[i],],pch=16,cex=.5)
# }
# 
# 
# plot(data$y,apply(y_prior.pred.jmExAg,2,mean));abline(a=0,b=1)
# plot(data$y,apply(y_prior.pred,2,mean));abline(a=0,b=1)
# 
# 
# plot(data$months+rnorm(length(data$months),0,0.01),
#      data$y,ylim=c(0,100),pch=16,cex=.5,
#      col=data$yearBags)
# 
# df.plot<-survivalData%>%
#   dplyr::group_by(months,yearBags) %>%
#   dplyr::mutate(p.mu=mean(y))
# points(df.plot$months,
#        df.plot$p.mu,ylim=c(0,100),
#      pch=16,cex=1.5,
#      col=df.plot$yearBags)
# 
# y_surv=MCMCchains(samples.rjags,params="y_surv")
# par(mfrow=c(2,3))
# for(i in 1:6){
#   hist.chain(y_surv[,i])
# }
# 
# par(mfrow=c(1,1))
# index=sample(1:dim(y_surv)[1],1000)
# plot(unique(data$months),y_surv[1,1:6],type='n',ylim=c(0,100))
# for(i in 1:length(index)){
#   points(unique(data$months)+rnorm(1,0,.01),
#          y_surv[index[i],1:6],
#          pch=16,cex=.5)
# }
# 
# 
# 
# # 
# y_surv = MCMCchains(samples.rjagsWrAg,params='y_surv')
# 
# 
# par(mfrow=c(1,1))
# plot(rep(data$y[1],100),sample(y_surv[,1],100),type='n',
#      xlim=c(0,100),
#      ylim=c(0,100))
# for(i in 1:dim(y_surv)[2]){
#   points(rep(data$y[i],100),sample(y_surv[,i],100),
#          pch=16,cex=.5,col='lightgray')
# }
# abline(a=0,b=1)
# 
# plot(apply(y_surv[,data$months==1],2,median),
#      data$y[data$months==1],
#      type='n',
#      xlim=c(0,100),
#      ylim=c(0,100))
# abline(a=0,b=1)
#   points(apply(y_surv[,data$months==1],2,median),
#          data$y[data$months==1],
#          pch=16,cex=1,col='lightgray')
# 
# 
# ## GERMINATION
#   
#   y_pred=MCMCchains(samples.rjagsWrAg,params="y_pred")
#   y_sim=MCMCchains(samples.rjagsWrAg,params="y_sim")
#   
#   
#   par(mfrow=c(1,2))
#   plot(data$seedlingJan,y_pred[sample(dim(y_pred)[1],1),],
#        type='n',xlim=c(0,100),ylim=c(0,100),
#        xlab="Observed data",ylab="Simulated data")
#   for(i in 1:100){
#     points(data$seedlingJan,y_pred[sample(dim(y_pred)[1],1),],pch=16,cex=.5,
#            col=data$yearGermination)
#   }
#   abline(a=0,b=1)
#   
#   plot(data$seedlingJan,y_sim[sample(dim(y_sim)[1],1),],
#        type='n',xlim=c(0,100),ylim=c(0,100),
#        xlab="Observed data",ylab="Simulated data")
#   for(i in 1:100){
#     points(data$seedlingJan,y_sim[sample(dim(y_sim)[1],1),],pch=16,cex=.5,
#            col=data$yearGermination)
#   }
#   abline(a=0,b=1)
#   
#   # Effect of survival data on germination estimate
#   
#   par(mfrow=c(3,1))
#   plot(density(MCMCchains(samples.rjagsWrAg,params="mu0_g")[,1]),xlim=c(-5,5),ylim=c(0,1.25))
#   lines(density(MCMCchains(samples.rjagsWrAg,params="mu0_pred")[,1]),col='red')
#   
#   lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,1]),col='black',lty='dotted')
#   lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,2]),col='black',lty='dotted')
#   lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,3]),col='black',lty='dotted')
#   
#   lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,1]),col='red',lty='dotted')
#   lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,2]),col='red',lty='dotted')
#   lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,3]),col='red',lty='dotted')
#   
#   plot(density(MCMCchains(samples.rjagsWrAg,params="mu0_g")[,2]),xlim=c(-5,5),ylim=c(0,1.25))
#   lines(density(MCMCchains(samples.rjagsWrAg,params="mu0_pred")[,2]),col='red')
#   
#   lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,4]),col='black',lty='dotted')
#   lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,5]),col='black',lty='dotted')
#   
#   lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,4]),col='red',lty='dotted')
#   lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,5]),col='red',lty='dotted')
#   
#   plot(density(MCMCchains(samples.rjagsWrAg,params="mu0_g")[,3]),xlim=c(-5,5),ylim=c(0,1.25))
#   lines(density(MCMCchains(samples.rjagsWrAg,params="mu0_pred")[,3]),col='red')
#   
#   lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,6]),col='black',lty='dotted')
#   
#   lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,6]),col='red',lty='dotted')
#   
#   
#   
#   
#   
#   
# y_pred <- MCMCchains(samples.rjagsWrAg,params="y_pred")
# 
# plot(data$seedlingJan,y_pred[sample(dim(y_pred)[1],1),],
#      type='n',xlim=c(0,100),ylim=c(0,100),
#      xlab="Observed data",ylab="Simulated data")
# for(i in 1:100){
#   points(data$seedlingJan,y_pred[sample(dim(y_pred)[1],1),],pch=16,cex=.5,
#          col=data$yearGermination)
# }
# abline(a=0,b=1)
# 
# #sigma_g=MCMCchains(samples.rjagsWrAg,params="sigma_g")[,1]
# #hist(sigma_g,breaks=100,freq=FALSE,ylim=c(0,2),xlim=c(0,5))
# 
# 
# 
# ## SIGMA
# 
# plot(density(MCMCchains(samples.rjagsWrAg,params="sigma0_g")[,1]),xlim=c(0,4),ylim=c(0,1.5))
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma0_pred")[,1]),col='red')
# 
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,1]),col='black',lty='dotted')
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,2]),col='black',lty='dotted')
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,3]),col='black',lty='dotted')
# 
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,1]),col='red',lty='dotted')
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,2]),col='red',lty='dotted')
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,3]),col='red',lty='dotted')
# 
# plot(density(MCMCchains(samples.rjagsWrAg,params="sigma0_g")[,2]),xlim=c(0,4),ylim=c(0,1.5))
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma0_pred")[,2]),col='red')
# 
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,4]),col='black',lty='dotted')
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,5]),col='black',lty='dotted')
# 
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,4]),col='red',lty='dotted')
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,5]),col='red',lty='dotted')
# 
# plot(density(MCMCchains(samples.rjagsWrAg,params="sigma0_g")[,3]),xlim=c(0,4),ylim=c(0,1.5))
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma0_pred")[,3]),col='red')
# 
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,6]),col='black',lty='dotted')
# 
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,6]),col='red',lty='dotted')
# 
# 
# ## SUMMARY
# MCMCsummary(samples.rjagsWrAg,params=c("tau0_g","tau_g"))
# MCMCsummary(samples.rjagsWrAg,params=c("sigma0_g","sigma_g"))
# MCMCsummary(samples.rjagsWrAg,params=c("mu0_g","mu_g"))
