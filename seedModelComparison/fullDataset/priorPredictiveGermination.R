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
  dplyr::filter(site=="KYE") %>%
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

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Weibull residence time, age-dependent germination
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# # tuning (n.adapt)
rjags::load.module("glm", quiet = TRUE)

jmWrAg = jags.model(paste0(dir,"jagsWrAgPriorPredictiveGermination.R"), data = data, n.adapt = n.adapt,
                    n.chains=2)

# burn-in (n.update)
update(jmWrAg, n.iter = n.update)

samples.rjagsWrAg = coda.samples(jmWrAg,
                                 variable.names = 
                                   c("mu0_g","tau0_g","mu_g","tau_g",
                                     "sigma0_g","sigma_g","g",
                                     "mu0_pred","tau0_pred","mu_pred","tau_pred",
                                      "sigma0_pred","sigma_pred","g_pred",
                                     "y_pred",'y_sim'),
                                 n.iter = n.iterations, thin = n.thin)

# MCMCsummary(samples.rjagsWrAg,params="y_g")
 MCMCsummary(samples.rjagsWrAg,params=c("tau0_g","tau_g"))
 MCMCsummary(samples.rjagsWrAg,params=c("mu0_g","mu_g"))
 
y_pred=MCMCchains(samples.rjagsWrAg,params="y_pred")
theta_g=MCMCchains(samples.rjagsWrAg,params="g")
theta_pred=MCMCchains(samples.rjagsWrAg,params="g_pred")

dim(y_pred)[1]

par(mfrow=c(3,3))
hist(y_pred[sample(dim(y_pred)[1],100),],breaks=100,freq=FALSE)

hist(theta_pred[,1:49],breaks=100,freq=FALSE)
#hist(theta_g[,1:49],breaks=100,freq=FALSE,add=TRUE)

#hist(data$seedlingJan,breaks=25,add=TRUE,col='red',freq=FALSE)

plot(data$seedlingJan,y_pred[sample(dim(y_pred)[1],1),],
     type='n',xlim=c(0,100),ylim=c(0,100),
     xlab="Observed data",ylab="Simulated data")
for(i in 1:100){
  points(data$seedlingJan,y_pred[sample(dim(y_pred)[1],1),],pch=16,cex=.5,
         col=data$yearGermination)
}
abline(a=0,b=1)

#sigma_g=MCMCchains(samples.rjagsWrAg,params="sigma_g")[,1]
#hist(sigma_g,breaks=100,freq=FALSE,ylim=c(0,2),xlim=c(0,5))



#par(mfrow=c(1,3))
plot(density(MCMCchains(samples.rjagsWrAg,params="mu0_g")[,1]),xlim=c(-4,4))
lines(density(MCMCchains(samples.rjagsWrAg,params="mu0_pred")[,1]),col='red')

lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,1]),col='black',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,2]),col='black',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,3]),col='black',lty='dotted')

lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,1]),col='red',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,2]),col='red',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,3]),col='red',lty='dotted')

plot(density(MCMCchains(samples.rjagsWrAg,params="mu0_g")[,2]),xlim=c(-4,4))
lines(density(MCMCchains(samples.rjagsWrAg,params="mu0_pred")[,2]),col='red')

lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,4]),col='black',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,5]),col='black',lty='dotted')

lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,4]),col='red',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,5]),col='red',lty='dotted')

plot(density(MCMCchains(samples.rjagsWrAg,params="mu0_g")[,3]),xlim=c(-4,4))
lines(density(MCMCchains(samples.rjagsWrAg,params="mu0_pred")[,3]),col='red')

lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,6]),col='black',lty='dotted')

lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,6]),col='red',lty='dotted')


## SIGMA

plot(density(MCMCchains(samples.rjagsWrAg,params="sigma0_g")[,1]))
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma0_pred")[,1]),col='red')

lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,1]),col='black',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,2]),col='black',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,3]),col='black',lty='dotted')

lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,1]),col='red',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,2]),col='red',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,3]),col='red',lty='dotted')

plot(density(MCMCchains(samples.rjagsWrAg,params="sigma0_g")[,2]))
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma0_pred")[,2]),col='red')

lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,4]),col='black',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,5]),col='black',lty='dotted')

lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,4]),col='red',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,5]),col='red',lty='dotted')

plot(density(MCMCchains(samples.rjagsWrAg,params="sigma0_g")[,3]))
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma0_pred")[,3]),col='red')

lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,6]),col='black',lty='dotted')

lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,6]),col='red',lty='dotted')


# 
# 
# hist(MCMCchains(samples.rjagsWrAg,params="mu0_g")[,2],breaks=100,freq=FALSE)
# lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,4]),col='red')
# lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,5]),col='red')
# 
# hist(MCMCchains(samples.rjagsWrAg,params="mu0_g")[,3],breaks=100,freq=FALSE)
# lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,6]),col='red')


#sigma0_g=MCMCchains(samples.rjagsWrAg,params="sigma0_g")[,1]

#par(mfrow=c(1,3))
# hist(MCMCchains(samples.rjagsWrAg,params="sigma0_g")[,1],breaks=100,freq=FALSE)
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,1]),col='red')
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,2]),col='red')
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,3]),col='red')
# 
# hist(MCMCchains(samples.rjagsWrAg,params="sigma0_g")[,2],breaks=100,freq=FALSE)
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,4]),col='red')
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,5]),col='red')
# 
# hist(MCMCchains(samples.rjagsWrAg,params="sigma0_g")[,3],breaks=100,freq=FALSE)
# lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,6]),col='red')
# 

## SUMMARY
MCMCsummary(samples.rjagsWrAg,params=c("tau0_g","tau_g"))
MCMCsummary(samples.rjagsWrAg,params=c("sigma0_g","sigma_g"))
MCMCsummary(samples.rjagsWrAg,params=c("mu0_g","mu_g"))


y_pred=MCMCchains(samples.rjagsWrAg,params="y_pred")
y_sim=MCMCchains(samples.rjagsWrAg,params="y_sim")


par(mfrow=c(1,2))
plot(data$seedlingJan,y_pred[sample(dim(y_pred)[1],1),],
     type='n',xlim=c(0,100),ylim=c(0,100),
     xlab="Observed data",ylab="Simulated data")
for(i in 1:100){
  points(data$seedlingJan,y_pred[sample(dim(y_pred)[1],1),],pch=16,cex=.5,
         col=data$yearGermination)
}
abline(a=0,b=1)

plot(data$seedlingJan,y_sim[sample(dim(y_sim)[1],1),],
     type='n',xlim=c(0,100),ylim=c(0,100),
     xlab="Observed data",ylab="Simulated data")
for(i in 1:100){
  points(data$seedlingJan,y_sim[sample(dim(y_sim)[1],1),],pch=16,cex=.5,
         col=data$yearGermination)
}
abline(a=0,b=1)
