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
 # dplyr::filter(site=="KYE"|site=="LO") %>%
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

jmWrAg = jags.model(paste0(dir,"jagsWrAgPriorPredictiveSurvival.R"), data = data, n.adapt = n.adapt,
                    n.chains=2)

# burn-in (n.update)
update(jmWrAg, n.iter = n.update)

samples.rjagsWrAg = coda.samples(jmWrAg,
                                 variable.names = 
                                   # GERMINATION
                                   c("mu0_g","tau0_g","mu_g","tau_g",
                                     "sigma0_g","sigma_g","g",
                                     "mu0_pred","tau0_pred","mu_pred","tau_pred",
                                      "sigma0_pred","sigma_pred","g_pred",
                                     #PRIOR PREDICTIVE
                                     "y_pred","y_prior.pred",
                                     # SURVIVAL
                                     "mu0_s","mu_s","b","theta_c",
                                     "sigma0_s","sigma_s",
                                     "beta","beta_mod","alpha_s","theta_s",
                                     "mu",
                                     # POSTERIOR PREDICTIVE
                                     "y_sim","y_surv"),
                                 n.iter = n.iterations, thin = n.thin)
 
fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/decayModel/")
dir.create(file.path(fileDirectory), showWarnings = FALSE)
# 
samples.rjags=samples.rjagsWrAg
saveRDS(samples.rjags,file=paste0(fileDirectory,"seedBurialSamples.rds"))
saveRDS(data,file=paste0(fileDirectory,"data.rds"))
#saveRDS(seedBagsData,file=paste0(fileDirectory,"seedBagExperiment.rds"))
#saveRDS(viabilityRawData,file=paste0(fileDirectory,"viabilityExperiment.rds"))



  
  theta_c=MCMCchains(samples.rjagsWrAg,params="theta_c")
  par(mfrow=c(4,6))
  for(i in 1:12){
  hist(theta_c[,i],breaks=100,freq=FALSE)
  }
  
  # mu=MCMCchains(samples.rjagsWrAg,params="mu")
  # dim(mu)
  # for(i in 1:6){
  #   hist(mu[,i],breaks=100,freq=FALSE)
  # }
  # 
  # for(i in 1:6){
  #   hist(theta_c[,i]*mu[,i],breaks=100,freq=FALSE)
  # }  
  # 
  # theta_s=MCMCchains(samples.rjagsWrAg,params="theta_s")
  # dim(theta_s)
  # for(i in 1:6){
  #   hist(theta_s[,i],breaks=100,freq=FALSE)
  # }
  
  ind=200
  par(mfrow=c(2,1))
  hist(MCMCchains(samples.rjagsWrAg,params="mu")[,ind],breaks=100)
  hist(MCMCchains(samples.rjagsWrAg,params="theta_s")[,ind],breaks=100)
  
  theta_s.mu=apply(MCMCchains(samples.rjagsWrAg,params="theta_s"),2,median)
  plot(data$y/data$seedStart,theta_s.mu);abline(a=0,b=1)
   
  
  
y_pred.g=MCMCchains(samples.rjagsWrAg,params="y_pred")
y_pred.s=MCMCchains(samples.rjagsWrAg,params="y_prior.pred")
# 
# theta_pred=MCMCchains(samples.rjagsWrAg,params="g_pred")
# theta_g1=MCMCchains(samples.rjagsWrAg,params="theta_g1")
# theta_s=MCMCchains(samples.rjagsWrAg,params="theta_s")
# theta_g3=MCMCchains(samples.rjagsWrAg,params="theta_g3")
# 
# par(mfrow=c(2,3))
# hist(y_pred[sample(dim(y_pred)[1],100),],breaks=100,freq=FALSE)
# 
# hist(theta_pred[,1:dim(y_pred)[2]],breaks=100,freq=FALSE)
# hist(theta_g1[,1:length(data$months)],breaks=100,freq=FALSE)
# hist(theta_s[,3],breaks=100,freq=FALSE)
# hist(theta_g3[,1:length(data$months)],breaks=100,freq=FALSE)
# 
# #hist(MCMCchains(samples.rjagsWrAg,params="b")[,1],breaks=100,freq=FALSE)
# hist(MCMCchains(samples.rjagsWrAg,params="beta")[,1],breaks=100,freq=FALSE)
# hist(MCMCchains(samples.rjagsWrAg,params="beta_mod")[,1],breaks=100,freq=FALSE)

# 
# par(mfrow=c(2,3))
# 
# hist(MCMCchains(samples.rjagsWrAg,params="mu")[,1],breaks=100)
# hist(MCMCchains(samples.rjagsWrAg,params="mu")[,2],breaks=100)
# hist(MCMCchains(samples.rjagsWrAg,params="mu")[,3],breaks=100)
# hist(MCMCchains(samples.rjagsWrAg,params="mu")[,4],breaks=100)
# hist(MCMCchains(samples.rjagsWrAg,params="mu")[,5],breaks=100)
# hist(MCMCchains(samples.rjagsWrAg,params="mu")[,6],breaks=100)
# 
# 
# par(mfrow=c(2,3))
# 
# hist(MCMCchains(samples.rjagsWrAg,params="theta_s")[,1],breaks=100)
# hist(MCMCchains(samples.rjagsWrAg,params="theta_s")[,2],breaks=100)
# hist(MCMCchains(samples.rjagsWrAg,params="theta_s")[,3],breaks=100)
# hist(MCMCchains(samples.rjagsWrAg,params="theta_s")[,4],breaks=100)
# hist(MCMCchains(samples.rjagsWrAg,params="theta_s")[,5],breaks=100)
# hist(MCMCchains(samples.rjagsWrAg,params="theta_s")[,6],breaks=100)


mu0_s=MCMCchains(samples.rjagsWrAg,params="mu0_s")
b0=exp(mu0_s)

mu_s=MCMCchains(samples.rjagsWrAg,params="mu_s")
hist(mu_s[,1],breaks=100,xlim=c(-2,2))
hist(mu_s[,2],breaks=100,xlim=c(-2,2))
hist(mu_s[,3],breaks=100,xlim=c(-2,2))

b=exp(mu_s)
par(mfrow=c(1,1))
plot(density(b[,1]),xlim=c(0,3),ylim=c(0,3))
lines(density(b[,2]),xlim=c(0,3))
lines(density(b[,3]),xlim=c(0,3))

b.parm=apply(b,2,quantile,probs=c(0.025,.5,.975))

t=seq(0,max(data$months),.01)

par(mfrow=c(1,1))
plot(t,exp(-b.parm[2,1]*t),
     xlim=c(0,max(data$months)),ylim=c(0,1),type='n')
colors=c('black',"purple","orange")
for(i in 1:3){
lines(t,exp(-b.parm[2,i]*t),
     xlim=c(0,max(data$months)),ylim=c(0,1),type='l',
     col=colors[i])
}
lines(t,exp(-median(b0)*t),lty='dotted')

par(mfrow=c(1,3))
## AGE 1
plot(t,exp(-b.parm[2,1]*t),
     xlim=c(0,max(data$months)),ylim=c(0,1),type='l')
for(i in 1:100){
  lines(t,exp(-sample(b[,1],1)*t),
        lwd=.5,col='lightgray')
}
lines(t,exp(-b.parm[2,1]*t),lwd=2)
lines(t,exp(-median(b0)*t),lty='dotted')

df<-survivalData[survivalData$yearSurvival==2006,]
points(df$months,df$y/df$seedStart,pch=16,cex=.75)

## AGE 2
plot(t,exp(-b.parm[2,2]*t),
     xlim=c(0,max(data$months)),ylim=c(0,1),type='l')
for(i in 1:100){
  lines(t,exp(-sample(b[,2],1)*t),
        lwd=.5,col='lightgray')
}
lines(t,exp(-b.parm[2,2]*t),lwd=2)
lines(t,exp(-median(b0)*t),lty='dotted')

df<-survivalData[survivalData$yearSurvival==2007,]
points(df$months,df$y/df$seedStart,pch=16,cex=.75)

## AGE 3
plot(t,exp(-b.parm[2,3]*t),
     xlim=c(0,max(data$months)),ylim=c(0,1),type='l')
for(i in 1:100){
  lines(t,exp(-sample(b[,3],1)*t),
        lwd=.5,col='lightgray')
}
lines(t,exp(-b.parm[2,3]*t),lwd=2)
lines(t,exp(-median(b0)*t),lty='dotted')

df<-survivalData[survivalData$yearSurvival==2008,]
points(df$months,df$y/df$seedStart,pch=16,cex=.75)

# 
# 
# plot(t,100*exp(-b.parm[2,]*t),
#      xlim=c(0,1),ylim=c(0,100),type='l')
# polygon(c(t,rev(t)),
#         c(100*exp(-b.parm[1,]*t),rev(100*exp(-b.parm[3,]*t))),
#         col="lightgray")




y_surv = MCMCchains(samples.rjagsWrAg,params='y_surv')


par(mfrow=c(1,1))
plot(rep(data$y[1],100),sample(y_surv[,1],100),type='n',
     xlim=c(0,100),
     ylim=c(0,100))
for(i in 1:dim(y_surv)[2]){
  points(rep(data$y[i],100),sample(y_surv[,i],100),
         pch=16,cex=.5,col='lightgray')
}
abline(a=0,b=1)

plot(apply(y_surv[,data$months==1],2,median),
     data$y[data$months==1],
     type='n',
     xlim=c(0,100),
     ylim=c(0,100))
abline(a=0,b=1)
  points(apply(y_surv[,data$months==1],2,median),
         data$y[data$months==1],
         pch=16,cex=1,col='lightgray')


## GERMINATION
  
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
  
  # Effect of survival data on germination estimate
  
  par(mfrow=c(3,1))
  plot(density(MCMCchains(samples.rjagsWrAg,params="mu0_g")[,1]),xlim=c(-5,5),ylim=c(0,1.25))
  lines(density(MCMCchains(samples.rjagsWrAg,params="mu0_pred")[,1]),col='red')
  
  lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,1]),col='black',lty='dotted')
  lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,2]),col='black',lty='dotted')
  lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,3]),col='black',lty='dotted')
  
  lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,1]),col='red',lty='dotted')
  lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,2]),col='red',lty='dotted')
  lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,3]),col='red',lty='dotted')
  
  plot(density(MCMCchains(samples.rjagsWrAg,params="mu0_g")[,2]),xlim=c(-5,5),ylim=c(0,1.25))
  lines(density(MCMCchains(samples.rjagsWrAg,params="mu0_pred")[,2]),col='red')
  
  lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,4]),col='black',lty='dotted')
  lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,5]),col='black',lty='dotted')
  
  lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,4]),col='red',lty='dotted')
  lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,5]),col='red',lty='dotted')
  
  plot(density(MCMCchains(samples.rjagsWrAg,params="mu0_g")[,3]),xlim=c(-5,5),ylim=c(0,1.25))
  lines(density(MCMCchains(samples.rjagsWrAg,params="mu0_pred")[,3]),col='red')
  
  lines(density(MCMCchains(samples.rjagsWrAg,params="mu_g")[,6]),col='black',lty='dotted')
  
  lines(density(MCMCchains(samples.rjagsWrAg,params="mu_pred")[,6]),col='red',lty='dotted')
  
  
  
  
  
  
y_pred <- MCMCchains(samples.rjagsWrAg,params="y_pred")

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



## SIGMA

plot(density(MCMCchains(samples.rjagsWrAg,params="sigma0_g")[,1]),xlim=c(0,4),ylim=c(0,1.5))
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma0_pred")[,1]),col='red')

lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,1]),col='black',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,2]),col='black',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,3]),col='black',lty='dotted')

lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,1]),col='red',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,2]),col='red',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,3]),col='red',lty='dotted')

plot(density(MCMCchains(samples.rjagsWrAg,params="sigma0_g")[,2]),xlim=c(0,4),ylim=c(0,1.5))
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma0_pred")[,2]),col='red')

lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,4]),col='black',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,5]),col='black',lty='dotted')

lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,4]),col='red',lty='dotted')
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,5]),col='red',lty='dotted')

plot(density(MCMCchains(samples.rjagsWrAg,params="sigma0_g")[,3]),xlim=c(0,4),ylim=c(0,1.5))
lines(density(MCMCchains(samples.rjagsWrAg,params="sigma0_pred")[,3]),col='red')

lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_g")[,6]),col='black',lty='dotted')

lines(density(MCMCchains(samples.rjagsWrAg,params="sigma_pred")[,6]),col='red',lty='dotted')


## SUMMARY
MCMCsummary(samples.rjagsWrAg,params=c("tau0_g","tau_g"))
MCMCsummary(samples.rjagsWrAg,params=c("sigma0_g","sigma_g"))
MCMCsummary(samples.rjagsWrAg,params=c("mu0_g","mu_g"))
