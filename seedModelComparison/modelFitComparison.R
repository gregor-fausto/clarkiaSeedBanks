# to fit distributions with shape parameter add offset
# https://sourceforge.net/p/mcmc-jags/discussion/610037/thread/c21ef62a/

# log likelihood
# https://sourceforge.net/p/mcmc-jags/discussion/610036/thread/8211df61/

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

set.seed(10)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Read data
# -------------------------------------------------------------------
# -------------------------------------------------------------------

data = readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/seedBagMaster.RDS")

data = data %>%
  dplyr::filter(site=="SM") %>%
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
                       months=c(3,12,15,24,27,36))

survivalData=data %>%
  dplyr::select(-seedlingJan) %>%
  tidyr::pivot_longer(cols=c(totalJan,intactOct)) %>%
  dplyr::left_join(refIndex,by=c("name","ageBags"))

refIndex <- data.frame(
                       ageBags=as.factor(c(1,2,3)),
                       gIndex=c(1,2,3))

germinationData=data %>%
  dplyr::select(-c(seedStart,intactOct)) %>%
  dplyr::mutate(gIndex=ageBags)

# octData = data %>%
#   dplyr::select(siteBags,yearBags,ageBags,intactOct) %>%
#   dplyr::mutate(
#     months = ifelse(ageBags==1 , 12, NA),
#     #
#     months = ifelse(ageBags==2 , 24, months),
#     #
#     months = ifelse(ageBags==3 , 36, months)
#   ) %>%
#   dplyr::mutate(seedStart=100)
# 
# data<-data %>%
#   dplyr::select(siteBags,yearBags,ageBags,totalJan,seedlingJan) %>%
#   dplyr::mutate(
#     months = ifelse(ageBags==1 , 3, NA),
#     #
#     months = ifelse(ageBags==2 , 15, months),
#     #
#     months = ifelse(ageBags==3 , 27, months)
#   ) %>%
#   dplyr::mutate(seedStart=100)
# 
# data1 <- data %>%
#   dplyr::filter(ageBags=="1")
# names(data1) = paste(names(data1),"1",sep="")
# 
# octData1 <- octData %>%
#   dplyr::filter(ageBags=="1")
# names(octData1) = paste(names(octData1),"2",sep="")
# 
# data2 <- data %>%
#   dplyr::filter(ageBags=="2")
# names(data2) = paste(names(data2),"3",sep="")
# 
# octData2 <- octData %>%
#   dplyr::filter(ageBags=="2")
# names(octData2) = paste(names(octData2),"4",sep="")
# 
# data3 <- data %>%
#   dplyr::filter(ageBags=="3")
# names(data3) = paste(names(data3),"5",sep="")
# 
# octData3 <- octData %>%
#   dplyr::filter(ageBags=="3")
# names(octData3) = paste(names(octData3),"6",sep="")

data <- tidybayes::compose_data(survivalData,germinationData)

data$n1 = dim(survivalData)[1]
data$n2 = dim(germinationData)[1]

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

dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/seedModelComparison/")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Negative exponential, constant germination
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# set inits for JAGS
inits <- list()

# model with constant age
# # tuning (n.adapt)
jmNeCg = jags.model(paste0(dir,"jagsNeCg.R"), data = data, n.adapt = n.adapt,
                    n.chains=2)
#  inits = inits,
#  n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jmNeCg, n.iter = n.update)

samples.rjagsNeCg = coda.samples(jmNeCg,
                                 variable.names = c("beta","k","g","logLik"),
                                 n.iter = n.iterations, thin = n.thin)

LLmat <- MCMCchains(samples.rjagsNeCg,params="logLik")
rel_n_eff <- relative_eff(exp(LLmat), chain_id = rep(1, each = 2*n.iterations))
looNeCg <- loo(LLmat, r_eff = rel_n_eff, cores = 2)

print(looNeCg)
plot(looNeCg)
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Negative exponential, age-dependent germination
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# Model with variable age
# # tuning (n.adapt)
jmNeAg = jags.model(paste0(dir,"jagsNeAg.R"), data = data, n.adapt = n.adapt)
#  inits = inits,
#  n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jmNeAg, n.iter = n.update)

samples.rjagsNeAg = coda.samples(jmNeAg,
                                 variable.names = c("beta","k","g","logLik"),
                                 n.iter = n.iterations, thin = n.thin)

MCMCsummary(samples.rjagsNeAg)

LLmat <- MCMCchains(samples.rjagsNeAg,params="logLik")
rel_n_eff <- relative_eff(exp(LLmat), chain_id = rep(1, each = n.iterations))
looNeAg <- loo(LLmat, r_eff = rel_n_eff, cores = 2)

print(looNeAg)
plot(looNeAg)
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Continuous exponential, constant germination
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# # tuning (n.adapt)
jmCeCg = jags.model(paste0(dir,"jagsCeCg.R"), data = data, n.adapt = n.adapt)
#  inits = inits,
#  n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jmCeCg, n.iter = n.update)

samples.rjagsCeCg = coda.samples(jmCeCg,
                                 variable.names = c("beta","a","b","g","logLik"),
                                 n.iter = n.iterations, thin = n.thin)

MCMCsummary(samples.rjagsCeCg)

LLmat <- MCMCchains(samples.rjagsCeCg,params="logLik")
rel_n_eff <- relative_eff(exp(LLmat), chain_id = rep(1, each = n.iterations))
looCeCg <- loo(LLmat, r_eff = rel_n_eff, cores = 2)

print(looCeCg)
plot(looCeCg)
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Continuous exponential, age-dependent germination
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# # tuning (n.adapt)
jmCeAg = jags.model(paste0(dir,"jagsCeAg.R"), data = data, n.adapt = n.adapt)
#  inits = inits,
#  n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jmCeAg, n.iter = n.update)

samples.rjagsCeAg = coda.samples(jmCeAg,
                                 variable.names = c("beta","a","b","g","logLik"),
                                 n.iter = n.iterations, thin = n.thin)

MCMCsummary(samples.rjagsCeAg)

LLmat <- MCMCchains(samples.rjagsCeAg,params="logLik")
rel_n_eff <- relative_eff(exp(LLmat), chain_id = rep(1, each = n.iterations))
looCeAg <- loo(LLmat, r_eff = rel_n_eff, cores = 2)

print(looCeAg)
plot(looCeAg)
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Weibull residence time, continuous germination
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# # tuning (n.adapt)
jmWrCg = jags.model(paste0(dir,"jagsWrCg.R"), data = data, n.adapt = n.adapt)
#  inits = inits,
#  n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jmWrCg, n.iter = n.update)

samples.rjagsWrCg = coda.samples(jmWrCg,
                                 variable.names = c("beta","a","b","g"),
                                 n.iter = n.iterations, thin = n.thin)

MCMCsummary(samples.rjagsWrCg)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Weibull residence time, age-dependent germination
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# # tuning (n.adapt)
jmWrAg = jags.model(paste0(dir,"jagsWrAg.R"), data = data, n.adapt = n.adapt)
#  inits = inits,
#  n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jmWrAg, n.iter = n.update)

samples.rjagsWrAg = coda.samples(jmWrAg,
                                 variable.names = c("beta","a","b","g"),
                                 n.iter = n.iterations, thin = n.thin)

MCMCsummary(samples.rjagsWrAg)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Extract parameters
# -------------------------------------------------------------------
# -------------------------------------------------------------------

kNeCg <- mean(MCMCchains(samples.rjagsNeCg,"k"))
kNeAg <- mean(MCMCchains(samples.rjagsNeAg,"k"))

parmsCeCg = c()
parmsCeCg = c(mean(MCMCchains(samples.rjagsCeCg,"a")),
              mean(MCMCchains(samples.rjagsCeCg,"b")))

funCeCg = function(t,par=parmsCeCg){
  a = par[1];
  b = par[2];
  theta = 1/((1+b*t)^a);
  return(theta)
}

parmsCeAg = c()
parmsCeAg = c(mean(MCMCchains(samples.rjagsCeAg,"a")),
              mean(MCMCchains(samples.rjagsCeAg,"b")))

funCeAg = function(t,par=parmsCeAg){
  a = par[1];
  b = par[2];
  theta = 1/((1+b*t)^a);
  return(theta)
}

parmsWrCg = c()
parmsWrCg = c(mean(MCMCchains(samples.rjagsWrCg,"a")),
              mean(MCMCchains(samples.rjagsWrCg,"b")))

funWrCg = function(t,par=parmsWrCg){
  a = par[1];
  b = par[2];
  theta = exp(-(t/b)^a);
  return(theta)
}

parmsWrAg = c()
parmsWrAg = c(mean(MCMCchains(samples.rjagsWrAg,"a")),
              mean(MCMCchains(samples.rjagsWrAg,"b")))

funWrAg = function(t,par=parmsWrAg){
  a = par[1];
  b = par[2];
  theta = exp(-(t/b)^a);
  return(theta)
}


gNeCg <- apply(MCMCchains(samples.rjagsNeCg,"g"),2,mean)
gNeAg <- apply(MCMCchains(samples.rjagsNeAg,"g"),2,mean)
gCeCg <- apply(MCMCchains(samples.rjagsCeCg,"g"),2,mean)
gCeAg <- apply(MCMCchains(samples.rjagsCeAg,"g"),2,mean)
gWrCg <- apply(MCMCchains(samples.rjagsWrCg,"g"),2,mean)
gWrAg <- apply(MCMCchains(samples.rjagsWrAg,"g"),2,mean)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Plotting
# -------------------------------------------------------------------
# -------------------------------------------------------------------
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(4, "PuOr"))(6)

# -------------------------------------------------------------------
# Models
# -------------------------------------------------------------------
par(mfrow=c(1,2))
# plot 1
plot(rep(3,length=data$n1),data$totalJan1/data$seedStart1,
     type='n',
     xlim=c(0,36),ylim=c(0,1),
     xlab="Months",ylab="P(survival)",
     main="Survival model")

lines(seq(0,36,by=0.1),exp(-kNeCg*seq(0,36,by=0.1)),lty='solid',col='red')
lines(seq(0,36,by=0.1),exp(-kNeAg*seq(0,36,by=0.1)),lty='dotted',col='red')
lines(seq(0,36,by=0.1),funCeCg(seq(0,36,by=0.1)),lty='solid',col='purple')
lines(seq(0,36,by=0.1),funCeAg(seq(0,36,by=0.1)),lty='dotted',col='purple')
lines(seq(0,36,by=0.1),funWrCg(seq(0,36,by=0.1)),lty='solid',col='orange')
lines(seq(0,36,by=0.1),funWrAg(seq(0,36,by=0.1)),lty='dotted',col='orange')

# plot 2
plot(rep(1,length=data$n1),data$seedlingJan1/data$totalJan1,
     type='n',
     xlim=c(0,3),ylim=c(0,1),
     xlab="Years",ylab="P(germination)",
     main="Germination model")

abline(h=gNeCg,lty=1,col="red")
abline(h=gCeCg,lty=1,col="purple")
abline(h=gWrCg,lty=1,col="orange")

points(x=c(1,2,3),y=gNeAg,col="red",pch=16,cex=1.5)
points(x=c(1,2,3),y=gCeAg,col="purple",pch=16,cex=1)
points(x=c(1,2,3),y=gWrAg,col="orange",pch=16,cex=.5)


# -------------------------------------------------------------------
# Model vs. data
# -------------------------------------------------------------------
par(mfrow=c(3,2))
# plot 1
plot(rep(3,length=data$n1),data$totalJan1/data$seedStart1,
     type='n',
     xlim=c(0,36),ylim=c(0,1),
     xlab="Months",ylab="P(survival)",
     main="Negative exponential\n Constant germination",cex.main=.75)

points(rep(3,length=data$n1),data$totalJan1/data$seedStart1,pch=16,cex=.5,col='gray')
points(rep(12,length=data$n1),data$intactOct2/data$seedStart2,pch=16,cex=.5,col='gray')
points(rep(15,length=data$n2),data$totalJan3/data$seedStart3,pch=16,cex=.5,col='gray')
points(rep(24,length=data$n2),data$intactOct4/data$seedStart4,pch=16,cex=.5,col='gray')
points(rep(27,length=data$n3),data$totalJan5/data$seedStart5,pch=16,cex=.5,col='gray')
points(rep(36,length=data$n3),data$intactOct6/data$seedStart6,pch=16,cex=.5,col='gray')

points(3,mean(data$totalJan1/data$seedStart1),col=colors[1],pch=16)
points(12,mean(data$intactOct2/data$seedStart2),col=colors[2],pch=16)
points(15,mean(data$totalJan3/data$seedStart3),col=colors[3],pch=16)
points(24,mean(data$intactOct4/data$seedStart4),col=colors[4],pch=16)
points(27,mean(data$totalJan5/data$seedStart5),col=colors[5],pch=16)
points(36,mean(data$intactOct6/data$seedStart6),col=colors[6],pch=16)

lines(seq(0,36,by=0.1),exp(-kNeCg*seq(0,36,by=0.1)),lty='solid',col='black',lwd=2)

points(3,exp(-kNeCg*3),pch=16,col='black')
points(12,(1-gNeCg)*exp(-kNeCg*12),pch=16,col='black')
points(15,(1-gNeCg)*exp(-kNeCg*15),pch=16,col='black')
points(24,(1-gNeCg)*(1-gNeCg)*exp(-kNeCg*24),pch=16,col='black')
points(27,(1-gNeCg)*(1-gNeCg)*exp(-kNeCg*27),pch=16,col='black')
points(36,(1-gNeCg)*(1-gNeCg)*(1-gNeCg)*exp(-kNeCg*36),pch=16,col='black')

# plot 2
plot(rep(3,length=data$n1),data$totalJan1/data$seedStart1,
     type='n',
     xlim=c(0,36),ylim=c(0,1),
     xlab="Months",ylab="P(survival)",
     main="Negative exponential\n Age-dependent germination",cex.main=.75)


points(rep(3,length=data$n1),data$totalJan1/data$seedStart1,pch=16,cex=.5,col='gray')
points(rep(12,length=data$n1),data$intactOct2/data$seedStart2,pch=16,cex=.5,col='gray')
points(rep(15,length=data$n2),data$totalJan3/data$seedStart3,pch=16,cex=.5,col='gray')
points(rep(24,length=data$n2),data$intactOct4/data$seedStart4,pch=16,cex=.5,col='gray')
points(rep(27,length=data$n3),data$totalJan5/data$seedStart5,pch=16,cex=.5,col='gray')
points(rep(36,length=data$n3),data$intactOct6/data$seedStart6,pch=16,cex=.5,col='gray')

points(3,mean(data$totalJan1/data$seedStart1),col=colors[1],pch=16)
points(12,mean(data$intactOct2/data$seedStart2),col=colors[2],pch=16)
points(15,mean(data$totalJan3/data$seedStart3),col=colors[3],pch=16)
points(24,mean(data$intactOct4/data$seedStart4),col=colors[4],pch=16)
points(27,mean(data$totalJan5/data$seedStart5),col=colors[5],pch=16)
points(36,mean(data$intactOct6/data$seedStart6),col=colors[6],pch=16)

lines(seq(0,36,by=0.1),exp(-kNeAg*seq(0,36,by=0.1)),lty='solid',col='black',lwd=2)

points(3,exp(-kNeAg*3),pch=16,col='black')
points(12,(1-gNeAg[1])*exp(-kNeAg*12),pch=16,col='black')
points(15,(1-gNeAg[1])*exp(-kNeAg*15),pch=16,col='black')
points(24,(1-gNeAg[1])*(1-gNeAg[2])*exp(-kNeAg*24),pch=16,col='black')
points(27,(1-gNeAg[1])*(1-gNeAg[2])*exp(-kNeAg*27),pch=16,col='black')
points(36,(1-gNeAg[1])*(1-gNeAg[2])*(1-gNeAg[3])*exp(-kNeAg*36),pch=16,col='black')

# plot 3
plot(rep(3,length=data$n1),data$totalJan1/data$seedStart1,
     type='n',
     xlim=c(0,36),ylim=c(0,1),
     xlab="Months",ylab="P(survival)",
     main="Continuous exponential\n Constant germination",cex.main=.75)

points(rep(3,length=data$n1),data$totalJan1/data$seedStart1,pch=16,cex=.5,col='gray')
points(rep(12,length=data$n1),data$intactOct2/data$seedStart2,pch=16,cex=.5,col='gray')
points(rep(15,length=data$n2),data$totalJan3/data$seedStart3,pch=16,cex=.5,col='gray')
points(rep(24,length=data$n2),data$intactOct4/data$seedStart4,pch=16,cex=.5,col='gray')
points(rep(27,length=data$n3),data$totalJan5/data$seedStart5,pch=16,cex=.5,col='gray')
points(rep(36,length=data$n3),data$intactOct6/data$seedStart6,pch=16,cex=.5,col='gray')

points(3,mean(data$totalJan1/data$seedStart1),col=colors[1],pch=16)
points(12,mean(data$intactOct2/data$seedStart2),col=colors[2],pch=16)
points(15,mean(data$totalJan3/data$seedStart3),col=colors[3],pch=16)
points(24,mean(data$intactOct4/data$seedStart4),col=colors[4],pch=16)
points(27,mean(data$totalJan5/data$seedStart5),col=colors[5],pch=16)
points(36,mean(data$intactOct6/data$seedStart6),col=colors[6],pch=16)

lines(seq(0,36,by=0.1),funCeCg(seq(0,36,by=0.1)),lty='solid',col='black',lwd=2)

points(3,funCeCg(3),pch=16,col='black')
points(12,(1-gNeCg)*funCeCg(12),pch=16,col='black')
points(15,(1-gNeCg)*funCeCg(15),pch=16,col='black')
points(24,(1-gNeCg)*(1-gNeCg)*funCeCg(24),pch=16,col='black')
points(27,(1-gNeCg)*(1-gNeCg)*funCeCg(27),pch=16,col='black')
points(36,(1-gNeCg)*(1-gNeCg)*(1-gNeCg)*funCeCg(36),pch=16,col='black')

# plot 4
plot(rep(3,length=data$n1),data$totalJan1/data$seedStart1,
     type='n',
     xlim=c(0,36),ylim=c(0,1),
     xlab="Months",ylab="P(survival)",
     main="Continuous exponential\n Age-dependent germination",cex.main=.75)

points(rep(3,length=data$n1),data$totalJan1/data$seedStart1,pch=16,cex=.5,col='gray')
points(rep(12,length=data$n1),data$intactOct2/data$seedStart2,pch=16,cex=.5,col='gray')
points(rep(15,length=data$n2),data$totalJan3/data$seedStart3,pch=16,cex=.5,col='gray')
points(rep(24,length=data$n2),data$intactOct4/data$seedStart4,pch=16,cex=.5,col='gray')
points(rep(27,length=data$n3),data$totalJan5/data$seedStart5,pch=16,cex=.5,col='gray')
points(rep(36,length=data$n3),data$intactOct6/data$seedStart6,pch=16,cex=.5,col='gray')

points(3,mean(data$totalJan1/data$seedStart1),col=colors[1],pch=16)
points(12,mean(data$intactOct2/data$seedStart2),col=colors[2],pch=16)
points(15,mean(data$totalJan3/data$seedStart3),col=colors[3],pch=16)
points(24,mean(data$intactOct4/data$seedStart4),col=colors[4],pch=16)
points(27,mean(data$totalJan5/data$seedStart5),col=colors[5],pch=16)
points(36,mean(data$intactOct6/data$seedStart6),col=colors[6],pch=16)

lines(seq(0,36,by=0.1),funCeAg(seq(0,36,by=0.1)),lty='solid',col='black',lwd=2)

points(3,funCeAg(3),pch=16,col='black')
points(12,(1-gCeAg[1])*funCeAg(12),pch=16,col='black')
points(15,(1-gCeAg[1])*funCeAg(15),pch=16,col='black')
points(24,(1-gCeAg[1])*(1-gCeAg[2])*funCeAg(24),pch=16,col='black')
points(27,(1-gCeAg[1])*(1-gCeAg[2])*funCeAg(27),pch=16,col='black')
points(36,(1-gCeAg[1])*(1-gCeAg[2])*(1-gCeAg[3])*funCeAg(36),pch=16,col='black')

# plot 5
plot(rep(3,length=data$n1),data$totalJan1/data$seedStart1,
     type='n',
     xlim=c(0,36),ylim=c(0,1),
     xlab="Months",ylab="P(survival)",
     main="Weibull residence time\n Constant germination",cex.main=.75)

points(rep(3,length=data$n1),data$totalJan1/data$seedStart1,pch=16,cex=.5,col='gray')
points(rep(12,length=data$n1),data$intactOct2/data$seedStart2,pch=16,cex=.5,col='gray')
points(rep(15,length=data$n2),data$totalJan3/data$seedStart3,pch=16,cex=.5,col='gray')
points(rep(24,length=data$n2),data$intactOct4/data$seedStart4,pch=16,cex=.5,col='gray')
points(rep(27,length=data$n3),data$totalJan5/data$seedStart5,pch=16,cex=.5,col='gray')
points(rep(36,length=data$n3),data$intactOct6/data$seedStart6,pch=16,cex=.5,col='gray')

points(3,mean(data$totalJan1/data$seedStart1),col=colors[1],pch=16)
points(12,mean(data$intactOct2/data$seedStart2),col=colors[2],pch=16)
points(15,mean(data$totalJan3/data$seedStart3),col=colors[3],pch=16)
points(24,mean(data$intactOct4/data$seedStart4),col=colors[4],pch=16)
points(27,mean(data$totalJan5/data$seedStart5),col=colors[5],pch=16)
points(36,mean(data$intactOct6/data$seedStart6),col=colors[6],pch=16)

lines(seq(0,36,by=0.1),funWrCg(seq(0,36,by=0.1)),lty='solid',col='black',lwd=2)

points(3,funWrCg(3),pch=16,col='black')
points(12,(1-gWrCg)*funWrCg(12),pch=16,col='black')
points(15,(1-gWrCg)*funWrCg(15),pch=16,col='black')
points(24,(1-gWrCg)*(1-gWrCg)*funWrCg(24),pch=16,col='black')
points(27,(1-gWrCg)*(1-gWrCg)*funWrCg(27),pch=16,col='black')
points(36,(1-gWrCg)*(1-gWrCg)*(1-gWrCg)*funWrCg(36),pch=16,col='black')

# plot 6
plot(rep(3,length=data$n1),data$totalJan1/data$seedStart1,
     type='n',
     xlim=c(0,36),ylim=c(0,1),
     xlab="Months",ylab="P(survival)",
     main="Weibull residence time\n Age-dependent germination",cex.main=.75)

points(rep(3,length=data$n1),data$totalJan1/data$seedStart1,pch=16,cex=.5,col='gray')
points(rep(12,length=data$n1),data$intactOct2/data$seedStart2,pch=16,cex=.5,col='gray')
points(rep(15,length=data$n2),data$totalJan3/data$seedStart3,pch=16,cex=.5,col='gray')
points(rep(24,length=data$n2),data$intactOct4/data$seedStart4,pch=16,cex=.5,col='gray')
points(rep(27,length=data$n3),data$totalJan5/data$seedStart5,pch=16,cex=.5,col='gray')
points(rep(36,length=data$n3),data$intactOct6/data$seedStart6,pch=16,cex=.5,col='gray')

points(3,mean(data$totalJan1/data$seedStart1),col=colors[1],pch=16)
points(12,mean(data$intactOct2/data$seedStart2),col=colors[2],pch=16)
points(15,mean(data$totalJan3/data$seedStart3),col=colors[3],pch=16)
points(24,mean(data$intactOct4/data$seedStart4),col=colors[4],pch=16)
points(27,mean(data$totalJan5/data$seedStart5),col=colors[5],pch=16)
points(36,mean(data$intactOct6/data$seedStart6),col=colors[6],pch=16)

lines(seq(0,36,by=0.1),funWrAg(seq(0,36,by=0.1)),lty='solid',col='black',lwd=2)

points(3,funWrAg(3),pch=16,col='black')
points(12,(1-gWrAg[1])*funWrAg(12),pch=16,col='black')
points(15,(1-gWrAg[1])*funWrAg(15),pch=16,col='black')
points(24,(1-gWrAg[1])*(1-gWrAg[2])*funWrAg(24),pch=16,col='black')
points(27,(1-gWrAg[1])*(1-gWrAg[2])*funWrAg(27),pch=16,col='black')
points(36,(1-gWrAg[1])*(1-gWrAg[2])*(1-gWrAg[3])*funWrAg(36),pch=16,col='black')


# -------------------------------------------------------------------
# Predicted vs. observed
# -------------------------------------------------------------------
par(mfrow=c(3,2))

t=c(3,12,15,24,27,36)

# plot 1

mu.pred.NeCg = c(exp(-kNeCg*3),
                 (1-gNeCg)*exp(-kNeCg*12),
                 (1-gNeCg)*exp(-kNeCg*15),
                 (1-gNeCg)*(1-gNeCg)*exp(-kNeCg*24),
                 (1-gNeCg)*(1-gNeCg)*exp(-kNeCg*27),
                 (1-gNeCg)*(1-gNeCg)*(1-gNeCg)*exp(-kNeCg*36))

mu.obs=c(mean(data$totalJan1/data$seedStart1),
         mean(data$intactOct2/data$seedStart2),
         mean(data$totalJan3/data$seedStart3),
         mean(data$intactOct4/data$seedStart4),
         mean(data$totalJan5/data$seedStart5),
         mean(data$intactOct6/data$seedStart6))

plot(mu.pred.NeCg,mu.obs,pch=16,
     xlim=c(0,1),ylim=c(0,1),
     xlab="Predicted intact",ylab="Observed intact",
     main="Negative exponential\n Constant germination",cex.main=.75)
abline(a=0,b=1)

# plot 2

mu.pred.NeAg = c(exp(-kNeAg*3),
                 (1-gNeAg[1])*exp(-kNeAg*12),
                 (1-gNeAg[1])*exp(-kNeAg*15),
                 (1-gNeAg[1])*(1-gNeAg[2])*exp(-kNeAg*24),
                 (1-gNeAg[1])*(1-gNeAg[2])*exp(-kNeAg*27),
                 (1-gNeAg[1])*(1-gNeAg[2])*(1-gNeAg[3])*exp(-kNeAg*36))

plot(mu.pred.NeAg,mu.obs,pch=16,
     xlim=c(0,1),ylim=c(0,1),
     xlab="Predicted intact",ylab="Observed intact",
     main="Negative exponential\n Age-dependent germination",cex.main=.75)
abline(a=0,b=1)

# plot 3
mu.pred.CeCg = c(funCeCg(3),
                 (1-gCeCg)*funCeCg(12),
                 (1-gCeCg)*funCeCg(15),
                 (1-gCeCg)*(1-gCeCg)*funCeCg(24),
                 (1-gCeCg)*(1-gCeCg)*funCeCg(27),
                 (1-gCeCg)*(1-gCeCg)*(1-gCeCg)*funCeCg(36))

plot(mu.pred.CeCg,mu.obs,pch=16,
     xlim=c(0,1),ylim=c(0,1),
     xlab="Predicted intact",ylab="Observed intact",
     main="Continuous exponential\n Constant germination",cex.main=.75)
abline(a=0,b=1)

# plot 4

mu.pred.CeAg = c(funCeAg(3),
                 (1-gCeAg[1])*funCeAg(12),
                 (1-gCeAg[1])*funCeAg(15),
                 (1-gCeAg[1])*(1-gCeAg[2])*funCeAg(24),
                 (1-gCeAg[1])*(1-gCeAg[2])*funCeAg(27),
                 (1-gCeAg[1])*(1-gCeAg[2])*(1-gCeAg[3])*funCeAg(36))

plot(mu.pred.CeAg,mu.obs,pch=16,
     xlim=c(0,1),ylim=c(0,1),
     xlab="Predicted intact",ylab="Observed intact",
     main="Continuous exponential\n Age-dependent germination",cex.main=.75)
abline(a=0,b=1)

# plot 5
mu.pred.WrCg = c(funWrCg(3),
                 (1-gWrCg)*funWrCg(12),
                 (1-gWrCg)*funWrCg(15),
                 (1-gWrCg)*(1-gWrCg)*funWrCg(24),
                 (1-gWrCg)*(1-gWrCg)*funWrCg(27),
                 (1-gWrCg)*(1-gWrCg)*(1-gWrCg)*funWrCg(36))

plot(mu.pred.WrCg,mu.obs,pch=16,
     xlim=c(0,1),ylim=c(0,1),
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

plot(mu.pred.WrAg,mu.obs,pch=16,
     xlim=c(0,1),ylim=c(0,1),
     xlab="Predicted intact",ylab="Observed intact",
     main="Weibull residence time\n Age-dependent germination",cex.main=.75)
abline(a=0,b=1)