# -------------------------------------------------------------------
# Models for joint estimates of year 1 below ground rates
# Seed survival, germination, and viability
# Models use log-odds parameterization
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

set.seed(10)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize viability trial data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# setwd and read data files
viabilityRawData <- readRDS("~/Dropbox/dataLibrary/postProcessingData/viabilityRawData.rds")

viabilityRawData <- viabilityRawData %>% 
  dplyr::select(-c(germPerc,germNot,viabPerc,viabPerc2,condTest))

viabilityRawData$bag <-as.integer(as.numeric(viabilityRawData$bagNo))

# one row is coded differently so that 
# viabStart=NA and viabStain=NA
# all others have viabStart=NA and viabStain=NA
# recode
viabilityRawData[is.na(viabilityRawData$viabStart),]$viabStart = 0

## data check
viabilityRawData %>% dplyr::filter(germStart - germCount - viabStart<0) 

# filter out rows with problems
# these need to be corrected
viabilityRawData<-viabilityRawData %>% 
  dplyr::filter(germStart - germCount - viabStart >= 0)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# 
# ## filter the dataset for testing purposes
# filterData<-function(x) {
#   x %>%
#     dplyr::filter(age==1) 
# }
# 
# seedBagsData<-filterData(seedBagsData)
# viabilityRawData<-filterData(viabilityRawData)
# 
# # assign variable that combines site and bag; unique id for each bag
# # for each dataset, create that unique identifier again and then
# # use that to link it to the reference identifier created above
# viabilityRawData<-viabilityRawData %>%
#   tidyr::unite(col='id', c(site,bagNo,round,age), sep="", remove=FALSE) %>%
#   tidyr::unite(col='siteBag', c(site,bagNo), sep="", remove=FALSE) %>%
#   dplyr::mutate(siteBag = as.factor(siteBag)) 
# 
# # once each identifier has been created and linked to the reference table
# # and the dataset filtered, the dataset needs to be re-indexed
# # this may be redundant?
# 
# # this line creates a unique id for the subsetted data that is then 
# # used to index each of the 2 datasets
# # and provides the reference set of bags that were included in the experiment
# 
# #https://community.rstudio.com/t/how-to-add-a-counter-to-each-group-in-dplyr/12986/2
# 
# referenceTable<-data.frame(id=union(seedBagsData$id, viabilityRawData$id)) %>%
#   dplyr::mutate(idNo = 1:length(id)) 
# 
# seedBagsData<-seedBagsData %>%
#   dplyr::left_join(referenceTable,by="id")
# 
# viabilityRawData<-viabilityRawData %>%
#   dplyr::left_join(referenceTable,by="id")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# seedBagsData = seedBagsData %>%
#   dplyr::mutate(year = as.factor(yearStart)) %>%
#   dplyr::select(site,year,totalJan,seedStart,seedlingJan,intactOct) %>%
#   dplyr::rename(siteBags = site,
#                 yearBags = year)

viabilityRawData = viabilityRawData %>%
  dplyr::mutate(year = as.factor(round)) %>%
  dplyr::select(site, year, age, germStart, germCount, viabStart, viabStain, bagNo) %>%
  dplyr::rename(siteViab = site,
                yearViab = year,
                ageViab = age,
                bag = bagNo) %>%
  dplyr::mutate(bag = as.factor(bag)) %>%
  dplyr::group_by(siteViab,yearViab, ageViab, bag) %>%
  # sum observations in each bag; this is ignoring some variation
  dplyr::summarise(germStart = sum(germStart),
                   germCount = sum(germCount),
                   viabStart = sum(viabStart),
                   viabStain = sum(viabStain))

data <- tidybayes::compose_data(viabilityRawData)

#data$n = dim(seedBagsData)[1]
# data$nViab = dim(viabilityExperiment)[1]
# data$n <- NULL

# pass data to list for JAGS
# data = list(
# 
#   # Seed burial experiment, year one
#   y = as.double(seedBagExperiment$totalJan),
#   n = as.double(seedBagExperiment$seedStart),
#   
#   N = nrow(seedBagExperiment)
#   
# )

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 3000
n.update = 5000
n.iterations = 10000
n.thin = 1

set.seed(10)
dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/priorChecking/jagsScripts/")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Complete pooling of germination and viability trials
# Partial pooling of seed burial experiment (site level)
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# set inits for JAGS
# inits = list(list(mu0_1 = rep(0,data$n_siteBags), sigma0_1 = rep(.5,data$n_siteBags),
#                   sigma_1 = matrix(rep(.5,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
#                   mu0_2 = rep(0,data$n_siteBags), sigma0_2 = rep(.5,data$n_siteBags),
#                   sigma_2 = matrix(rep(.5,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
#                   mu0_3 = rep(0,data$n_siteBags), sigma0_3 = rep(.5,data$n_siteBags),
#                   sigma_3 = matrix(rep(.5,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
#                   mu0_g = rep(0,data$n_siteBags), sigma0_g = rep(.5,data$n_siteBags),
#                   sigma_g = matrix(rep(.5,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
#                   mu0_v = rep(0,data$n_siteBags), sigma0_v = rep(.5,data$n_siteBags),
#                   sigma_v = matrix(rep(.5,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags)),
#              list(mu0_1 = rep(-1,data$n_siteBags), sigma0_1 = rep(1,data$n_siteBags),
#                   sigma_1 = matrix(rep(1,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
#                   mu0_2 = rep(-1,data$n_siteBags), sigma0_2 = rep(1,data$n_siteBags),
#                   sigma_2 = matrix(rep(1,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
#                   mu0_3 = rep(-1,data$n_siteBags), sigma0_3 = rep(1,data$n_siteBags),
#                   sigma_3 = matrix(rep(1,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
#                   mu0_g = rep(-1,data$n_siteBags), sigma0_g = rep(1,data$n_siteBags),
#                   sigma_g = matrix(rep(1,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
#                   mu0_v = rep(-1,data$n_siteBags), sigma0_v = rep(1,data$n_siteBags),
#                   sigma_v = matrix(rep(1,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags)),
#              list(mu0_1 = rep(1,data$n_siteBags), sigma0_1 = rep(1.25,data$n_siteBags),
#                   sigma_1 = matrix(rep(1.25,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
#                   mu0_2 = rep(1,data$n_siteBags), sigma0_2 = rep(1.25,data$n_siteBags),
#                   sigma_2 = matrix(rep(1.25,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
#                   mu0_3 = rep(1,data$n_siteBags), sigma0_3 = rep(1.25,data$n_siteBags),
#                   sigma_3 = matrix(rep(1.25,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
#                   mu0_g = rep(1,data$n_siteBags), sigma0_g = rep(1.25,data$n_siteBags),
#                   sigma_g = matrix(rep(1.25,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags),
#                   mu0_v = rep(1,data$n_siteBags), sigma0_v = rep(1.25,data$n_siteBags),
#                   sigma_v = matrix(rep(1.25,data$n_siteBags*data$n_yearBags),nrow=data$n_siteBags,ncol=data$n_yearBags)))

# # Call to JAGS
# 
# # tuning (n.adapt)
jm = jags.model(paste0(dir,"jagsViabilityTrials.R"), data = data,
                n.chains = 3, n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

# parsToMonitor_1 = c("theta_1","mu0_1","sigma0_1","mu_1","sigma_1","p_1")
# parsToMonitor_2 = c("theta_2","mu0_2","sigma0_2","mu_2","sigma_2","p_2")
# parsToMonitor_3 = c("theta_3","mu0_3","sigma0_3","mu_3","sigma_3","p_3")
 parsToMonitor_g = c("theta_g","mu0_g","sigma0_g","mu_g","sigma_g","p_g")
 parsToMonitor_v = c("theta_v","mu0_v","sigma0_v","mu_v","sigma_v","p_v")
# parsToMonitor_deriv = c("nu_1","s1","g1","s2")
# parsToMonitor_deriv2 = c("nu0_1","s1.0","g1.0","s2.0")


# chain (n.iter)
samples.rjags = coda.samples(jm, 
                             variable.names = c(parsToMonitor_g,parsToMonitor_v),
                             n.iter = n.iterations, thin = n.thin)

fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBurial/")
dir.create(file.path(fileDirectory), showWarnings = FALSE)
#
saveRDS(samples.rjags,file=paste0(fileDirectory,"viabilityTrialSamples.rds"))
saveRDS(data,file=paste0(fileDirectory,"viabilityData.rds"))
#saveRDS(seedBagsData,file=paste0(fileDirectory,"seedBagExperiment.rds"))
#saveRDS(viabilityRawData,file=paste0(fileDirectory,"viabilityExperiment.rds"))


MCMCsummary(samples.rjags, params = c("mu0_g","mu0_v"))

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBurial/"
modelFittingFiles <- paste0(directory,list.files(directory))

samples.rjags <- readRDS(modelFittingFiles[[6]])
data <- readRDS(modelFittingFiles[[4]])


##
g_pop=apply(MCMCchains(samples.rjags, params = "mu_g"),2,boot::inv.logit)
v_pop=apply(MCMCchains(samples.rjags, params = "mu_v"),2,boot::inv.logit)

nu_pop = g_pop + v_pop*(1-g_pop)
nu1.sum = apply(nu_pop,2,quantile,c(.025,.5,.975))

#plots summarizing germination and viability experiments
g.sum = apply(g_pop,2,quantile,c(.025,.5,.975))
v.sum = apply(v_pop,2,quantile,c(.025,.5,.975))

siteName = unique(viabilityRawData$siteViab)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  plot(NA,NA,type='n',xlim=c(.5,3.5),ylim=c(0,1),ylab='',xlab='',xaxt='n',yaxt='n')
    tmp.g = g.sum[,c(i,i+60,i+100)]
    tmp.v = v.sum[,c(i,i+60,i+100)]
    tmp.nu = nu.sum[,c(i,i+60,i+100)]
    
   # segments(y0=tmp.g[1,],y1=tmp.g[3,],x0=c(1,2,3),col='gray')
   # points(x=c(1,2,3),y=tmp.g[2,],pch=19,cex=.6,col='gray')
    lines(x=c(1,2,3),y=tmp.g[2,],col='gray')
    
   # segments(y0=tmp.v[1,],y1=tmp.v[3,],x0=c(1,2,3),col='orange')
   # points(x=c(1,2,3),y=tmp.v[2,],pch=19,cex=.6,col='orange')
    lines(x=c(1,2,3),y=tmp.v[2,],col='orange')
    
    segments(y0=tmp.nu[1,],y1=tmp.nu[3,],x0=c(1,2,3))
    points(x=c(1,2,3),y=tmp.nu[2,],pch=19,cex=.6)
    lines(x=c(1,2,3),y=tmp.nu[2,],col='black')
    
    text(x=1,y=.05,siteName[i])
    ifelse(i%in%c(16:20),axis(1L),NA)
    ifelse(i%in%c(1,6,11,16),axis(2L),NA)
    ifelse(i==1,
           legend(x=2,y=0.25,
                  c("Germination","Viability","Total"),
                  lty=c(1,1,1),col=c('gray','orange','black'),
                  cex=.5,box.lty=0),
           NA)
}
mtext("Age (years)", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of viability", side = 2, outer = TRUE, line = 2.2)
mtext("Round 1", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  plot(NA,NA,type='n',xlim=c(.5,3.5),ylim=c(0,1),ylab='',xlab='',xaxt='n',yaxt='n')
  tmp.g = g.sum[,c(i+20,i+80)]
  tmp.v = v.sum[,c(i+20,i+80)]
  tmp.nu = nu.sum[,c(i+20,i+80)]
  
  # segments(y0=tmp.g[1,],y1=tmp.g[3,],x0=c(1,2,3),col='gray')
  # points(x=c(1,2,3),y=tmp.g[2,],pch=19,cex=.6,col='gray')
  lines(x=c(1,2),y=tmp.g[2,],col='gray')
  
  # segments(y0=tmp.v[1,],y1=tmp.v[3,],x0=c(1,2,3),col='orange')
  # points(x=c(1,2,3),y=tmp.v[2,],pch=19,cex=.6,col='orange')
  lines(x=c(1,2),y=tmp.v[2,],col='orange')
  
  segments(y0=tmp.nu[1,],y1=tmp.nu[3,],x0=c(1,2))
  points(x=c(1,2),y=tmp.nu[2,],pch=19,cex=.6)
  lines(x=c(1,2),y=tmp.nu[2,],col='black')
  
  text(x=1,y=.05,siteName[i])
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
  ifelse(i==1,
         legend(x=2,y=0.25,
                c("Germination","Viability","Total"),
                lty=c(1,1,1),col=c('gray','orange','black'),
                cex=.5,box.lty=0),
         NA)
}
mtext("Age (years)", side = 1, outer = TRUE, line = 2.2)
mtext("Probability", side = 2, outer = TRUE, line = 2.2)
mtext("Round 2", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)


par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  plot(NA,NA,type='n',xlim=c(.5,3.5),ylim=c(0,1),ylab='',xlab='',xaxt='n',yaxt='n')
  tmp.g = as.matrix(g.sum[,c(i+40)])
  tmp.v = as.matrix(v.sum[,c(i+40)])
  tmp.nu = as.matrix(nu.sum[,c(i+40)])
  
  # segments(y0=tmp.g[1,],y1=tmp.g[3,],x0=c(1,2,3),col='gray')
   points(x=c(1),y=tmp.g[2,],pch=19,cex=.6,col='gray')
  lines(x=c(1),y=tmp.g[2,],col='gray')
  
  # segments(y0=tmp.v[1,],y1=tmp.v[3,],x0=c(1,2,3),col='orange')
   points(x=c(1),y=tmp.v[2,],pch=19,cex=.6,col='orange')
  lines(x=c(1),y=tmp.v[2,],col='orange')
  
  segments(y0=tmp.nu[1,],y1=tmp.nu[3,],x0=c(1))
  points(x=c(1),y=tmp.nu[2,],pch=19,cex=.6)
  lines(x=c(1),y=tmp.nu[2,],col='black')
  
  text(x=1,y=.05,siteName[i])
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
  ifelse(i==1,
         legend(x=2,y=0.25,
                c("Germination","Viability","Total"),
                lty=c(1,1,1),col=c('gray','orange','black'),
                cex=.5,box.lty=0),
         NA)
}
mtext("Age (years)", side = 1, outer = TRUE, line = 2.2)
mtext("Probability", side = 2, outer = TRUE, line = 2.2)
mtext("Round 3", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)



plot(g.sum[2,],v.sum[2,],xlim=c(0,1),ylim=c(0,1))

par(mfrow=c(3,3))
plot(g.sum[2,1:60],v.sum[2,1:60],xlim=c(0,1),ylim=c(0,1))
plot(g.sum[2,61:100],v.sum[2,61:100],xlim=c(0,1),ylim=c(0,1))
plot(g.sum[2,101:120],v.sum[2,101:120],xlim=c(0,1),ylim=c(0,1))

plot(g.sum[2,c(1:20,61:80,101:120)],v.sum[2,c(1:20,61:80,101:120)],xlim=c(0,1),ylim=c(0,1))
plot(g.sum[2,c(21:40,81:100)],v.sum[2,c(21:40,81:100)],xlim=c(0,1),ylim=c(0,1))
plot(g.sum[2,41:60],v.sum[2,41:60],xlim=c(0,1),ylim=c(0,1))

plot(nu.sum[2,c(1:20,61:80,101:120)],g.sum[2,c(1:20,61:80,101:120)],xlim=c(0,1),ylim=c(0,1))
plot(nu.sum[2,c(21:40,81:100)],g.sum[2,c(21:40,81:100)],xlim=c(0,1),ylim=c(0,1))
plot(nu.sum[2,41:60],g.sum[2,41:60],xlim=c(0,1),ylim=c(0,1))

plot(nu.sum[2,c(1:20,61:80,101:120)],v.sum[2,c(1:20,61:80,101:120)],xlim=c(0,1),ylim=c(0,1))
plot(nu.sum[2,c(21:40,81:100)],v.sum[2,c(21:40,81:100)],xlim=c(0,1),ylim=c(0,1))
plot(nu.sum[2,41:60],v.sum[2,41:60],xlim=c(0,1),ylim=c(0,1))

plot(g.sum[2,c(1:20,61:80,101:120)],nu.sum[2,c(1:20,61:80,101:120)],xlim=c(0,1),ylim=c(0,1))
plot(g.sum[2,c(21:40,81:100)],nu.sum[2,c(21:40,81:100)],xlim=c(0,1),ylim=c(0,1))
plot(g.sum[2,41:60],nu.sum[2,41:60],xlim=c(0,1),ylim=c(0,1))

plot(v.sum[2,c(1:20,61:80,101:120)],nu.sum[2,c(1:20,61:80,101:120)],xlim=c(0,1),ylim=c(0,1))
plot(v.sum[2,c(21:40,81:100)],nu.sum[2,c(21:40,81:100)],xlim=c(0,1),ylim=c(0,1))
plot(v.sum[2,41:60],nu.sum[2,41:60],xlim=c(0,1),ylim=c(0,1))

### START YEAR BY YEAR PLOT


#par(mfrow=c(1,2))
mu_g<-MCMCchains(samples.rjags, params = "mu_g")
mu_g.inv <- apply(mu_g,2,boot::inv.logit)

mu_v<-MCMCchains(samples.rjags, params = "mu_v")
mu_v.inv <- apply(mu_v,2,boot::inv.logit)

nu=exp(mu_g)/(1+exp(mu_g))+(exp(-mu_g+mu_v)/(1+exp(-mu_g)+exp(mu_v)+exp(-mu_g+mu_v)))

# 
# hist(nu[,1],breaks=25,xlim=c(0,1)) 
# hist(nu_pop[,1],breaks=25,add=TRUE,col='red') 

nu_year1=nu[,c(1:20,61:80,101:120)]

nu_ratio1=(nu_year1[,c(1:20)])^(1/3)
nu_ratio_sum1 = apply(nu_ratio1,2,quantile,c(.025,.5,.975))

nu_ratio2=nu_year1[,c(1:20)]*(nu_year1[,c(21:40)]/nu_year1[,c(1:20)])^(1/3)
nu_ratio_sum2 = apply(nu_ratio2,2,quantile,c(.025,.5,.975))

nu_ratio3=nu_year1[,c(21:40)]*(nu_year1[,c(41:60)]/nu_year1[,c(21:40)])^(1/3)
nu_ratio_sum3 = apply(nu_ratio3,2,quantile,c(.025,.5,.975))

g_popyear=apply(MCMCchains(samples.rjags, params = "mu_g"),2,boot::inv.logit)
v_popyear=apply(MCMCchains(samples.rjags, params = "mu_v"),2,boot::inv.logit)


g_popyear=g_popyear[,c(1:20,61:80,101:120)]
v_popyear=v_popyear[,c(1:20,61:80,101:120)]

nu_pop = g_popyear + v_popyear*(1-g_popyear)
nu1.sum = apply(nu_pop,2,quantile,c(.025,.5,.975))

pdf("~/Dropbox/clarkiaSeedBanks/products/figures/viability-estimates.pdf",width=8,height=6)

# AGE 1
siteName = unique(viabilityRawData$siteViab)
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  plot(NA,NA,type='n',xlim=c(0,3.5),ylim=c(0,1),ylab='',xlab='',xaxt='n',yaxt='n')
  tmp = nu1.sum[,c(i,i+20,i+40)]
  segments(y0=tmp[1,],y1=tmp[3,],x0=c(1,2,3))
  points(x=c(1,2,3),y=tmp[2,],pch=19)
  
  # x0=seq(0,1,by=.01)
  # lines(x=x0,tmp[2,1]^(x0),lty='dotted')
  # 
  # x1=seq(1,2,by=.01)
  # lines(x=x1,tmp[2,1]*(tmp[2,2]/tmp[2,1])^(x1-1),lty='dotted')
  # 
  # x2=seq(2,3,by=.01)
  # lines(x=x2,tmp[2,2]*(tmp[2,3]/tmp[2,2])^(x2-2),lty='dotted')
  # 
  # points(x=1/3,tmp[2,1]^(1/3),pch=21)
  # points(x=1+1/3,tmp[2,1]*(tmp[2,2]/tmp[2,1])^(1/3),pch=21,bg='white')
 
   segments(y0=nu_ratio_sum1[1,i],y1=nu_ratio_sum1[3,i],x0=c(1/3),col='black',lty='dotted')
  points(x=1/3,(nu_ratio_sum1[2,i]),pch=21,bg='white')
    
  segments(y0=nu_ratio_sum2[1,i],y1=nu_ratio_sum2[3,i],x0=c(1+1/3),col='black',lty='dotted')
  points(x=1+1/3,(nu_ratio_sum2[2,i]),pch=21,bg='white')
  
  segments(y0=nu_ratio_sum3[1,i],y1=nu_ratio_sum3[3,i],x0=c(2+1/3),col='black',lty='dotted')
  points(x=2+1/3,(nu_ratio_sum3[2,i]),pch=21,bg='white')
  # points(x=2+1/3,tmp[2,2]*(tmp[2,3]/tmp[2,2])^(1/3),pch=21,bg='white')

  text(x=.5,y=.05,siteName[i])
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Age (years)", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of viability", side = 2, outer = TRUE, line = 2.2)
mtext("Round 1", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)

#dev.off()

mu_g<-MCMCchains(samples.rjags, params = "mu_g")
mu_g.inv <- apply(mu_g,2,boot::inv.logit)

mu_v<-MCMCchains(samples.rjags, params = "mu_v")
mu_v.inv <- apply(mu_v,2,boot::inv.logit)

nu=exp(mu_g)/(1+exp(mu_g))+(exp(-mu_g+mu_v)/(1+exp(-mu_g)+exp(mu_v)+exp(-mu_g+mu_v)))

nu_year2=nu[,c(21:40,81:100)]

nu_ratio4=(nu_year2[,c(1:20)])^(1/3)
nu_ratio_sum4 = apply(nu_ratio4,2,quantile,c(.025,.5,.975))

nu_ratio5=nu_year2[,c(1:20)]*(nu_year2[,c(21:40)]/nu_year2[,c(1:20)])^(1/3)
nu_ratio_sum5 = apply(nu_ratio5,2,quantile,c(.025,.5,.975))

g_popyear=apply(MCMCchains(samples.rjags, params = "mu_g"),2,boot::inv.logit)
v_popyear=apply(MCMCchains(samples.rjags, params = "mu_v"),2,boot::inv.logit)

g_popyear=g_popyear[,c(21:40,81:100)]
v_popyear=v_popyear[,c(21:40,81:100)]

nu_pop = g_popyear + v_popyear*(1-g_popyear)
nu2.sum = apply(nu_pop,2,quantile,c(.025,.5,.975))


siteName = unique(viabilityRawData$siteViab)
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  plot(NA,NA,type='n',xlim=c(0,3.5),ylim=c(0,1),ylab='',xlab='',xaxt='n',yaxt='n')
  tmp = nu2.sum[,c(i,i+20)]
  segments(y0=tmp[1,],y1=tmp[3,],x0=c(1,2))
  points(x=c(1,2),y=tmp[2,],pch=19)
  
  # x0=seq(0,1,by=.01)
  # lines(x=x0,tmp[2,1]^(x0),lty='dotted')
  # 
  # x1=seq(1,2,by=.01)
  # lines(x=x1,tmp[2,1]*(tmp[2,2]/tmp[2,1])^(x1-1),lty='dotted')
  # 
  # points(x=1/3,tmp[2,1]^(1/3),pch=21)
  # points(x=1+1/3,tmp[2,1]*(tmp[2,2]/tmp[2,1])^(1/3),pch=21,bg='white')

  segments(y0=nu_ratio_sum4[1,i],y1=nu_ratio_sum4[3,i],x0=c(1/3),col='black',lty='dotted')
  points(x=1/3,(nu_ratio_sum4[2,i]),pch=21,bg='white')
  
  segments(y0=nu_ratio_sum5[1,i],y1=nu_ratio_sum5[3,i],x0=c(1+1/3),col='black',lty='dotted')
  points(x=1+1/3,(nu_ratio_sum5[2,i]),pch=21,bg='white')
  
  text(x=.5,y=.05,siteName[i])
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
  #ifelse(i==1,mtext("Round 2", side = 3, line = 0),NA)
}
mtext("Age (years)", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of viability", side = 2, outer = TRUE, line = 2.2)
mtext("Round 2", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)


mu_g<-MCMCchains(samples.rjags, params = "mu_g")
mu_g.inv <- apply(mu_g,2,boot::inv.logit)

mu_v<-MCMCchains(samples.rjags, params = "mu_v")
mu_v.inv <- apply(mu_v,2,boot::inv.logit)

nu=exp(mu_g)/(1+exp(mu_g))+(exp(-mu_g+mu_v)/(1+exp(-mu_g)+exp(mu_v)+exp(-mu_g+mu_v)))

nu_year3=nu[,c(41:60)]

nu_ratio6=(nu_year3[,c(1:20)])^(1/3)
nu_ratio_sum6 = apply(nu_ratio6,2,quantile,c(.025,.5,.975))

g_popyear=apply(MCMCchains(samples.rjags, params = "mu_g"),2,boot::inv.logit)
v_popyear=apply(MCMCchains(samples.rjags, params = "mu_v"),2,boot::inv.logit)

g_popyear=g_popyear[,c(41:60)]
v_popyear=v_popyear[,c(41:60)]

nu_pop = g_popyear + v_popyear*(1-g_popyear)
nu3.sum = apply(nu_pop,2,quantile,c(.025,.5,.975))

siteName = unique(viabilityRawData$siteViab)
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  plot(NA,NA,type='n',xlim=c(0,3.5),ylim=c(0,1),ylab='',xlab='',xaxt='n',yaxt='n')
  tmp = as.matrix(nu3.sum[,c(i)])
  segments(y0=tmp[1,],y1=tmp[3,],x0=c(1))
  points(x=c(1),y=tmp[2,],pch=19)
  
  # x0=seq(0,1,by=.01)
  # lines(x=x0,tmp[2,1]^(x0),lty='dotted')
  # 
  # points(x=1/3,tmp[2,1]^(1/3),pch=21)
  
  segments(y0=nu_ratio_sum6[1,i],y1=nu_ratio_sum6[3,i],x0=c(1/3),col='black',lty='dotted')
  points(x=1/3,(nu_ratio_sum6[2,i]),pch=21,bg='white')
  
  text(x=.5,y=.05,siteName[i])
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Age (years)", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of viability", side = 2, outer = TRUE, line = 2.2)
mtext("Round 3", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)

# COMBINE ALL AGE SPECIFIC ESTIMATES
# AGE 1
siteName = unique(viabilityRawData$siteViab)
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  plot(NA,NA,type='n',xlim=c(0,3.5),ylim=c(0,1),ylab='',xlab='',xaxt='n',yaxt='n')
  tmp = nu1.sum[,c(i,i+20,i+40)]
  segments(y0=tmp[1,],y1=tmp[3,],x0=c(c(1,2)-.1,3))
  points(x=c(c(1,2)-.1,3),y=tmp[2,],pch=19,cex=.6)
  
  segments(y0=nu_ratio_sum1[1,i],y1=nu_ratio_sum1[3,i],x0=c(1/3)-.1,col='black',lty='dotted')
  points(x=1/3-.1,(nu_ratio_sum1[2,i]),pch=21,bg='white',cex=.6)
  
  segments(y0=nu_ratio_sum2[1,i],y1=nu_ratio_sum2[3,i],x0=c(1+1/3)-.1,col='black',lty='dotted')
  points(x=1+1/3-.1,(nu_ratio_sum2[2,i]),pch=21,bg='white',cex=.6)
  
  segments(y0=nu_ratio_sum3[1,i],y1=nu_ratio_sum3[3,i],x0=c(2+1/3),col='black',lty='dotted')
  points(x=2+1/3,(nu_ratio_sum3[2,i]),pch=21,bg='white',cex=.6)

  # round 2
  tmp = nu2.sum[,c(i,i+20)]
  segments(y0=tmp[1,],y1=tmp[3,],x0=c(1,2))
  points(x=c(1,2),y=tmp[2,],pch=19,cex=.6)
  
  segments(y0=nu_ratio_sum4[1,i],y1=nu_ratio_sum4[3,i],x0=c(1/3),col='black',lty='dotted')
  points(x=1/3,(nu_ratio_sum4[2,i]),pch=21,bg='white',cex=.6)
  
  segments(y0=nu_ratio_sum5[1,i],y1=nu_ratio_sum5[3,i],x0=c(1+1/3),col='black',lty='dotted')
  points(x=1+1/3,(nu_ratio_sum5[2,i]),pch=21,bg='white',cex=.6)
  
  # round 3
  tmp = as.matrix(nu3.sum[,c(i)])
  segments(y0=tmp[1,],y1=tmp[3,],x0=c(1)+.1)
  points(x=c(1)+.1,y=tmp[2,],pch=19, cex=.6)

  segments(y0=nu_ratio_sum6[1,i],y1=nu_ratio_sum6[3,i],x0=c(1/3)+.1,col='black',lty='dotted')
  points(x=1/3+.1,(nu_ratio_sum6[2,i]),pch=21,bg='white', cex=.6)
  
  
  text(x=.5,y=.05,siteName[i])
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Age (years)", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of viability", side = 2, outer = TRUE, line = 2.2)
mtext("Round 1-3", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)


dev.off()


## POPULATION PLOT



#par(mfrow=c(1,2))
mu_g<-MCMCchains(samples.rjags, params = "mu_g")
mu_g.inv <- apply(mu_g,2,boot::inv.logit)

mu_v<-MCMCchains(samples.rjags, params = "mu_v")
mu_v.inv <- apply(mu_v,2,boot::inv.logit)

nu=exp(mu_g)/(1+exp(mu_g))+(exp(-mu_g+mu_v)/(1+exp(-mu_g)+exp(mu_v)+exp(-mu_g+mu_v)))


mu0_g<-MCMCchains(samples.rjags, params = "mu0_g")
mu0_g.inv <- apply(mu0_g,2,boot::inv.logit)

mu0_v<-MCMCchains(samples.rjags, params = "mu0_v")
mu0_v.inv <- apply(mu0_v,2,boot::inv.logit)

nu0=exp(mu0_g)/(1+exp(mu0_g))+(exp(-mu0_g+mu0_v)/(1+exp(-mu0_g)+exp(mu0_v)+exp(-mu0_g+mu0_v)))

nu0_ratio1=(nu0[,c(1:20)])^(1/3)
nu0_ratio_sum1 = apply(nu0_ratio1,2,quantile,c(.025,.5,.975))

nu0_ratio2=nu0[,c(1:20)]*(nu0[,c(21:40)]/nu0[,c(1:20)])^(1/3)
nu0_ratio_sum2 = apply(nu0_ratio2,2,quantile,c(.025,.5,.975))

nu0_ratio3=nu0[,c(21:40)]*(nu[,c(101:120)]/nu0[,c(21:40)])^(1/3)
nu0_ratio_sum3 = apply(nu0_ratio3,2,quantile,c(.025,.5,.975))

g_popyear=apply(MCMCchains(samples.rjags, params = "mu_g"),2,boot::inv.logit)
v_popyear=apply(MCMCchains(samples.rjags, params = "mu_v"),2,boot::inv.logit)

nu_pop = g_popyear + v_popyear*(1-g_popyear)
nu.sum = apply(nu_pop,2,quantile,c(.025,.5,.975))

g_pop=apply(MCMCchains(samples.rjags, params = "mu0_g"),2,boot::inv.logit)
v_pop=apply(MCMCchains(samples.rjags, params = "mu0_v"),2,boot::inv.logit)

nu0_pop = g_pop + v_pop*(1-g_pop)
nu0.sum = apply(nu0_pop,2,quantile,c(.025,.5,.975))

pdf("~/Dropbox/clarkiaSeedBanks/products/figures/viability-estimates-population.pdf",width=8,height=6)

# AGE 1
siteName = unique(viabilityRawData$siteViab)
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  plot(NA,NA,type='n',xlim=c(0,3.5),ylim=c(0,1),ylab='',xlab='',xaxt='n',yaxt='n')
  tmp = nu0.sum[,c(i,i+20)]
  segments(y0=tmp[1,],y1=tmp[3,],x0=c(1,2))
  points(x=c(1,2),y=tmp[2,],pch=19)
  
  tmp2 = as.matrix(nu.sum[,c(i+100)])
  segments(y0=tmp2[1,],y1=tmp2[3,],x0=c(3))
  points(x=c(3),y=tmp2[2,],pch=19)
  
  segments(y0=nu0_ratio_sum1[1,i],y1=nu0_ratio_sum1[3,i],x0=c(1/3),col='black',lty='dotted')
  points(x=1/3,(nu0_ratio_sum1[2,i]),pch=21,bg='white')
  
  segments(y0=nu0_ratio_sum2[1,i],y1=nu0_ratio_sum2[3,i],x0=c(1+1/3),col='black',lty='dotted')
  points(x=1+1/3,(nu0_ratio_sum2[2,i]),pch=21,bg='white')
  
  segments(y0=nu0_ratio_sum3[1,i],y1=nu0_ratio_sum3[3,i],x0=c(2+1/3),col='black',lty='dotted')
  points(x=2+1/3,(nu0_ratio_sum3[2,i]),pch=21,bg='white')
  # points(x=2+1/3,tmp[2,2]*(tmp[2,3]/tmp[2,2])^(1/3),pch=21,bg='white')
  
  text(x=.5,y=.05,siteName[i])
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Age (years)", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of viability", side = 2, outer = TRUE, line = 2.2)
mtext("Population-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)






# COMBINE ALL AGE SPECIFIC ESTIMATES
# AGE 1
siteName = unique(viabilityRawData$siteViab)
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  plot(NA,NA,type='n',xlim=c(0,3.5),ylim=c(0,1),ylab='',xlab='',xaxt='n',yaxt='n')
  
  tmp = nu0.sum[,c(i,i+20)]
  segments(y0=tmp[1,],y1=tmp[3,],x0=c(1,2),col='gray')
  points(x=c(1,2),y=tmp[2,],pch=19,col='gray')
  
  tmp2 = as.matrix(nu.sum[,c(i+100)])
  segments(y0=tmp2[1,],y1=tmp2[3,],x0=c(3),col='gray')
  points(x=c(3),y=tmp2[2,],pch=19,col='gray')
  
  segments(y0=nu0_ratio_sum1[1,i],y1=nu0_ratio_sum1[3,i],x0=c(1/3),col='gray',lty='dotted')
  points(x=1/3,(nu0_ratio_sum1[2,i]),pch=21,bg='white',col='gray')
  
  segments(y0=nu0_ratio_sum2[1,i],y1=nu0_ratio_sum2[3,i],x0=c(1+1/3),col='gray',lty='dotted')
  points(x=1+1/3,(nu0_ratio_sum2[2,i]),pch=21,bg='white',col='gray')
  
  segments(y0=nu0_ratio_sum3[1,i],y1=nu0_ratio_sum3[3,i],x0=c(2+1/3),col='black',lty='dotted')
  points(x=2+1/3,(nu0_ratio_sum3[2,i]),pch=21,bg='white',col='gray')

  
  
  
  # YEARLY ESTIMATES
  
  tmp = nu1.sum[,c(i,i+20,i+40)]
  segments(y0=tmp[1,],y1=tmp[3,],x0=c(c(1,2)-.1,3))
  points(x=c(c(1,2)-.1,3),y=tmp[2,],pch=19,cex=.6)
  
  segments(y0=nu_ratio_sum1[1,i],y1=nu_ratio_sum1[3,i],x0=c(1/3)-.1,col='black',lty='dotted')
  points(x=1/3-.1,(nu_ratio_sum1[2,i]),pch=21,bg='white',cex=.6)
  
  segments(y0=nu_ratio_sum2[1,i],y1=nu_ratio_sum2[3,i],x0=c(1+1/3)-.1,col='black',lty='dotted')
  points(x=1+1/3-.1,(nu_ratio_sum2[2,i]),pch=21,bg='white',cex=.6)
  
  segments(y0=nu_ratio_sum3[1,i],y1=nu_ratio_sum3[3,i],x0=c(2+1/3),col='black',lty='dotted')
  points(x=2+1/3,(nu_ratio_sum3[2,i]),pch=21,bg='white',cex=.6)
  
  # round 2
  tmp = nu2.sum[,c(i,i+20)]
  segments(y0=tmp[1,],y1=tmp[3,],x0=c(1,2))
  points(x=c(1,2),y=tmp[2,],pch=19,cex=.6)
  
  segments(y0=nu_ratio_sum4[1,i],y1=nu_ratio_sum4[3,i],x0=c(1/3),col='black',lty='dotted')
  points(x=1/3,(nu_ratio_sum4[2,i]),pch=21,bg='white',cex=.6)
  
  segments(y0=nu_ratio_sum5[1,i],y1=nu_ratio_sum5[3,i],x0=c(1+1/3),col='black',lty='dotted')
  points(x=1+1/3,(nu_ratio_sum5[2,i]),pch=21,bg='white',cex=.6)
  
  # round 3
  tmp = as.matrix(nu3.sum[,c(i)])
  segments(y0=tmp[1,],y1=tmp[3,],x0=c(1)+.1)
  points(x=c(1)+.1,y=tmp[2,],pch=19, cex=.6)
  
  segments(y0=nu_ratio_sum6[1,i],y1=nu_ratio_sum6[3,i],x0=c(1/3)+.1,col='black',lty='dotted')
  points(x=1/3+.1,(nu_ratio_sum6[2,i]),pch=21,bg='white', cex=.6)
  
  
  text(x=.5,y=.05,siteName[i])
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Age (years)", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of viability", side = 2, outer = TRUE, line = 2.2)
mtext("Round 1-3", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)




dev.off()
