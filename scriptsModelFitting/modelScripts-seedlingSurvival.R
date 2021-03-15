# list of JAGS scripts used in this file

# -------------------------------------------------------------------
# Models for ...
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

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# setwd and read data files
censusSeedlingsFruitingPlants <- readRDS("~/Dropbox/dataLibrary/postProcessingData/censusSeedlingsFruitingPlants.RDS")

# -------------------------------------------------------------------
# Clean and organize data
# REMOVE DATA WHERE NUMBER OF FRUITING PLANTS IS GREATER THAN NUMBER OF SEEDLINGS
# -------------------------------------------------------------------
censusSeedlingsFruitingPlants <- censusSeedlingsFruitingPlants %>% 
  dplyr::filter(!(fruitplNumber>seedlingNumber)) #%>%
  # FOR PRIOR PREDICTIVE, SELECT ONE SITE
  # dplyr::filter(site%in%c("BG","LCE"))

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------
# -------------------------------------------------------------------
censusSeedlingsFruitingPlants$year <- as.character(censusSeedlingsFruitingPlants$year)

censusSeedlingsFruitingPlants=censusSeedlingsFruitingPlants %>%
  dplyr::filter(site=="LCW")

data <- tidybayes::compose_data(censusSeedlingsFruitingPlants)

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
n.thin = 10

#set.seed(10)
dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/priorChecking/jagsScripts/")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Describe pooling structure
# -------------------------------------------------------------------
# -------------------------------------------------------------------

initsMu0 <- function(samps = data$n_site){
  rnorm(n = samps, mean = 0, sd = 1)
}

initsSigma0 <- function(samps = data$n_site){
  extraDistr::rhnorm(n = samps, sigma = 1)
}

initsSigma <- function(rows = data$n_site, cols = data$n_year){
  matrix(extraDistr::rhnorm(n = rows*cols, sigma = 1), rows, cols)
}

# set inits for JAGS
inits <- list()
for(i in 1:3){
  inits[[i]] <- list(initsMu0(), initsSigma0(), initsSigma() )
  
  names(inits[[i]]) = c("mu0","sigma0","sigma")
  
}

# # Call to JAGS
# 
# # tuning (n.adapt)
jm = jags.model(paste0(dir,"binomialLikelihood-logitLink-normalHierarchical-3.R"), 
                data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

#parsToMonitor = c("mu0","sigma0","mu","sigma","theta")
parsToCheck = c(# LOG LIKELIHOODS
                # "logLik",
                # POSTERIOR PREDICTIVE
                "fruitplNumber_sim",
                "e.mu","mu",
                # MODEL CHECKING 
               "chi2.obs","chi2.sim")


# chain (n.iter)
samples.rjags = coda.samples(jm, 
                             variable.names = c(parsToCheck), 
                             n.iter = n.iterations, thin = n.thin)

fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedlingSurvival/")
dir.create(file.path(fileDirectory), showWarnings = FALSE)
#
saveRDS(samples.rjags,file=paste0(fileDirectory,"seedSurvivalSamplesChecksPooling.rds"))
#saveRDS(samples.rjags,file=paste0("/Users/Gregor/Dropbox/dataLibrary/posteriors/seedSurvivalSamplesCheck.rds"))
#saveRDS(data,file=paste0(fileDirectory,"data.rds"))
#saveRDS(censusSeedlingsFruitingPlants,file=paste0(fileDirectory,"censusSeedlingsFruitingPlants.rds"))
# 
# # -------------------------------------------------------------------
# # -------------------------------------------------------------------
# # Evaluate priors
# # -------------------------------------------------------------------
# # -------------------------------------------------------------------
# 
# hist.chain=function(x){
#   name<-colnames(x)
#   tmp = hist(x,breaks=100,plot=FALSE)
#   tmp$xname = name
#   plot(tmp,freq=FALSE)
# }
# 
# ## Priors on hyperparameters
# par(mfrow=c(1,3))
# hist.chain(MCMCchains(samples.rjags,params="mu0"))
# hist.chain(MCMCchains(samples.rjags,params="sigma0"))
# hist.chain(MCMCchains(samples.rjags,params="sigma")[,1,drop=FALSE])
# hist.chain(MCMCchains(samples.rjags,params="mu")[,1,drop=FALSE])
# hist.chain(boot::inv.logit(MCMCchains(samples.rjags,params="mu")[,1,drop=FALSE]))
# 
# for(i in 1){
#   out=apply(MCMCchains(samples.rjags,params="mu"),2,boot::inv.logit)
#   index=grep(paste0("\\[",2,","),colnames(out))
#   tmp=out[,index]
#   
#   par(mfrow=c(3,5))
#  for(i in 1:14){ hist.chain(tmp[,i]); 
#    abline(v = median(tmp[,i]),col='red',lwd=2); 
#    abline(v=rethinking::chainmode(tmp[,i]),col='blue',lwd=2);
#    abline(v=coda::HPDinterval(as.mcmc(tmp[,i]), prob = 0),col='orange',lwd=2,lty='dotted')}
#   
#   lines(x=2006:2019,y=tmp[3,],type='b')
#   segments(y0=tmp[1,],y1=tmp[5,],x0=2006:2019)
#   segments(y0=tmp[2,],y1=tmp[4,],x0=2006:2019,lwd=3)
#   points(y=tmp[3,],x=2006:2019,pch=21,bg='white')
#   
#   
#   text(x=2007,y=.9,siteNames[i],pos=4)
#   ifelse(i%in%c(16:20),axis(1L, at = c(2006,2008,2010,2012,2014,2016,2018,2020)),NA)
#   ifelse(i%in%c(1,6,11,16),axis(2L, at = c(0, .2, .4, .6, .8, 1)),NA)
# }
# 
# 
# ## Priors on transformed parameters
# n_params=dim(MCMCchains(samples.rjags,params=c("theta")))[2]
# par(mfrow=c(1,3))
# hist.chain(MCMCchains(samples.rjags,params=c("theta"))[,sample(1:n_params,1),drop=FALSE])
# hist.chain(MCMCchains(samples.rjags,params=c("theta"))[,sample(1:n_params,1),drop=FALSE])
# hist.chain(MCMCchains(samples.rjags,params=c("theta"))[,sample(1:n_params,1),drop=FALSE])
# 
# ## Prior predictive
# y_prior=MCMCchains(samples.rjags,params="y_prior")
# 
# par(mfrow=c(1,2))
# plot(data$fruitplNumber,y_prior[sample(dim(y_prior)[1],1),],
#      type='n',ylim=c(0,200),
#      xlab="Observed data",ylab="Simulated data")
# for(i in 1:100){
#   points(data$fruitplNumber,y_prior[sample(dim(y_prior)[1],1),],
#          pch=16,cex=.5)
# }
# abline(a=0,b=1)
# 
# # easier to see the full range from the prior in a proportion plot
# # certain values are only observed in certain combinations because of denominator (number of trials)
# plot(data$fruitplNumber/data$seedlingNumber,
#      y_prior[sample(dim(y_prior)[1],1),]/data$seedlingNumber,
#      type='n',xlim=c(0,1),ylim=c(0,1),
#      xlab="Observed data",ylab="Simulated data")
# for(i in 1:100){
#   points(data$fruitplNumber/data$seedlingNumber,
#          y_prior[sample(dim(y_prior)[1],1),]/data$seedlingNumber,
#          pch=16,cex=.25)
# }
# abline(a=0,b=1)
# 
# # only 3 plots with >10 seedlings have survival>.8
# censusSeedlingsFruitingPlants %>%
#   dplyr::mutate(p = fruitplNumber/seedlingNumber) %>%
#   dplyr::filter(p>0.8)
# # 
# # ## SUMMARY
# # MCMCsummary(samples.rjags,params=c("mu0","sigma0"))
# # MCMCsummary(samples.rjags,params=c("sigma"))
# # 
# MCMCsummary(samples.rjags,params=c("mu"))
# # 
# # mu<-MCMCchains(samples.rjags,params=c("mu"))
# # mu.med=apply(mu,2,median)
# # mu.mat = matrix(NA,nrow=14,ncol=20)
# # for(i in 1:20){
# #   index=grep(paste0("\\[",i,","),names(mu.med))
# #   mu.mat[,i]=mu.med[index]
# # }
# # mu.mat=boot::inv.logit(mu.mat)
# # 
# # 
# directory2 = "/Users/Gregor/Dropbox/dataLibrary/workflow/tidyData/"
# dataFiles <- paste0(directory2,list.files(directory2))
# 
# data2 <- readRDS(dataFiles[[1]])
# siteNames = unique(data2$siteBags)
# 
# position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>%
#   dplyr::select(site,easting,dominant.surface.rock.type) %>%
#   dplyr::mutate(easting=easting/1000)
# # 
# # siteSurvival=MCMCchains(samples.rjags,params=c("mu0"))
# # siteSurvival=apply(siteSurvival,2,boot::inv.logit)
# # 
# # f.param = function(x.v,parm="sigma"){
# #   
# #   x.sum=apply(x.v,2,quantile,probs=c(0.025,.25,.5,.75,.975))
# #   
# #   plot(NA,NA,type='n',xlim=c(0,1),ylim=c(0,20),
# #        axes=FALSE,frame=FALSE,
# #        xlab="",ylab="")
# #   #abline(v=0,col='gray')
# #   y.pt = 20:1
# #   for(i in 20:1){
# #     tmp<-x.sum[,i]
# #     segments(x0=tmp[1],x1=tmp[5],y0=y.pt[i])
# #     segments(x0=tmp[2],x1=tmp[4],y0=y.pt[i],lwd=3)
# #     points(x=tmp[3],y=y.pt[i],pch=21,bg='white')
# #   }
# #   
# #   axis(1,  seq(0,1,by=.2), col.ticks = 1)
# #   axis(2, (1:20),
# #        labels = rev(siteNames), las = 1, 
# #        col = NA, col.ticks = 1, cex.axis = 1)
# #   mtext("Probability",
# #         side=1,line=2.5,adj=.5,col='black',cex=1)
# #   mtext(paste("Parameter:",parm),
# #         side=3,line=0,adj=0,col='black',cex=.65)
# #   
# #   x.sum.df<-data.frame(t(x.sum),position)
# #   names(x.sum.df)[1:5] = c("ci.lolo","ci.lo","ci.med","ci.hi","ci.hihi")
# #   
# #   plot(NA,NA,type='n',ylim=seq(0,1),xlim=c(340,375),
# #        axes=FALSE,frame=FALSE,
# #        xlab="",ylab="")
# #   #abline(h=0,col='gray')
# #   
# #   segments(x0=x.sum.df$easting,y0=x.sum.df$ci.lolo,y1=x.sum.df$ci.hihi,lwd=1)
# #   segments(x0=x.sum.df$easting,y0=x.sum.df$ci.lo,y1=x.sum.df$ci.hi,lwd=3)
# #   points(x=x.sum.df$easting,y=x.sum.df$ci.med,pch=21,bg='white')
# #   
# #   axis(2, seq(0,1,by=.2), col.ticks = 1)
# #   axis(1, seq(340,375,by=5),
# #        labels = seq(340,375,by=5), las = 1, 
# #        col.ticks = 1, cex.axis = 1)
# #   mtext("Probability",
# #         side=2,line=2.5,adj=.5,col='black',cex=1)
# #   mtext("Easting (km)",
# #         side=1,line=2.5,adj=.5,col='black',cex=1)
# # }
# # 
# # dev.off()
# # par(mfrow=c(1,2))
# # f.param(siteSurvival,parm="sigma")
# # 
# # 
# siteYearSurvival = MCMCchains(samples.rjags,params=c("mu"))
# siteYearSurvival = apply(siteYearSurvival,2,boot::inv.logit)
# 
# x.sum=apply(siteYearSurvival,2,quantile,probs=c(0.025,.25,.5,.75,.975))
# 
# 
# 
# plot(NA,NA,type='n',xlim=c(0,1),ylim=c(0,2),
#      axes=FALSE,frame=FALSE,
#      xlab="",ylab="")
# #abline(v=0,col='gray')
# y.pt = 2:1
# for(i in 2:1){
#   tmp<-x.sum[,i]
#   segments(x0=tmp[1],x1=tmp[5],y0=y.pt[i])
#   segments(x0=tmp[2],x1=tmp[4],y0=y.pt[i],lwd=3)
#   points(x=tmp[3],y=y.pt[i],pch=21,bg='white')
# }
# 
# df.summary = censusSeedlingsFruitingPlants %>%
#   dplyr::group_by(site,year) %>%
#   dplyr::summarise(y=sum(fruitplNumber,na.rm=TRUE),
#                    n=sum(seedlingNumber,na.rm=TRUE)) %>%
#   dplyr::mutate(p = y/n)
# 
# df.summary=df.summary %>%
#   dplyr::left_join(data.frame(site=siteNames,index=1:20)  ) %>%
#   dplyr::left_join(data.frame(year=as.character(2006:2019),index2=1:14))
# 
# 
# 
# par(mfrow = c(1,2),
#     oma = c(5,4,0,0) + 0.1,
#     mar = c(0,0,1,1) + 0.1)
# for(i in 1:2){
#   plot(NA,NA,xlim=c(2006,2020),ylim=c(0,1),ylab='',xlab='',xaxt='n',yaxt='n')
# 
#   i=ifelse(i==1,i,11)
#   tmp=df.summary[df.summary$index==i,]
#   lines(x=2006:2019,y=tmp$p,type='b',col='orange',pch=16)
# 
#   i=ifelse(i==11,2,i)
#   index=grep(paste0("\\[",i,","),colnames(x.sum))
#   tmp=x.sum[,index]
# 
#   lines(x=2006:2019,y=tmp[3,],type='b')
#   segments(y0=tmp[1,],y1=tmp[5,],x0=2006:2019)
#   segments(y0=tmp[2,],y1=tmp[4,],x0=2006:2019,lwd=3)
#   points(y=tmp[3,],x=2006:2019,pch=21,bg='white')
# 
# 
#   text(x=2007,y=.9,siteNames[i],pos=4)
#   ifelse(i%in%c(16:20),axis(1L, at = c(2006,2008,2010,2012,2014,2016,2018,2020)),NA)
#   ifelse(i%in%c(1,6,11,16),axis(2L, at = c(0, .2, .4, .6, .8, 1)),NA)
# }
# #   
# # 
# # par(mfrow = c(4,5),
# #     oma = c(5,4,0,0) + 0.1,
# #     mar = c(0,0,1,1) + 0.1)
# # for(i in 1:20){
# # 
# #   tmp.freq=df.summary[df.summary$index==i,]
# #  # lines(x=2006:2019,y=tmp.freq$p,type='b',col='orange',pch=16)
# #   
# #   index=grep(paste0("\\[",i,","),colnames(x.sum))
# #   tmp=x.sum[,index]
# #   
# #   plot(NA,NA,xlim=c(0,1),ylim=c(0,1),ylab='',xlab='',xaxt='n',yaxt='n')
# #   abline(a=0,b=1,col='gray')
# #   points(tmp.freq$p,tmp[3,])
# # 
# #   
# #   text(x=.1,y=.9,siteNames[i],pos=4)
# #   ifelse(i%in%c(16:20),axis(1L, at = c(0, .2, .4, .6, .8, 1)),NA)
# #   ifelse(i%in%c(1,6,11,16),axis(2L, at = c(0, .2, .4, .6, .8, 1)),NA)
# # }
# # 
# # 
# # 
# # 
# # counts = censusSeedlingsFruitingPlants %>%
# #   dplyr::group_by(site,year) %>%
# #   dplyr::summarise(number = sum(seedlingNumber,na.rm=TRUE)) %>%
# #   dplyr::left_join(data.frame(site=siteNames,index=1:20)  ) %>%
# #   dplyr::left_join(data.frame(year=as.character(2006:2019),index2=1:14))
# # 
# # par(mfrow = c(4,5),
# #     oma = c(5,4,0,0) + 0.1,
# #     mar = c(0,0,1,1) + 0.1)
# # for(i in 1:20){
# #   
# #   
# #   index=grep(paste0("\\[",i,","),colnames(x.sum))
# #   tmp=x.sum[,index]
# #   
# #   tmp.freq=df.summary[df.summary$index==i,]
# #   del=tmp[3,]-tmp.freq$p
# #   
# #   tmp.freq=counts[counts$index==i,]
# #   
# #   plot(NA,NA,xlim=c(0,max(tmp.freq$number)),ylim=c(-.2,.2),ylab='',xlab='')
# #   abline(h=0,col='gray')
# #  # points(tmp.freq$number,tmp[3,])
# #  # segments(y0=tmp[1,],y1=tmp[5,],x0=tmp.freq$number)
# #  # segments(y0=tmp[2,],y1=tmp[4,],x0=tmp.freq$number,lwd=3)
# #   points(y=del,x=tmp.freq$number,pch=21,bg='white')
# #   
# #   text(x=.1,y=.9,siteNames[i],pos=4)
# #  # ifelse(i%in%c(16:20),axis(1L, at = c(0, .2, .4, .6, .8, 1)),NA)
# #  # ifelse(i%in%c(1,6,11,16),axis(2L, at = c(0, .2, .4, .6, .8, 1)),NA)
# # }
# # 
# # 
# # dev.off()
# # del=x.sum[3,]-df.summary$p
# # 
# # plot(NA,NA,xlim=c(0,max(counts$number)),ylim=c(-1,1),ylab='',xlab='')
# # abline(h=0,col='gray')
# # points(y=del,x=counts$number,pch=16,cex=.5)
# # 
