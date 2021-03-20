rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)

library(MCMCvis)
library(tidybayes)
library(tidyverse)
library(magrittr)
library(bayesplot)
library(rethinking)

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/decayModel/"
modelFittingFiles <- paste0(directory,list.files(directory))

mcmcSamples <- readRDS(modelFittingFiles[[grep("seedBurialSamples.rds",modelFittingFiles)]])
data <- readRDS(modelFittingFiles[[1]])

censusSeedlingsFruitingPlants <- readRDS("~/Dropbox/dataLibrary/postProcessingData/censusSeedlingsFruitingPlants.RDS")
siteNames = unique(censusSeedlingsFruitingPlants$site)

# ---
# Convergence diagnostics -------------------------------------------------------------------
# ---

# MCMCsummary(mcmcSamples, params = c("mu0_g"))
# MCMCsummary(mcmcSamples, params = c("sigma0_g"))
# MCMCsummary(mcmcSamples, params = c("sigma_g"))
# 
# MCMCsummary(mcmcSamples, params = c("mu0_s"))
# MCMCsummary(mcmcSamples, params = c("sigma0_s"))
# MCMCsummary(mcmcSamples, params = c("sigma_s"))
# 
# MCMCsummary(mcmcSamples, params = c("a"))
# 
# alpha<-MCMCchains(mcmcSamples, params = c("a"))
# alpha.sum<-apply(alpha,2,quantile,c(.025,.5,.975))
# 
# par(mfrow=c(1,1))
# 
# plot(NA,NA,type='n',xlim=c(0,2),ylim=c(0,20))
# for(i in 1:20){
#   tmp<-alpha.sum[,i]
#   segments(x0=tmp[1],x1=tmp[3],y0=i)
#   points(x=tmp[2],y=i,pch=19)
# }
# 
# for(i in 1:20){
#   hist(MCMCchains(mcmcSamples,params="a")[,i],
#        breaks=100,freq=FALSE,xlim=c(0,2));
#   abline(v=1,col='red')
# }
# 
# diag.obj = gelman.diag(mcmcSamples)
# plot(diag.obj$psrf[,1]);abline(v=c(21,61,81,101))
# names(diag.obj$psrf[,1])
# (diag.obj$psrf[,1])[order(diag.obj$psrf[,1])]

# ---
# Graphical checks ----
# ---

f = function(x="parm",model="belowground",jpeg.quality = 75){
  # get chains for parameter as a list
  
  chains <- MCMCchains(mcmcSamples, params=x, mcmc.list=TRUE)
  chains.list <- lapply(chains, as.matrix)
  
  n = (dim(chains.list[[1]])[2])/20
  counter <- 1
  parm = x
  while(counter <= n) {
    jpeg(filename=paste0("~/Dropbox/clarkiaSeedBanks/products/figures/convergence/convergence-",model,"-",parm,"-",counter,".jpeg"),quality=75)
    
    par.names =  colnames(chains[[1]])
    par(mfrow = c(4,5), oma = c(5,4,0,0) + 0.1, mar = c(0,0,1,1) + 0.1)
    for(i in 1:20){
      plot(as.vector(chains[[1]][,(20*(counter-1)+i)]),type='n',axes=FALSE,frame=FALSE,xaxt='n',yaxt='n',xlab="",ylab="")
      lines(as.vector(chains[[1]][,(20*(counter-1)+i)]),col=rgb(0.9882353, 0.5529412, 0.3490196 ,.5))
      lines(as.vector(chains[[2]][,(20*(counter-1)+i)]),col=rgb(1,1,.75,.5))
      lines(as.vector(chains[[3]][,(20*(counter-1)+i)]),col=rgb(0.5686275, 0.7490196, 0.8588235,.5))
      text(.5*length(as.vector(chains[[1]][,(20*(counter-1)+i)])),.9*max(as.vector(chains[[1]][,(20*(counter-1)+i)])),par.names[(20*(counter-1)+i)])
    }
    dev.off()
    counter <- counter + 1
  }
  
}

f(x="mu0_g",model="belowground")
f(x="sigma0_g",model="belowground")
f(x="sigma_g",model="belowground")

all.chains=MCMCchains(mcmcSamples,params=c("mu0_g","sigma0_g","sigma_g"),mcmc.list=TRUE)

par(mfrow=c(1,3))
hist(coda::gelman.diag(all.chains, confidence = 0.95)$psrf[,1],
     col='black',border='white',breaks=25,main="Distribution of R-hat: germination")
hd=coda::heidel.diag(all.chains)
p = c()
for(i in 1:3){
  p[i] = sum(hd[[i]][,1])/length(hd[[i]][,1])
}
p=signif(p,2)
dt=hist(coda::gelman.diag(all.chains, confidence = 0.95)$psrf[,1],
    breaks=25,plot=FALSE)
text(.999*max(dt$mids),.9*max(dt$counts),p[1])
text(.999*max(dt$mids),.85*max(dt$counts),p[2])
text(.999*max(dt$mids),.8*max(dt$counts),p[3])


f(x="a",model="belowground")
f(x="mu0_s",model="belowground")
f(x="sigma0_s",model="belowground")
f(x="sigma_s",model="belowground")

all.chains=MCMCchains(mcmcSamples,params=c("a","mu0_s","sigma0_s","sigma_s"),mcmc.list=TRUE)

hist(coda::gelman.diag(all.chains, confidence = 0.95)$psrf[,1],
     col='black',border='white',breaks=25,main="Distribution of R-hat: seed survival")
hd=coda::heidel.diag(all.chains)
p = c()
for(i in 1:3){
  p[i] = sum(hd[[i]][,1])/length(hd[[i]][,1])
}
p=signif(p,2)
dt=hist(coda::gelman.diag(all.chains, confidence = 0.95)$psrf[,1],
        breaks=25,plot=FALSE)
text(.999*max(dt$mids),.9*max(dt$counts),p[1])
text(.999*max(dt$mids),.85*max(dt$counts),p[2])
text(.999*max(dt$mids),.8*max(dt$counts),p[3])

f(x="mu0_s0",model="belowground")
f(x="sigma0_s0",model="belowground")
f(x="sigma_s0",model="belowground")

all.chains=MCMCchains(mcmcSamples,params=c("mu0_s0","sigma0_s0","sigma_s0"),mcmc.list=TRUE)

hist(coda::gelman.diag(all.chains, confidence = 0.95)$psrf[,1],
     col='black',border='white',breaks=25,main="Distribution of R-hat: s0")
hd=coda::heidel.diag(all.chains)
p = c()
for(i in 1:3){
  p[i] = sum(hd[[i]][,1])/length(hd[[i]][,1])
}
p=signif(p,2)
dt=hist(coda::gelman.diag(all.chains, confidence = 0.95)$psrf[,1],
        breaks=25,plot=FALSE)
text(.999*max(dt$mids),.9*max(dt$counts),p[1])
text(.999*max(dt$mids),.85*max(dt$counts),p[2])
text(.999*max(dt$mids),.8*max(dt$counts),p[3])


            