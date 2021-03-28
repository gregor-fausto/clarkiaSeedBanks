# rm(list=ls(all=TRUE)) # clear R environment
rm(list=setdiff(ls(all=TRUE),c("scriptCheckDirectory","fileDirectory","outputDirectory"))) # if using in source(script)
options(stringsAsFactors = FALSE)

library(MCMCvis)
library(tidybayes)
library(tidyverse)
library(magrittr)
library(bayesplot)
library(rethinking)

mcmcSampleDirectory <- paste0(fileDirectory,list.files(fileDirectory))
mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("seedsPerFruitSamples.rds",mcmcSampleDirectory)]])
data <- readRDS(mcmcSampleDirectory[[grep("seedsPerFruitData.rds",mcmcSampleDirectory)]])

# -------------------------------------------------------------------
# Import and organize data
# -------------------------------------------------------------------
dataDirectory = "/Users/Gregor/Dropbox/dataLibrary/postProcessingData-2021/"

countSeedPerFruit <-readRDS(paste0(dataDirectory,"countSeedPerFruit.RDS"))
siteNames = unique(countSeedPerFruit[countSeedPerFruit$demography==1,]$site)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Graphical Posterior Predictive Checks
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# ---
# Graphical checks: Seeds per undamaged fruit ----
# ---

y.sim=MCMCchains(mcmcSamples, params = "y_sd.sim")
y.obs=data$sdno
n.iter=dim(y.sim)[1]
years=2006:2020

pdf(file=paste0(outputDirectory,"seedsPerUndamagedFruit-ppc-population.pdf"),height=6,width=6)

for(i in 1:length(years)){
  par(mfrow = c(4,5),
      oma = c(4,5,0,0) + 0.1,
      mar = c(1,0,1,1) + 0.1)
  
  n.samples = 25
  iter.ind = sample(1:n.iter,n.samples)
  
  for(j in 1:20){
    
    index=data$site3==j&data$year3==i
    
    p.obs=y.obs[index]
    p.sim=as.matrix(y.sim[,index])
    p.sim=as.matrix(p.sim[iter.ind,])
    
    if (dim(p.sim)[2]<2) {NA} else {
      
      list.dens=apply(p.sim,1,density,na.rm=TRUE)
      all.max.y=max(unlist(lapply(list.dens, "[", "y")))
      all.max.x=max(unlist(lapply(list.dens, "[", "x")))
      
      m.max=max(p.obs,na.rm=TRUE)
      m.min=min(p.obs,na.rm=TRUE)
      dens.obs = density(p.obs,from=m.min,to=m.max,na.rm=TRUE)
      
      all.max.y = max(all.max.y,max(dens.obs$y))
      all.max.x = max(all.max.x,max(dens.obs$x))
      
      plot(NA,NA,
           ylim=c(0,all.max.y),xlim=c(0,all.max.x),
           xaxt='n',xlab='',ylab='',yaxt='n')
      for(h in 1:n.samples){
        m.max=max(p.sim[h,],na.rm=TRUE)
        m.min=min(p.sim[h,],na.rm=TRUE)
        
        dens.x = density(p.sim[h,],from=m.min,to=m.max,na.rm=TRUE)
        lines(x=dens.x$x,y=dens.x$y,lwd=0.25,col='orange')
      }
      
      
      lines(x=dens.obs$x,y=dens.obs$y)
      
      ifelse(i%in%c(1),axis(2L),NA)
      
      axis(2,cex=.25,tick=FALSE,line=-1)
      axis(1,cex=.25,tick=FALSE,line=-1)
      
      text(x=all.max.x*.7,y=all.max.y*.9,siteNames[j],pos=4)
    }
  }
  
  mtext(paste0("Counts of seeds per undamaged frut (",years[i],")"), side = 1, outer = TRUE, line = 1.5)
  mtext("Density", side = 2, outer = TRUE, line = 2.2)
}
dev.off()

# ---
# Graphical checks: Seeds per damaged fruit ----
# ---

y.sim=MCMCchains(mcmcSamples, params = "y_sd_dam.sim")
y.obs=data$sdno_dam
n.iter=dim(y.sim)[1]
years=2013:2020

pdf(file=paste0(outputDirectory,"seedsPerDamagedFruit-ppc-population.pdf"),height=6,width=6)

for(i in 1:length(years)){
  par(mfrow = c(4,5),
      oma = c(4,5,0,0) + 0.1,
      mar = c(1,0,1,1) + 0.1)
  
  n.samples = 25
  iter.ind = sample(1:n.iter,n.samples)
  
  for(j in 1:20){
    
    index=data$site4==j&data$year4==i
    
    p.obs=y.obs[index]
    p.sim=as.matrix(y.sim[,index])
    p.sim=as.matrix(p.sim[iter.ind,])
    
    if (dim(p.sim)[2]<2) {NA} else {
      
      list.dens=apply(p.sim,1,density,na.rm=TRUE)
      all.max.y=max(unlist(lapply(list.dens, "[", "y")))
      all.max.x=max(unlist(lapply(list.dens, "[", "x")))
      
      m.max=max(p.obs,na.rm=TRUE)
      m.min=min(p.obs,na.rm=TRUE)
      dens.obs = density(p.obs,from=m.min,to=m.max,na.rm=TRUE)
      
      all.max.y = max(all.max.y,max(dens.obs$y))
      all.max.x = max(all.max.x,max(dens.obs$x))
      
      plot(NA,NA,
           ylim=c(0,all.max.y),xlim=c(0,all.max.x),
           xaxt='n',xlab='',ylab='',yaxt='n')
      for(h in 1:n.samples){
        m.max=max(p.sim[h,],na.rm=TRUE)
        m.min=min(p.sim[h,],na.rm=TRUE)
        
        dens.x = density(p.sim[h,],from=m.min,to=m.max,na.rm=TRUE)
        lines(x=dens.x$x,y=dens.x$y,lwd=0.25,col='orange')
      }
      
      lines(x=dens.obs$x,y=dens.obs$y)
      
      ifelse(i%in%c(1),axis(2L),NA)
      
      axis(2,cex=.25,tick=FALSE,line=-1)
      axis(1,cex=.25,tick=FALSE,line=-1)
      
      text(x=all.max.x*.7,y=all.max.y*.9,siteNames[j],pos=4)
    }
  }
  
  mtext(paste0("Counts of seeds per damaged frut (",years[i],")"), side = 1, outer = TRUE, line = 1.5)
  mtext("Density", side = 2, outer = TRUE, line = 2.2)
}
dev.off()


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Posterior Predictive Checks
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# ---
# Test statistics: Seeds per undamaged fruit ----
# ---
y.sim=MCMCchains(mcmcSamples, params = "y_sd.sim")
y.obs=data$sdno
n.iter=dim(y.sim)[1]
years=2006:2020

chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.sd.obs"))
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.sd.sim"))

p.chi.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=length(years))
for(j in 1:20){
  
  for(i in 1:length(years)){
    index = data$site3==j&data$year3==i
    tmp.chi2.obs=as.matrix(chi2.obs[,index])
    tmp.chi2.sim=as.matrix(chi2.sim[,index])
    fit.obs=apply(tmp.chi2.obs,1,sum)
    fit.sim=apply(tmp.chi2.sim,1,sum)
    p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
    p.pop[,i] = p.chi2.calc
  }
  p.chi.list[[j]] = apply(p.pop,2,mean,na.rm=TRUE)
}
p.chi.mat=do.call(rbind,p.chi.list)


f=function(y.sim=chains,y.obs=data,model.fun=mean){
  
  n.iter=dim(y.sim)[1]
  
  p.test.list = list()
  p.pop = matrix(NA,nrow=n.iter,ncol=length(years))
  for(j in 1:20){
    
    for(i in 1:length(years)){
      index = data$site3==j&data$year3==i
      tmp.obs=as.matrix(y.obs[index])
      tmp.sim=as.matrix(y.sim[,index])
      test.obs=model.fun(tmp.obs,na.rm=TRUE)
      test.sim=apply(tmp.sim,1,model.fun,na.rm=TRUE)
      p.test.calc=ifelse(test.sim-test.obs>=0,1,0)
      p.pop[,i] = p.test.calc
    }
    p.test.list[[j]] = apply(p.pop,2,mean,na.rm=TRUE)
  }
  p.test.mat=do.call(rbind,p.test.list)
  return(p.test.mat)
}

sims=y.sim
df=data$sdno

sdno.min=f(y.sim=sims,y.obs=df,model.fun=min)
sdno.max=f(y.sim=sims,y.obs=df,model.fun=max)
sdno.mean=f(y.sim=sims,y.obs=df,model.fun=mean)
sdno.sd=f(y.sim=sims,y.obs=df,model.fun=sd)

colfunc <- colorRampPalette(c("white", "black"))
col.vec=colfunc(length(years))

pdf(file=paste0(outputDirectory,"seedsPerUndamagedFruit-pvals.pdf"),height=6,width=6)

par(mfrow=c(1,1))
time.sample = 1:length(years)
plot(NA,NA,
     ylim=c(0,1),pch=16,xlim=c(-.5,24.5),
     xlab="",ylab="",
     main="Total fruit equivalents per plant",type='n',
     xaxt='n',yaxt='n',frame=FALSE)

start = c(0,5,10,15,20)
offset = seq(0,4,length.out=length(years))#c(-.125,.125)
polygon(x=c(4.5,4.5,9.5,9.5),y=c(-1,2,2,-1),col='gray90',border=0)
polygon(x=c(14.5,14.5,19.5,19.5),y=c(-1,2,2,-1),col='gray90',border=0)
box(which="plot",bty="l",col='black')
abline(h=0.5,lty='dotted')

f=function(x,y,width=.2){
  d=boxplot(y,plot=FALSE)
  segments(x0=x,y0=d$stats[1],y1=d$stats[5],lty='solid')
  rect(xleft=x-width/2,xright=x+width/2,
       ybottom=d$stats[2],ytop=d$stats[4],col='white')
  segments(x0=x-width/2,x1=x+width/2,
           y0=d$stats[3],lwd=2)
  
}

for(j in 1:length(years)){
  
  f(x=start[1]+offset[j],y=sdno.min[,j],width=.2)
  
  f(x=start[2]+offset[j],y=sdno.max[,j],width=.2)
  
  f(x=start[3]+offset[j],y=sdno.mean[,j],width=.2)
  
  f(x=start[4]+offset[j],y=sdno.sd[,j],width=.2)
  
  f(x=start[5]+offset[j],y=p.chi.mat[,j],width=.2)
  
}

axis(1, c(2,7,12,17,22),
     labels = c("min","max","mean","sd","Chi-2"), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(2, seq(0,1,by=.2),
     seq(0,1,by=.2), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)

dev.off()


# ---
# Test statistics: Seeds per damaged fruit ----
# ---
y.sim=MCMCchains(mcmcSamples, params = "y_sd_dam.sim")
y.obs=data$sdno_dam
n.iter=dim(y.sim)[1]
years=2013:2020

chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.sd_dam.obs"))
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.sd_dam.sim"))

p.chi.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=length(years))
for(j in 1:20){
  
  for(i in 1:length(years)){
    index = data$site4==j&data$year4==i
    tmp.chi2.obs=as.matrix(chi2.obs[,index])
    tmp.chi2.sim=as.matrix(chi2.sim[,index])
    fit.obs=apply(tmp.chi2.obs,1,sum)
    fit.sim=apply(tmp.chi2.sim,1,sum)
    p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
    p.pop[,i] = p.chi2.calc
  }
  p.chi.list[[j]] = apply(p.pop,2,mean,na.rm=TRUE)
}
p.chi.mat=do.call(rbind,p.chi.list)


f=function(y.sim=chains,y.obs=data,model.fun=mean){
  
  n.iter=dim(y.sim)[1]
  
  p.test.list = list()
  p.pop = matrix(NA,nrow=n.iter,ncol=length(years))
  for(j in 1:20){
    
    for(i in 1:length(years)){
      index = data$site4==j&data$year4==i
      tmp.obs=as.matrix(y.obs[index])
      tmp.sim=as.matrix(y.sim[,index])
      test.obs=model.fun(tmp.obs,na.rm=TRUE)
      test.sim=apply(tmp.sim,1,model.fun,na.rm=TRUE)
      p.test.calc=ifelse(test.sim-test.obs>=0,1,0)
      p.pop[,i] = p.test.calc
    }
    p.test.list[[j]] = apply(p.pop,2,mean,na.rm=TRUE)
  }
  p.test.mat=do.call(rbind,p.test.list)
  return(p.test.mat)
}

sims=y.sim
df=data$sdno

sdno.min=f(y.sim=sims,y.obs=df,model.fun=min)
sdno.max=f(y.sim=sims,y.obs=df,model.fun=max)
sdno.mean=f(y.sim=sims,y.obs=df,model.fun=mean)
sdno.sd=f(y.sim=sims,y.obs=df,model.fun=sd)

colfunc <- colorRampPalette(c("white", "black"))
col.vec=colfunc(length(years))

pdf(file=paste0(outputDirectory,"seedsPerDamagedFruit-pvals.pdf"),height=6,width=6)

par(mfrow=c(1,1))
time.sample = 1:length(years)
plot(NA,NA,
     ylim=c(0,1),pch=16,xlim=c(-.5,24.5),
     xlab="",ylab="",
     main="Total fruit equivalents per plant",type='n',
     xaxt='n',yaxt='n',frame=FALSE)

start = c(0,5,10,15,20)
offset = seq(0,4,length.out=length(years))#c(-.125,.125)
polygon(x=c(4.5,4.5,9.5,9.5),y=c(-1,2,2,-1),col='gray90',border=0)
polygon(x=c(14.5,14.5,19.5,19.5),y=c(-1,2,2,-1),col='gray90',border=0)
box(which="plot",bty="l",col='black')
abline(h=0.5,lty='dotted')

f=function(x,y,width=.2){
  d=boxplot(y,plot=FALSE)
  segments(x0=x,y0=d$stats[1],y1=d$stats[5],lty='solid')
  rect(xleft=x-width/2,xright=x+width/2,
       ybottom=d$stats[2],ytop=d$stats[4],col='white')
  segments(x0=x-width/2,x1=x+width/2,
           y0=d$stats[3],lwd=2)
  
}

for(j in 1:length(years)){
  
  f(x=start[1]+offset[j],y=sdno.min[,j],width=.2)
  
  f(x=start[2]+offset[j],y=sdno.max[,j],width=.2)
  
  f(x=start[3]+offset[j],y=sdno.mean[,j],width=.2)
  
  f(x=start[4]+offset[j],y=sdno.sd[,j],width=.2)
  
  f(x=start[5]+offset[j],y=p.chi.mat[,j],width=.2)
  
}

axis(1, c(2,7,12,17,22),
     labels = c("min","max","mean","sd","Chi-2"), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(2, seq(0,1,by=.2),
     seq(0,1,by=.2), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)

dev.off()

