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
mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("fruitsSamples.rds",mcmcSampleDirectory)]])
data <- readRDS(mcmcSampleDirectory[[grep("fruitsData.rds",mcmcSampleDirectory)]])

# -------------------------------------------------------------------
# Import and organize data
# -------------------------------------------------------------------
dataDirectory = "/Users/Gregor/Dropbox/dataLibrary/postProcessingData-2021/"

countFruitsPerPlantAllPlots <- readRDS(paste0(dataDirectory,"countFruitsPerPlantAllPlots.RDS"))
siteNames <- unique(countFruitsPerPlantAllPlots$site)

countFruitsPerPlantAllPlots <- countFruitsPerPlantAllPlots %>%
  dplyr::rename(y_tfe = countFruitNumberPerPlant) %>%
  dplyr::select(site,year,y_tfe)

countFruitsPerPlantAllPlots$year <- as.character(countFruitsPerPlantAllPlots$year)

countUndamagedDamagedFruitsPerPlantAllPlots <- readRDS(paste0(dataDirectory,"countUndamagedDamagedFruitsPerPlantAllPlots.RDS"))

countUndamagedDamagedFruitsPerPlantAllPlots <- countUndamagedDamagedFruitsPerPlantAllPlots %>%
  dplyr::rename(y_und = countUndamagedFruitNumberPerPlant) %>%
  dplyr::rename(y_dam = countDamagedFruitNumberPerPlant) %>%
  dplyr::rename(site2 = site) %>%
  dplyr::rename(year2 = year) %>%
  dplyr::select(site2,year2,y_und,y_dam) 

countUndamagedDamagedFruitsPerPlantAllPlots$year2 <- as.character(countUndamagedDamagedFruitsPerPlantAllPlots$year2)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Graphical Posterior Predictive Checks
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# ---
# Graphical checks: Total fruit equivalents per plant ----
# ---

z.sim=MCMCchains(mcmcSamples, params = "z_tfe")

n.rows = dim(z.sim)[1];n.cols = dim(z.sim)[2]
y.sim = matrix( rpois(n.rows*n.cols,z.sim), n.rows, n.cols) 
y.obs=data$y_tfe
n.iter=dim(y.sim)[1]
years=2006:2012

pdf(file=paste0(outputDirectory,"totalFruitEquivalentsPerPlant-ppc-population.pdf"),height=6,width=6)

for(i in 1:length(years)){
  par(mfrow = c(4,5),
      oma = c(4,5,0,0) + 0.1,
      mar = c(1,0,1,1) + 0.1)
  
  n.samples = 50
  iter.ind = sample(1:n.iter,n.samples)
  
  for(j in 1:20){
    
    index=data$site==j&data$year==i
    
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
  
  mtext(paste0("Counts of total fruit equivalents per plant (",years[i],")"), side = 1, outer = TRUE, line = 1.5)
  mtext("Density", side = 2, outer = TRUE, line = 2.2)
}
dev.off()


# ---
# Graphical checks: Undamaged fruits per plant ----
# ---

z.sim=MCMCchains(mcmcSamples, params = "z_und")

n.rows = dim(z.sim)[1];n.cols = dim(z.sim)[2]
y.sim = matrix( rpois(n.rows*n.cols,z.sim), n.rows, n.cols) 
y.obs=data$y_und
n.iter=dim(y.sim)[1]
years=2013:2020

pdf(file=paste0(outputDirectory,"undamagedFruitsPerPlant-ppc-population.pdf"),height=6,width=6)

for(i in 1:length(years)){
  par(mfrow = c(4,5),
      oma = c(4,5,0,0) + 0.1,
      mar = c(1,0,1,1) + 0.1)
  
  n.samples = 50
  iter.ind = sample(1:n.iter,n.samples)
  
  for(j in 1:20){
    
    index=data$site2==j&data$year2==i
    
    p.obs=y.obs[index]
    p.sim=as.matrix(y.sim[,index])
    p.sim=as.matrix(p.sim[iter.ind,])
    
    if (dim(p.sim)[2]<2) {NA} else {
      
      list.dens=apply(p.sim,1,density,na.rm=TRUE)
      all.max.y=max(unlist(lapply(list.dens, "[", "y")))
      all.max.x=max(unlist(lapply(list.dens, "[", "x")))
      
      m.max=max(p.obs,na.rm=TRUE)
      m.min=min(p.obs,na.rm=TRUE)
      if(sum(is.na(p.obs))!=length(p.obs)) {dens.obs = density(p.obs,from=m.min,to=m.max,na.rm=TRUE)} else {NA}
      
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
  
  mtext(paste0("Counts of undamaged fruits per plant (",years[i],")"), side = 1, outer = TRUE, line = 1.5)
  mtext("Density", side = 2, outer = TRUE, line = 2.2)
}
dev.off()

# ---
# Graphical checks: Damaged fruits per plant ----
# ---

z.sim=MCMCchains(mcmcSamples, params = "z_dam")

n.rows = dim(z.sim)[1];n.cols = dim(z.sim)[2]
y.sim = matrix( rpois(n.rows*n.cols,z.sim), n.rows, n.cols) 
y.obs=data$y_dam
n.iter=dim(y.sim)[1]
years=2013:2020

pdf(file=paste0(outputDirectory,"damagedFruitsPerPlant-ppc-population.pdf"),height=6,width=6)

for(i in 1:length(years)){
  par(mfrow = c(4,5),
      oma = c(4,5,0,0) + 0.1,
      mar = c(1,0,1,1) + 0.1)
  
  n.samples = 50
  iter.ind = sample(1:n.iter,n.samples)
  
  for(j in 1:20){
    
    index=data$site2==j&data$year2==i
    
    p.obs=y.obs[index]
    p.sim=as.matrix(y.sim[,index])
    p.sim=as.matrix(p.sim[iter.ind,])
    
    if (dim(p.sim)[2]<2) {NA} else {
      
      list.dens=apply(p.sim,1,density,na.rm=TRUE)
      all.max.y=max(unlist(lapply(list.dens, "[", "y")))
      all.max.x=max(unlist(lapply(list.dens, "[", "x")))
      
      m.max=max(p.obs,na.rm=TRUE)
      m.min=min(p.obs,na.rm=TRUE)
      if(sum(is.na(p.obs))!=length(p.obs)) {dens.obs = density(p.obs,from=m.min,to=m.max,na.rm=TRUE)} else {NA}
      
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
  
  mtext(paste0("Counts of damaged fruits per plant (",years[i],")"), side = 1, outer = TRUE, line = 1.5)
  mtext("Density", side = 2, outer = TRUE, line = 2.2)
}
dev.off()


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Posterior Predictive Checks
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# ---
# Test statistics: Total fruit equivalents per plant ----
# ---
z.sim=MCMCchains(mcmcSamples, params = "z_tfe")
n.rows = dim(z.sim)[1];n.cols = dim(z.sim)[2]
y.sim = matrix( rpois(n.rows*n.cols,z.sim), n.rows, n.cols) 
y.obs=data$y_tfe
n.iter=dim(y.sim)[1]
years=2006:2012


chi2.obs=(sweep(z.sim, 2, y.obs, FUN = '-')^2)/z.sim
chi2.sim=((y.sim-z.sim)^2)/z.sim

n.iter = dim(chi2.obs)[1]

p.chi.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=length(years))
for(j in 1:20){
  
  for(i in 1:length(years)){
    index = data$site==j&data$year==i
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
      index = data$site==j&data$year==i
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
df=data$y_tfe

tfe.min=f(y.sim=sims,y.obs=df,model.fun=min)
tfe.max=f(y.sim=sims,y.obs=df,model.fun=max)
tfe.mean=f(y.sim=sims,y.obs=df,model.fun=mean)
tfe.sd=f(y.sim=sims,y.obs=df,model.fun=sd)

colfunc <- colorRampPalette(c("white", "black"))
col.vec=colfunc(length(years))

pdf(file=paste0(outputDirectory,"totalFruitEquivalentsPerPlant-pvals.pdf"),height=6,width=6)

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
  
  f(x=start[1]+offset[j],y=tfe.min[,j],width=.2)
  
  f(x=start[2]+offset[j],y=tfe.max[,j],width=.2)
  
  f(x=start[3]+offset[j],y=tfe.mean[,j],width=.2)
  
  f(x=start[4]+offset[j],y=tfe.sd[,j],width=.2)
  
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
# Test statistics: Undamaged fruits ----
# ---
z.sim=MCMCchains(mcmcSamples, params = "z_und")
n.rows = dim(z.sim)[1];n.cols = dim(z.sim)[2]
y.sim = matrix( rpois(n.rows*n.cols,z.sim), n.rows, n.cols) 
y.obs=data$y_und
n.iter=dim(y.sim)[1]
years=2013:2020


chi2.obs=(sweep(z.sim, 2, y.obs, FUN = '-')^2)/z.sim
chi2.sim=((y.sim-z.sim)^2)/z.sim

n.iter = dim(chi2.obs)[1]

p.chi.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=length(years))
for(j in 1:20){
  
  for(i in 1:length(years)){
    index = data$site2==j&data$year2==i
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
      index = data$site2==j&data$year2==i
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
df=data$y_und

und.min=f(y.sim=sims,y.obs=df,model.fun=min)
und.max=f(y.sim=sims,y.obs=df,model.fun=max)
und.mean=f(y.sim=sims,y.obs=df,model.fun=mean)
und.sd=f(y.sim=sims,y.obs=df,model.fun=sd)

colfunc <- colorRampPalette(c("white", "black"))
col.vec=colfunc(length(years))

pdf(file=paste0(outputDirectory,"undamagedFruitsPerPlant-pvals.pdf"),height=6,width=6)

par(mfrow=c(1,1))
time.sample = 1:length(years)
plot(NA,NA,
     ylim=c(0,1),pch=16,xlim=c(-.5,24.5),
     xlab="",ylab="",
     main="Undamaged fruits per plant",type='n',
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
  
  f(x=start[1]+offset[j],y=und.min[,j],width=.2)
  
  f(x=start[2]+offset[j],y=und.max[,j],width=.2)
  
  f(x=start[3]+offset[j],y=und.mean[,j],width=.2)
  
  f(x=start[4]+offset[j],y=und.sd[,j],width=.2)
  
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
# Test statistics: Damaged fruits ----
# ---
z.sim=MCMCchains(mcmcSamples, params = "z_dam")
n.rows = dim(z.sim)[1];n.cols = dim(z.sim)[2]
y.sim = matrix( rpois(n.rows*n.cols,z.sim), n.rows, n.cols) 
y.obs=data$y_dam
n.iter=dim(y.sim)[1]
years=2013:2020

chi2.obs=(sweep(z.sim, 2, y.obs, FUN = '-')^2)/z.sim
chi2.sim=((y.sim-z.sim)^2)/z.sim

n.iter = dim(chi2.obs)[1]

p.chi.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=length(years))
for(j in 1:20){
  
  for(i in 1:length(years)){
    index = data$site2==j&data$year2==i
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
      index = data$site2==j&data$year2==i
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
df=data$y_dam

dam.min=f(y.sim=sims,y.obs=df,model.fun=min)
dam.max=f(y.sim=sims,y.obs=df,model.fun=max)
dam.mean=f(y.sim=sims,y.obs=df,model.fun=mean)
dam.sd=f(y.sim=sims,y.obs=df,model.fun=sd)

colfunc <- colorRampPalette(c("white", "black"))
col.vec=colfunc(length(years))

pdf(file=paste0(outputDirectory,"damagedFruitsPerPlant-pvals.pdf"),height=6,width=6)

par(mfrow=c(1,1))
time.sample = 1:length(years)
plot(NA,NA,
     ylim=c(0,1),pch=16,xlim=c(-.5,24.5),
     xlab="",ylab="",
     main="Damaged fruits per plant",type='n',
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
  
  f(x=start[1]+offset[j],y=dam.min[,j],width=.2)
  
  f(x=start[2]+offset[j],y=dam.max[,j],width=.2)
  
  f(x=start[3]+offset[j],y=dam.mean[,j],width=.2)
  
  f(x=start[4]+offset[j],y=dam.sd[,j],width=.2)
  
  f(x=start[5]+offset[j],y=p.chi.mat[,j],width=.2)
  
}

axis(1, c(2,7,12,17,22),
     labels = c("min","max","mean","sd","Chi-2"), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(2, seq(0,1,by=.2),
     seq(0,1,by=.2), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)

dev.off()

