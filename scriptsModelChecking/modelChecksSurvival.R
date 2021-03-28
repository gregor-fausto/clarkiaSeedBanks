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
mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("seedlingSurvivalSamples.rds",mcmcSampleDirectory)]])
data <- readRDS(mcmcSampleDirectory[[grep("seedlingSurvivalData.rds",mcmcSampleDirectory)]])

censusSeedlingsFruitingPlants <- readRDS("~/Dropbox/dataLibrary/postProcessingData-2021/censusSeedlingsFruitingPlants.RDS")
siteNames = unique(censusSeedlingsFruitingPlants$site)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Graphical Posterior Predictive Checks
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# ---
# Graphical checks ----
# ---

y.sim=MCMCchains(mcmcSamples, params = "fruitplNumber_sim")
y.obs=data$fruitplNumber
n.obs=data$seedlingNumber
n.iter=dim(y.sim)[1]
years=2006:2020

pdf(file=paste0(outputDirectory,"seedlingSurvivalFruiting-ppc-population.pdf"),height=6,width=6)

for(i in 1:length(years)){
  par(mfrow = c(4,5),
      oma = c(4,5,0,0) + 0.1,
      mar = c(1,0,1,1) + 0.1)
  
  n.samples = 50
  iter.ind = sample(1:n.iter,n.samples)
  
  for(j in 1:20){
   
    index=data$site==j&data$year==i
   
    p.obs=y.obs[index]
    p.sim=y.sim[,index]
    p.sim=p.sim[iter.ind,]
    
    list.dens=apply(p.sim,1,density,na.rm=TRUE)
    all.max.y=max(unlist(lapply(list.dens, "[", "y")))
    all.max.x=max(unlist(lapply(list.dens, "[", "x")))
    
    plot(NA,NA,
         ylim=c(0,all.max.y),xlim=c(0,all.max.x),
         xaxt='n',xlab='',ylab='',yaxt='n')
    for(h in 1:n.samples){
      m.max=max(p.sim[h,],na.rm=TRUE)
      m.min=min(p.sim[h,],na.rm=TRUE)
      
      dens.x = density(p.sim[h,],from=m.min,to=m.max,na.rm=TRUE)
      lines(x=dens.x$x,y=dens.x$y,lwd=0.25,col='orange')
    }
    
    p = data$fruitplNumber[index]
    m.max=max(p,na.rm=TRUE)
    m.min=min(p,na.rm=TRUE)
    
    dens.x = density(p,from=m.min,to=m.max,na.rm=TRUE)
    
    lines(x=dens.x$x,y=dens.x$y)
    
    ifelse(i%in%c(1),axis(2L),NA)
    
    axis(2,cex=.25,tick=FALSE,line=-1)
    axis(1,cex=.25,tick=FALSE,line=-1)
    
    text(x=all.max.x*.7,y=all.max.y*.9,siteNames[j],pos=4)
  }
  
  mtext(paste0("Probability of seedling survival (",years[i],")"), side = 1, outer = TRUE, line = 1.5)
  mtext("Density", side = 2, outer = TRUE, line = 2.2)
}
dev.off()


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Posterior Predictive Checks
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# ---
# Test statistics ----
# ---

chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.obs"))
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.sim"))

n.iter = dim(chi2.obs)[1]

p.chi.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=length(years))
for(j in 1:20){
  
  for(i in 1:length(years)){
    index = data$site==j&data$year==i
    tmp.chi2.obs=chi2.obs[,index]
    tmp.chi2.sim=chi2.sim[,index]
    fit.obs=apply(tmp.chi2.obs,1,sum)
    fit.sim=apply(tmp.chi2.sim,1,sum)
    p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
    p.pop[,i] = p.chi2.calc
  }
  p.chi.list[[j]] = apply(p.pop,2,mean,na.rm=TRUE)
}
p.chi.mat=do.call(rbind,p.chi.list)


f=function(y.sim=chains,y.obs=data,n.obs=data2,model.fun=mean){
  
  n.iter=dim(y.sim)[1]
  
  p.test.list = list()
  p.pop = matrix(NA,nrow=n.iter,ncol=length(years))
  for(j in 1:20){
    
    for(i in 1:length(years)){
      index = data$site==j&data$year==i
      tmp.obs=y.obs[index]/n.obs[index]
      tmp.sim=sweep(y.sim[,index], 2, n.obs[index], FUN = '/')
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

sims=MCMCchains(mcmcSamples,params=c("fruitplNumber_sim"))
df=data$fruitplNumber
df2=data$seedlingNumber

seedlings.min=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=min)
seedlings.max=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=max)
seedlings.mean=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=mean)
seedlings.sd=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=sd)

colfunc <- colorRampPalette(c("white", "black"))
col.vec=colfunc(15)

pdf(file=paste0(outputDirectory,"seedlingSurvivalFruiting-pvals.pdf"),height=6,width=6)

par(mfrow=c(1,1))
time.sample = 1:length(years)
plot(NA,NA,
     ylim=c(0,1),pch=16,xlim=c(-.5,24.5),
     xlab="",ylab="",
     main="Seedling survival to fruiting",type='n',
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

  f(x=start[1]+offset[j],y=seedlings.min[,j],width=.2)
  
  f(x=start[2]+offset[j],y=seedlings.max[,j],width=.2)
  
  f(x=start[3]+offset[j],y=seedlings.mean[,j],width=.2)
  
  f(x=start[4]+offset[j],y=seedlings.sd[,j],width=.2)
  
  f(x=start[5]+offset[j],y=p.chi.mat[,j],width=.2)
  
}

segments(x0=start[2]+offset[8]-.025,
         x1=start[2]+offset[11]+.025,
         y0=1.02,lty='solid',lwd=2)
segments(x0=start[3]+offset[8]-.025,
         x1=start[3]+offset[11]+.025,
         y0=1.02,lty='solid',lwd=2)
segments(x0=start[4]+offset[8]-.025,
         x1=start[4]+offset[11]+.025,
         y0=1.02,lty='solid',lwd=2)
segments(x0=start[5]+offset[8]-.025,
         x1=start[5]+offset[11]+.025,
         y0=1.02,lty='solid',lwd=2)

axis(1, c(2,7,12,17,22),
     labels = c("min","max","mean","sd","Chi-2"), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(2, seq(0,1,by=.2),
     seq(0,1,by=.2), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)

dev.off()

# -------
# Chi-2 ----
# -------

chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.obs"))
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.sim"))
# calculations are rowwise
fit.obs=apply(chi2.obs,1,sum)
fit.sim=apply(chi2.sim,1,sum)
p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
mean(p.chi2.calc)

p.pop = matrix(NA,nrow=dim(chi2.obs)[1],ncol=20)
for(i in 1:20){
  tmp.chi2.obs=chi2.obs[,data$site==i]
  tmp.chi2.sim=chi2.sim[,data$site==i]
  fit.obs=apply(tmp.chi2.obs,1,sum)
  fit.sim=apply(tmp.chi2.sim,1,sum)
  p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
  p.pop[,i] = p.chi2.calc
}
apply(p.pop,2,mean)


pdf(file=paste0(outputDirectory,"seedlingSurvivalFruiting-population.pdf"),height=6,width=8)

par(mfrow=c(1,2))
pop.sample = 1:20
plot(y=(pop.sample),x=rev(apply(p.pop,2,mean)),
     xlim=c(0,1),pch=16,
     ylab="Population",xlab="p-Value (chi-squared)",
     main="Fruiting plant counts",    
     axes=FALSE,frame=FALSE,xaxt='n',yaxt='n')
abline(v=c(.1,.9),lty='dotted')
abline(v=c(.2,.8),lty='dotted',col='gray')

axis(2, (1:20),
     labels = rev(siteNames), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(1,  seq(0,1,by=.2), col.ticks = 1)

# population-wide, this seems to pass checks okay

n.iter = dim(chi2.obs)[1]

p.chi.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=length(years))
for(j in 1:20){
  for(i in 1:length(years)){
    tmp.chi2.obs=chi2.obs[,data$site==j& data$year==i]
    tmp.chi2.sim=chi2.sim[,data$site==j& data$year==i]
    fit.obs=apply(tmp.chi2.obs,1,sum)
    fit.sim=apply(tmp.chi2.sim,1,sum)
 #   exclude=fit.obs!=fit.sim
 #   p.chi2.calc=ifelse(fit.sim[exclude]-fit.obs[exclude]>=0,1,0)
      p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
    p.pop[,i] = p.chi2.calc
  }
  p.chi.list[[j]] = apply(p.pop,2,mean,na.rm=TRUE)
}
p.chi.mat=do.call(rbind,p.chi.list)

time.sample = 1:length(years)+2005
plot(time.sample,p.chi.mat[1,],
     ylim=c(0,1),pch=16,xlim=c(2006,2019),
     xlab="Year",ylab="p-Value (chi-squared)",
     main="Fruiting plant counts",type='n')
for(i in 1:20){
  points(time.sample+rnorm(1,0,sd=.05),
         p.chi.mat[i,],pch=19,cex=.5,col='black')
}
abline(h=c(.1,.9),lty='dotted')
abline(h=c(.2,.8),lty='dotted',col='gray')
dev.off()

cbind(time.sample,apply(ifelse(p.chi.mat>.9,1,0),2,sum)/20)

# per year*population shows some values that are not well modeled this way

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

time.sample = 1:length(years)+2005
for(i in 1:20){
  
  tot.seeds=c()
  for(j in 1:length(years)){
    tot.seeds[j]=sum(data$seedlingNumber[data$site==i&data$year==j],na.rm=TRUE)
  }
  
  plot(NA,NA,
       ylim=c(0,1),pch=16,xlim=c(2006,2019),
       ylab='',xlab='',xaxt='n',yaxt='n')
  points(time.sample+rnorm(1,0,sd=.05),
         p.chi.mat[i,],pch=19,cex=1,
         col="black")
  
  abline(h=c(.1,.9),lty='dotted')
  abline(h=c(.2,.8),lty='dotted',col='gray')
  
  text(x=2005.5,y=.04,siteNames[i],pos=4)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
  ifelse(i%in%c(5), legend(x = 15, y = 1,
                           col = c('gray','orange'),
                           lty = c(1,1),
                           legend = c("Persistence only","Persistence & viability"),
                           cex=.55,
                           box.lty=0), NA)
}
mtext("Year", side = 1, outer = TRUE, line = 2.2)
mtext("Bayesian p-value (chi-2)", side = 2, outer = TRUE, line = 2.2)


par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1,
    mgp=c(3,0,0))

for(i in 1:20){

  tot.fruitpl=c()
  for(j in 1:length(years)){
    tot.fruitpl[j]=sum(data$fruitplNumber[data$site==i&data$year==j],na.rm=TRUE)
  }
  
  plot(NA,NA,
       ylim=c(0,1),pch=16,xlim=c(0,max(tot.fruitpl)),
       ylab='',xlab='',yaxt='n',xaxt='n')
  points(tot.fruitpl,
         p.chi.mat[i,],pch=19,cex=1,
         col=ifelse(p.chi.mat[i,]>.9|p.chi.mat[i,]<.1,"purple","black"))
  
  abline(h=c(.1,.9),lty='dotted')
  
  text(x=5,y=.04,siteNames[i],pos=4)
  axis(1,cex.axis=.75,tick=FALSE)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Total number of fruiting plants", side = 1, outer = TRUE, line = 2.2)
mtext("Bayesian p-value (chi-2)", side = 2, outer = TRUE, line = 2.2)


par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1,
    mgp=c(3,0,0))

for(i in 1:20){
  
  var.fruitpl=c()
  for(j in 1:length(years)){
    var.fruitpl[j]=var(data$fruitplNumber[data$site==i&data$year==j],na.rm=TRUE)
  }
  
  plot(NA,NA,
       ylim=c(0,1),pch=16,xlim=c(0,max(var.fruitpl)),
       ylab='',xlab='',yaxt='n',xaxt='n')
  points(var.fruitpl,
         p.chi.mat[i,],pch=19,cex=1,
         col=ifelse(p.chi.mat[i,]>.9|p.chi.mat[i,]<.1,"purple","black"))
  
  abline(h=c(.1,.9),lty='dotted')
  
  text(x=5,y=.04,siteNames[i],pos=4)
  axis(1,cex.axis=.75,tick=FALSE)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Variance in fruiting plant number", side = 1, outer = TRUE, line = 2.2)
mtext("Bayesian p-value (chi-2)", side = 2, outer = TRUE, line = 2.2)

# Years in which almost no plants survived, this model does not capture that


pdf(file=paste0(outputDirectory,"seedlingSurvivalFruiting-populationyear.pdf"),height=6,width=6)


par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1,
    mgp=c(3,0,0))

for(i in 1:20){
  
  tot.seeds=c()
  for(j in 1:length(years)){
    tot.seeds[j]=sum(data$seedlingNumber[data$site==i&data$year==j],na.rm=TRUE)
  }
  
  plot(NA,NA,
       ylim=c(0,1),pch=16,xlim=c(0,max(tot.seeds)),
       ylab='',xlab='',yaxt='n',xaxt='n')
  points(tot.seeds,
         p.chi.mat[i,],pch=19,cex=1,
         col=ifelse(p.chi.mat[i,]>.9|p.chi.mat[i,]<.1,"purple","black"))
  
  abline(h=c(.1,.9),lty='dotted')
  
  text(x=5,y=.04,siteNames[i],pos=4)
  axis(1,cex.axis=.75,tick=FALSE)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Total number of seedlings", side = 1, outer = TRUE, line = 2.2)
mtext("Bayesian p-value (chi-2)", side = 2, outer = TRUE, line = 2.2)
#mtext("Population*year-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)


par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1,
    mgp=c(3,0,0))

for(i in 1:20){
  
  var.fruitpl=c()
  for(j in 1:length(years)){
    var.fruitpl[j]=sum(ifelse(data$seedlingNumber[data$site==i&data$year==j]>0,1,0))
  }
  
  plot(NA,NA,
       ylim=c(0,1),pch=16,xlim=c(0,max(var.fruitpl)),
       ylab='',xlab='',yaxt='n',xaxt='n')
  points(var.fruitpl,
         p.chi.mat[i,],pch=19,cex=1,
         col=ifelse(p.chi.mat[i,]>.9|p.chi.mat[i,]<.1,"purple","black"))
  
  abline(h=c(.1,.9),lty='dotted')
  
  text(x=5,y=.04,siteNames[i],pos=4)
  axis(1,cex.axis=.75,tick=FALSE)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Plots with nonzero numbers of seedlings", side = 1, outer = TRUE, line = 2.2)
mtext("Bayesian p-value (chi-2)", side = 2, outer = TRUE, line = 2.2)
#mtext("Population*year-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)



# years with low number of seedlings (trials) and fruiting plants (successes)
# pull the mean towards the population-level average
# due to pooling


dev.off()


# -------------------------------------------------------------------
# Bayesian p-values: omnibus Chi-squared for counts of fruiting plants
# entire dataset, population-level
# removing data with 0 trials
# -------------------------------------------------------------------
chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.obs"))
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.sim"))
trialsIndex=data$seedlingNumber!=0
chi2.obs=chi2.obs[,trialsIndex]
chi2.sim=chi2.sim[,trialsIndex]

# calculations are rowwise
fit.obs=apply(chi2.obs,1,sum)
fit.sim=apply(chi2.sim,1,sum)
p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
mean(p.chi2.calc)

data2 = list()
data2$site=data$site[trialsIndex]
data2$year=data$year[trialsIndex]

n.iter = dim(chi2.obs)[1]

p.chi.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=length(years))
for(j in 1:20){
  for(i in 1:length(years)){
    tmp.chi2.obs=as.matrix(chi2.obs[,data2$site==j& data2$year==i])
    tmp.chi2.sim=as.matrix(chi2.sim[,data2$site==j& data2$year==i])
    
    fit.obs=    if (sum(tmp.chi2.obs)==0) {NA} else {apply(tmp.chi2.obs,1,sum)}
    fit.sim=    if (sum(tmp.chi2.sim)==0) {NA} else {apply(tmp.chi2.sim,1,sum)}

    #   exclude=fit.obs!=fit.sim
    #   p.chi2.calc=ifelse(fit.sim[exclude]-fit.obs[exclude]>=0,1,0)
    
    p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
    p.pop[,i] = p.chi2.calc
  }
  p.chi.list[[j]] = apply(p.pop,2,mean,na.rm=TRUE)
}
p.chi.mat=do.call(rbind,p.chi.list)

# per year*population shows some values that are not well modeled this way

pdf(file=paste0(outputDirectory,"seedlingSurvivalFruiting-zerotrials.pdf"),height=6,width=6)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

time.sample = 1:length(years)+2005
for(i in 1:20){
  
  plot(NA,NA,
       ylim=c(0,1),pch=16,xlim=c(2006,2019),
       ylab='',xlab='',xaxt='n',yaxt='n')
  
  index=is.na(p.chi.mat[i,])
  if (sum(index)>0) {for(j in 1:length(time.sample[index])){
    ts=time.sample[index]
    polygon(x=c(ts[j]-.5,ts[j]-.5,
                ts[j]+.5, ts[j]+.5),
            y=c(-.1,1.1,1.1,-.1),col='gray90',border=0)
  }} else {NA}
  
  points(time.sample,
         p.chi.mat[i,],pch=19,cex=1,
         col="black")

  
  abline(h=c(.1,.9),lty='dotted')
  abline(h=c(.2,.8),lty='dotted',col='gray')
  
  text(x=2005.5,y=.04,siteNames[i],pos=4)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)

}
mtext("Year", side = 1, outer = TRUE, line = 2.2)
mtext("Bayesian p-value (chi-2)", side = 2, outer = TRUE, line = 2.2)
dev.off()