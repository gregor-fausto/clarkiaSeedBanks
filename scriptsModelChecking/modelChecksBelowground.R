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

mcmcSamples <- readRDS(modelFittingFiles[[grep("seedBurialSamplesChecks",modelFittingFiles)]])
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

# ---
# *Germination ----
# ---

par(mfrow=c(1,3))

y.sim=MCMCchains(mcmcSamples, params = "seedlingJan_sim")

pdf(file="~/Dropbox/clarkiaSeedBanks/products/figures/modelChecks/germination-population.pdf",height=6,width=6)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  index=data$siteGermination==i&data$gIndex==1
  # hist(data$seedlingJan[index]/data$totalJan[index], breaks = 10, 
  #      freq = FALSE, main='',
  #      ylab='',xlab='',xaxt='n',yaxt='n') 
  plot(NA,NA,xlim=c(0,1),ylim=c(0,10), main='',
             ylab='',xlab='',xaxt='n',yaxt='n')
  tmp=sweep(y.sim[,index], 2, data$totalJan[index], FUN = '/')
  
  for(j in 1:20){
    pg=ggplot_build(ggplot(data.frame(d=tmp[j,]))+stat_density(aes(x=d)))
   # plot(pg$data[[1]]$x,pg$data[[1]]$y,type='l')
   # lines(density(tmp[j,]), col = "gray",lwd=.05)
    lines(pg$data[[1]]$x,pg$data[[1]]$y,lwd=0.25,col='gray')
  }
  pg=ggplot_build(ggplot(data.frame(d=data$seedlingJan[index]/data$totalJan[index]))+stat_density(aes(x=d)))
  
  lines(pg$data[[1]]$x,pg$data[[1]]$y,lwd=1,col='black')
  
  text(x=.05,y=9.5,siteNames[i],pos=4)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Probability of germination (1)", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)



par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  index=data$siteGermination==i&data$gIndex==2
  # hist(data$seedlingJan[index]/data$totalJan[index], breaks = 10, 
  #      freq = FALSE, main='',
  #      ylab='',xlab='',xaxt='n',yaxt='n') 
  
  plot(NA,NA,xlim=c(0,1),ylim=c(0,10), main='',
       ylab='',xlab='',xaxt='n',yaxt='n')
  tmp=sweep(y.sim[,index], 2, data$totalJan[index], FUN = '/')
  
  for(j in 1:20){
    pg=ggplot_build(ggplot(data.frame(d=tmp[j,]))+stat_density(aes(x=d)))
   # plot(pg$data[[1]]$x,pg$data[[1]]$y,type='l')
   # lines(density(tmp[j,]), col = "gray",lwd=.05)
    lines(pg$data[[1]]$x,pg$data[[1]]$y,lwd=0.25,col='gray')
  }
  pg=ggplot_build(ggplot(data.frame(d=data$seedlingJan[index]/data$totalJan[index]))+stat_density(aes(x=d)))
  
  lines(pg$data[[1]]$x,pg$data[[1]]$y,lwd=1,col='black')
  
  text(x=.05,y=9.5,siteNames[i],pos=4)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Probability of germination (2)", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)



par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  index=data$siteGermination==i&data$gIndex==3
  # hist(data$seedlingJan[index]/data$totalJan[index], breaks = 10, 
  #      freq = FALSE, main='',
  #      ylab='',xlab='',xaxt='n',yaxt='n') 
  
  plot(NA,NA,xlim=c(0,1),ylim=c(0,10), main='',
       ylab='',xlab='',xaxt='n',yaxt='n')
  tmp=sweep(y.sim[,index], 2, data$totalJan[index], FUN = '/')
  
  for(j in 1:20){
    pg=ggplot_build(ggplot(data.frame(d=tmp[j,]))+stat_density(aes(x=d)))
   # plot(pg$data[[1]]$x,pg$data[[1]]$y,type='l')
   # lines(density(tmp[j,]), col = "gray",lwd=.05)
    lines(pg$data[[1]]$x,pg$data[[1]]$y,lwd=0.25,col='gray')
  }
  pg=ggplot_build(ggplot(data.frame(d=data$seedlingJan[index]/data$totalJan[index]))+stat_density(aes(x=d)))
  
  lines(pg$data[[1]]$x,pg$data[[1]]$y,lwd=1,col='black')
  
  text(x=.05,y=9.5,siteNames[i],pos=4)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Probability of germination (3)", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)

dev.off()

# ---
# *Intact seed counts ----
# ---

par(mfrow=c(2,3))

y.sim=MCMCchains(mcmcSamples, params = "y_sim")

for(i in c(1,7,11)){
  hist(data$y[data$siteSurvival==1&data$compIndex==i], breaks = 10, 
       freq = FALSE, main = "Simulated and real data for germination", 
       xlab = expression(paste("germinant count")), cex.lab = 1.2) 
  y_sim=y.sim[,data$siteSurvival==1&data$compIndex==i]
  hist(y_sim, col = "red",add=TRUE,freq=FALSE)
}

dev.off()



pdf(file="~/Dropbox/clarkiaSeedBanks/products/figures/modelChecks/decay-population.pdf",height=8,width=6)

par(mfrow=c(5,6),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

group = list(1:5,6:10,11:15,16:20)

for(g in 1:4){
  
  group.tmp=group[[g]]
  
  for(j in 1:5){
    
    j.tmp=group.tmp[j]
    
    n.samples = 25
    iter.ind = sample(1:n.iter,n.samples)
    
    for(i in 1:6){
      
      index = data$siteSurvival==j.tmp&data$compIndex==i
      y_sim=y.sim[,index]
      
      tmp=sweep(y_sim, 2, data$seedStart[index], FUN = '/')
      tmp=tmp[iter.ind,]
      
      list.dens=apply(tmp,1,density)
      all.max=max(unlist(lapply(list.dens, "[", "y")))

      plot(NA,NA,
           xlim=c(0,all.max),ylim=c(0,1),
           main=ifelse(j==1,time.sample[i],''),xaxt='n',xlab='',ylab='',yaxt='n')
      for(h in 1:n.samples){
        m.max=max(tmp[h,])
        m.min=min(tmp[h,])
        
        dens.x = density(tmp[h,],from=m.min,to=m.max)
        
        lines(y=dens.x$x,x=dens.x$y,lwd=0.25,col='orange')
      }
      
      p = data$y[index]/data$seedStart[index]
      m.max=max(p)
      m.min=min(p)
      
      dens.x = density(p,from=m.min,to=m.max)
      
      lines(y=dens.x$x,x=dens.x$y)
      
      ifelse(i%in%c(1),axis(2L),NA)
      axis(1L)
      if(i==6) text(x=.05,y=1,siteNames[j.tmp],pos=4) else NA
      
    }
    
    mtext("Probability of persistence", side = 2, outer = TRUE, line = 2.2)
    mtext("Density", side = 1, outer = TRUE, line = 2.2)
  }
}
dev.off()



# ---
# *s0 ----
# ---

y.sim=MCMCchains(mcmcSamples, params = "plotSeedlings_sim")

pdf(file="~/Dropbox/clarkiaSeedBanks/products/figures/modelChecks/s0-population.pdf",height=6,width=6)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  n.samples = 25
  iter.ind = sample(1:n.iter,n.samples)
  
  index=data$sitePlot==i
  y_sim=y.sim[,index]
  
  tmp=sweep(y_sim, 2, data$fec[index], FUN = '/')
  tmp=tmp[iter.ind,]
  
  list.dens=apply(tmp,1,density)
  all.max=max(unlist(lapply(list.dens, "[", "y")))

  plot(NA,NA,
       xlim=c(0,.5),ylim=c(0,all.max),
       main='',xaxt='n',xlab='',ylab='',yaxt='n')
 
   for(h in 1:n.samples){
    m.max=max(tmp[h,])
    m.min=min(tmp[h,])
    
    dens.x = density(tmp[h,],from=m.min,to=m.max)
    
    lines(x=dens.x$x,y=dens.x$y,lwd=0.25,col='orange')
  }
  
  p = data$plotSeedlings[index]/data$fec[index]
  m.max=max(p)
  m.min=min(p)
  
  dens.x = density(p,from=m.min,to=m.max)
  
  lines(x=dens.x$x,y=dens.x$y)
  
  text(x=.3,y=all.max-1.5,siteNames[i],pos=4)
#  text(x=.7,y=all.max-4,sum(data$fec[index]),pos=4)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Probability of seedling emerging from seed in plot", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)
dev.off()


# ---
# Posterior predictive checks ----
# ---

# ---
# *Germination ----
# ---

chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.obs"))
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.sim"))
# calculations are rowwise
fit.obs=apply(chi2.obs,1,sum)
fit.sim=apply(chi2.sim,1,sum)
p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
mean(p.chi2.calc)

p.pop = matrix(NA,nrow=dim(chi2.obs)[1],ncol=20)
for(i in 1:20){
  tmp.chi2.obs=chi2.obs[,data$siteGermination==i]
  tmp.chi2.sim=chi2.sim[,data$siteGermination==i]
  fit.obs=apply(tmp.chi2.obs,1,sum)
  fit.sim=apply(tmp.chi2.sim,1,sum)
  p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
  p.pop[,i] = p.chi2.calc
}
apply(p.pop,2,mean)


par(mfrow=c(1,2))
pop.sample = 1:20
plot(y=(pop.sample),x=rev(apply(p.pop,2,mean)),
     xlim=c(0,1),pch=16,
     ylab="Population",xlab="p-Value",
     main="Germination",    
     axes=FALSE,frame=FALSE,xaxt='n',yaxt='n')
abline(v=c(.1,.9),lty='dotted')
abline(v=c(.2,.8),lty='dotted',col='gray')

axis(2, (1:20),
     labels = rev(siteNames), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(1,  seq(0,1,by=.2), col.ticks = 1)

n.iter = dim(chi2.obs)[1]

p.chi.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=3)
for(j in 1:20){
  for(i in 1:3){
    tmp.chi2.obs=chi2.obs[,data$siteGermination==j&data$gIndex==i]
    tmp.chi2.sim=chi2.sim[,data$siteGermination==j&data$gIndex==i]
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

age.sample = 1:3
plot(age.sample,p.chi.mat[1,],
     ylim=c(0,1),pch=16,xlim=c(1,3),
     xlab="Year",ylab="p-Value",
     main="Germination",type='n')
for(i in 1:20){
  points(age.sample+rnorm(1,0,sd=.05),
         p.chi.mat[i,],pch=19,cex=.5,col='black')
}
abline(h=c(.1,.9),lty='dotted')
abline(h=c(.2,.8),lty='dotted',col='gray')



y.sim=MCMCchains(mcmcSamples,params=c("seedlingJan_sim"))
y.obs=data$seedlingJan
test.sim=apply(y.sim,1,mean)
test.obs=mean(y.obs)
p.test = ifelse(test.sim-test.obs>=0,1,0)
mean(p.test)

n.iter=dim(y.sim)[1]

p.test.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=3)
for(j in 1:20){
  for(i in 1:3){
    tmp.obs=y.obs[data$siteGermination==j&data$gIndex==i]
    tmp.sim=y.sim[,data$siteGermination==j&data$gIndex==i]
    test.obs=mean(tmp.obs)
    test.sim=apply(tmp.sim,1,mean)
    p.test.calc=ifelse(test.sim-test.obs>=0,1,0)
    p.pop[,i] = p.test.calc
  }
  p.test.list[[j]] = apply(p.pop,2,mean,na.rm=TRUE)
}
p.test.mat=do.call(rbind,p.test.list)


par(mfrow=c(1,1))
age.sample = 1:3
plot(age.sample,p.test.mat[1,],
     ylim=c(0,1),pch=16,xlim=c(0,4),
     xlab="Months",ylab="p-Value",
     main="Germinant counts",type='n')
for(i in 1:20){
  points(age.sample+rnorm(1,0,sd=.05),
         p.test.mat[i,],pch=1,cex=.5)
}
abline(h=c(.1,.9),lty='dotted')



f=function(y.sim=chains,y.obs=data,n.obs=data2,model.fun=mean){
  
  n.iter=dim(y.sim)[1]
  
  p.test.list = list()
  p.pop = matrix(NA,nrow=n.iter,ncol=3)
  for(j in 1:20){
    for(i in 1:3){
      index=data$siteGermination==j&data$gIndex==i
      tmp.obs=y.obs[index]/n.obs[index]
      tmp.sim=sweep(y.sim[,index], 2, n.obs[index], FUN = '/')
      test.obs=model.fun(tmp.obs)
      test.sim=apply(tmp.sim,1,model.fun)
      p.test.calc=ifelse(test.sim-test.obs>=0,1,0)
      p.pop[,i] = p.test.calc
    }
    p.test.list[[j]] = apply(p.pop,2,mean,na.rm=TRUE)
  }
  p.test.mat=do.call(rbind,p.test.list)
  return(p.test.mat)
}

sims=MCMCchains(mcmcSamples,params=c("seedlingJan_sim"))
df=data$seedlingJan
df2=data$totalJan

germ.min=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=min)
germ.max=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=max)
germ.mean=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=mean)
germ.sd=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=sd)

pdf(file="~/Dropbox/clarkiaSeedBanks/products/figures/modelChecks/germination-ppc.pdf",height=6,width=6)

par(mfrow=c(1,1))
age.sample = 1:3
plot(age.sample,germ.sd[1,],
     ylim=c(0,1),pch=16,xlim=c(.5,5.5),
     xlab="Months",ylab="p-Value",
     main="Germinant counts",type='n',
     xaxt='n',yaxt='n')
abline(h=0.5,lty='dotted')

col.vec = c("white","gray","black")
offset = c(-.25,0,.25)
  for(j in 1:3){
  points(1+offset[j]+rnorm(20,0,sd=.025),
         germ.min[,j],pch=21,cex=.5,bg=col.vec[j])
    
    points(2+offset[j]+rnorm(20,0,sd=.025),
           germ.max[,j],pch=21,cex=.5,bg=col.vec[j])
    
    points(3+offset[j]+rnorm(20,0,sd=.025),
           germ.mean[,j],pch=21,cex=.5,bg=col.vec[j])
    
    points(4+offset[j]+rnorm(20,0,sd=.025),
           germ.sd[,j],pch=21,cex=.5,bg=col.vec[j])
    
    points(5+offset[j]+rnorm(20,0,sd=.025),
           p.chi.mat[,j],pch=21,cex=.5,bg=col.vec[j])
  }
abline(v=c(1.5,2.5,3.5,4.5))
axis(1, (1:5),
     labels = c("min","max","mean","sd","Chi-2"), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(2, seq(0,1,by=.2),
     seq(0,1,by=.2), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
legend(x=.5,y=.1,legend=c("Age 1", "Age 2", "Age 3"),
       pch =c(21,21,21), cex =.5, pt.bg=col.vec)


f2=function(y.sim=chains,y.obs=data,n.obs=data2,age=1,model.fun=mean){
  
  n.iter=dim(y.sim)[1]
  
  test.list = list()
  for(j in 1:20){
      index=data$siteGermination==j&data$gIndex==age
      tmp.obs=y.obs[index]/n.obs[index]
      tmp.sim=sweep(y.sim[,index], 2, n.obs[index], FUN = '/')
      test.obs=model.fun(tmp.obs)
      test.sim=apply(tmp.sim,1,model.fun)
      test.list[[j]] = list(test.obs,test.sim)
  }
  return(test.list)
}

sims=MCMCchains(mcmcSamples,params=c("seedlingJan_sim"))
df=data$seedlingJan
df2=data$totalJan

germ.min=f2(y.sim=sims,y.obs=df,n.obs=df2,age=1,model.fun=min)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(1,0,1,1) + 0.1)

for(i in 1:20){
  tmp = germ.min[[i]]
  hist( tmp[[2]], breaks=100,main='',freq=FALSE,ylab='',yaxt='n')
  abline(v=tmp[[1]],col='red',lwd=2,lty='dotted')
  y.max = max(hist( tmp[[2]], breaks=100,plot=FALSE)$density)*.75
  x.max = max(hist( tmp[[2]], breaks=100,plot=FALSE)$mids)*.5
  
  text(x=x.max,y=y.max,
       siteNames[i],pos=4)
}

mtext("Probability of germination (1)", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)
dev.off()

# ---
# *Intact seed counts ----
# ---

chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.yobs"))
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.ysim"))

n.iter = dim(chi2.obs)[1]

p.chi.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=6)
for(j in 1:20){
  
  classes = list(c(1,7,11),c(2,8,12),c(3,9),c(4,10),c(5),c(6))
  
  for(i in 1:6){
    cIndex=classes[[i]]
    index = data$siteSurvival==j&data$compIndex%in%cIndex
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
  p.pop = matrix(NA,nrow=n.iter,ncol=6)
  for(j in 1:20){
    
    classes = list(c(1,7,11),c(2,8,12),c(3,9),c(4,10),c(5),c(6))
    
    for(i in 1:6){
      cIndex=classes[[i]]
      index = data$siteSurvival==j&data$compIndex%in%cIndex
      
      tmp.obs=y.obs[index]/n.obs[index]
      tmp.sim=sweep(y.sim[,index], 2, n.obs[index], FUN = '/')
      test.obs=model.fun(tmp.obs)
      test.sim=apply(tmp.sim,1,model.fun)
      p.test.calc=ifelse(test.sim-test.obs>=0,1,0)
      p.pop[,i] = p.test.calc
    }
    p.test.list[[j]] = apply(p.pop,2,mean,na.rm=TRUE)
  }
  p.test.mat=do.call(rbind,p.test.list)
  return(p.test.mat)
}

sims=MCMCchains(mcmcSamples,params=c("y_sim"))
df=data$y
df2=data$seedStart

seeds.min=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=min)
seeds.max=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=max)
seeds.mean=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=mean)
seeds.sd=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=sd)

colfunc <- colorRampPalette(c("white", "black"))
col.vec=colfunc(6)

pdf(file="~/Dropbox/clarkiaSeedBanks/products/figures/modelChecks/decay-ppc.pdf",height=6,width=6)

par(mfrow=c(1,1))
time.sample = 1:6
plot(NA,NA,
     ylim=c(0,1),pch=16,xlim=c(.5,13),
     xlab="Months",ylab="p-Value",
     main="Intact seed counts",type='n',
     xaxt='n',yaxt='n')
abline(h=0.5,lty='dotted')

start = c(1,3.5,6,8.5,11)
offset = seq(0,1.5,length.out=6)
for(j in 1:6){
  points(start[1]+offset[j]+rnorm(20,0,sd=.025),
         seeds.min[,j],pch=21,cex=.5,bg=col.vec[j])

  points(start[2]+offset[j]+rnorm(20,0,sd=.025),
         seeds.max[,j],pch=21,cex=.5,bg=col.vec[j])
  
  points(start[3]+offset[j]+rnorm(20,0,sd=.025),
         seeds.mean[,j],pch=21,cex=.5,bg=col.vec[j])
  
  points(start[4]+offset[j]+rnorm(20,0,sd=.025),
         seeds.sd[,j],pch=21,cex=.5,bg=col.vec[j])
  
  points(start[5]+offset[j]+rnorm(20,0,sd=.025),
         p.chi.mat[,j],pch=21,cex=.5,bg=col.vec[j])
}
abline(v=c(3,5.5,8,10.5))
axis(1, (start+.75),
     labels = c("min","max","mean","sd","Chi-2"), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(2, seq(0,1,by=.2),
     seq(0,1,by=.2), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
legend(x=11.5,y=1,legend=paste0(c(4,12,16,24,28,36)),
       title="Time (months)",
       pch=21, cex =.5, pt.bg=col.vec)

f2=function(y.sim=chains,y.obs=data,n.obs=data2,class=4,model.fun=mean){
  
  n.iter=dim(y.sim)[1]
  
  test.list = list()
  for(j in 1:20){
    
    index = data$siteSurvival==j&data$compIndex==class
    
    tmp.obs=y.obs[index]/n.obs[index]
    tmp.sim=sweep(y.sim[,index], 2, n.obs[index], FUN = '/')
    test.obs=model.fun(tmp.obs)
    test.sim=apply(tmp.sim,1,model.fun)
    test.list[[j]] = list(test.obs,test.sim)

  }
  return(test.list)
}

sims=MCMCchains(mcmcSamples,params=c("y_sim"))
df=data$y
df2=data$seedStart

seeds.min=f2(y.sim=sims,y.obs=df,n.obs=df2,class=1,model.fun=max)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(1,0,1,1) + 0.1)

for(i in 1:20){
  tmp = seeds.min[[i]]
  hist( tmp[[2]], breaks=100,main='',freq=FALSE,ylab='',yaxt='n')
  abline(v=tmp[[1]],col='red',lwd=2,lty='dotted')
  y.max = max(hist( tmp[[2]], breaks=100,plot=FALSE)$density)*.75
  x.max = min(hist( tmp[[2]], breaks=100,plot=FALSE)$mids)*1
  
  text(x=x.max,y=y.max,
       siteNames[i],pos=4)
}

mtext("Maximum (first January)", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)
dev.off()


# ---
# *s0 ----
# ---

chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.plot.obs"))
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.plot.sim"))

n.iter = dim(chi2.obs)[1]

p.chi.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=2)
for(j in 1:20){
  
  for(i in 1:2){
    index = data$sitePlot==j&data$yearPlot==i
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
  p.pop = matrix(NA,nrow=n.iter,ncol=2)
  for(j in 1:20){
    
    for(i in 1:2){
      index = data$sitePlot==j&data$yearPlot==i
      tmp.obs=y.obs[index]/n.obs[index]
      tmp.sim=sweep(y.sim[,index], 2, n.obs[index], FUN = '/')
      test.obs=model.fun(tmp.obs)
      test.sim=apply(tmp.sim,1,model.fun)
      p.test.calc=ifelse(test.sim-test.obs>=0,1,0)
      p.pop[,i] = p.test.calc
    }
    p.test.list[[j]] = apply(p.pop,2,mean,na.rm=TRUE)
  }
  p.test.mat=do.call(rbind,p.test.list)
  return(p.test.mat)
}

sims=MCMCchains(mcmcSamples,params=c("plotSeedlings_sim"))
df=data$plotSeedlings
df2=data$fec

seeds.min=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=min)
seeds.max=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=max)
seeds.mean=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=mean)
seeds.sd=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=sd)

colfunc <- colorRampPalette(c("white", "black"))
col.vec=colfunc(2)

pdf(file="~/Dropbox/clarkiaSeedBanks/products/figures/modelChecks/s0-ppc.pdf",height=6,width=6)

par(mfrow=c(1,1))
time.sample = 1:2
plot(NA,NA,
     ylim=c(0,1),pch=16,xlim=c(.5,5.5),
     xlab="",ylab="",
     main="Emerging seedlings",type='n',
     xaxt='n',yaxt='n')
abline(h=0.5,lty='dotted')

start = c(1,2,3,4,5)
offset = c(-.125,.125)
for(j in 1:2){
  points(start[1]+offset[j]+rnorm(20,0,sd=.025),
         seeds.min[,j],pch=21,cex=.5,bg=col.vec[j])
  
  points(start[2]+offset[j]+rnorm(20,0,sd=.025),
         seeds.max[,j],pch=21,cex=.5,bg=col.vec[j])
  
  points(start[3]+offset[j]+rnorm(20,0,sd=.025),
         seeds.mean[,j],pch=21,cex=.5,bg=col.vec[j])
  
  points(start[4]+offset[j]+rnorm(20,0,sd=.025),
         seeds.sd[,j],pch=21,cex=.5,bg=col.vec[j])
  
  points(start[5]+offset[j]+rnorm(20,0,sd=.025),
         p.chi.mat[,j],pch=21,cex=.5,bg=col.vec[j])
}
abline(v=c(1.5,2.5,3.5,4.5))
axis(1, (1:5),
     labels = c("min","max","mean","sd","Chi-2"), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(2, seq(0,1,by=.2),
     seq(0,1,by=.2), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
legend(x=.5,y=.1,legend=paste0(c(2008,2009)),
       title="Year",
       pch =21, cex =.5, pt.bg=col.vec)
dev.off()

f2=function(y.sim=chains,y.obs=data,n.obs=data2,class=1,model.fun=mean){
  
  n.iter=dim(y.sim)[1]
  
  test.list = list()
  for(j in 1:20){
    
    index = data$sitePlot==j&data$yearPlot==class
    
    tmp.obs=y.obs[index]/n.obs[index]
    tmp.sim=sweep(y.sim[,index], 2, n.obs[index], FUN = '/')
    test.obs=model.fun(tmp.obs)
    test.sim=apply(tmp.sim,1,model.fun)
    test.list[[j]] = list(test.obs,test.sim)
    
  }
  return(test.list)
}

sims=MCMCchains(mcmcSamples,params=c("plotSeedlings_sim"))
df=data$plotSeedlings
df2=data$fec

seeds.min=f2(y.sim=sims,y.obs=df,n.obs=df2,class=2,model.fun=min)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(1,0,1,1) + 0.1)

for(i in 1:20){
  tmp = seeds.min[[i]]
  mat=ifelse(is.finite(tmp[[2]]),tmp[[2]],NA)
  if(is.na(mat)[1]){
    plot(0,type='n',axes=FALSE,ann=FALSE,xlim=c(0,1),ylim=c(0,1))
    text(x=0.05,y=.95,
         siteNames[i],pos=4)
  } else {
      hist( mat, breaks=100,main='',freq=FALSE,ylab='',yaxt='n')
    abline(v=tmp[[1]],col='red',lwd=2,lty='dotted')
    y.max = max(hist( mat, breaks=100,plot=FALSE)$density)*.75
    x.max = min(hist( mat, breaks=100,plot=FALSE)$mids)*1
    
    text(x=x.max,y=y.max,
         siteNames[i],pos=4)
    }

}

mtext("Minimum (2008)", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)
dev.off()


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Model comparison
# -------------------------------------------------------------------
# -------------------------------------------------------------------

LLmat.wb <- MCMCchains(mcmcSamples,params="logLik_y")
rel_n_eff.wb <- relative_eff(exp(LLmat.wb), chain_id = rep(1, each = n.iter))
looWb <- loo(LLmat.wb, r_eff = rel_n_eff.wb, cores = 2)

print(looWb)
plot(looWb)


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Model comparison
# -------------------------------------------------------------------
# -------------------------------------------------------------------


# WEIBULL

a.wb=MCMCchains(mcmcSamples,params="a")
b.wb=MCMCchains(mcmcSamples,params="mu_s")
b0.wb=MCMCchains(mcmcSamples,params="mu0_s")


  inv.b0=cbind(exp(-b0.wb/a.wb))
  inv.b0.parm=apply(inv.b0,2,quantile,probs=c(0.025,.5,.975))
  
  inv.b=cbind(exp(-b.wb[,1]/a.wb),exp(-b.wb[,2]/a.wb),exp(-b.wb[,3]/a.wb))
  inv.b.parm=apply(inv.b,2,quantile,probs=c(0.025,.5,.975))
  
  a.parm=apply(a.wb,2,quantile,probs=c(0.025,.5,.975))


  par(mfrow = c(4,5),
      oma = c(5,4,0,0) + 0.1,
      mar = c(0,0,1,1) + 0.1)

  t=seq(0,max(data$months),.01)
# 
# par(mfrow=c(1,1))
# plot(t,exp(-(t/inv.b0.parm[2,1])^a.parm[2,1]),
#      xlim=c(0,max(data$months)),ylim=c(0,1),type='l')
# lines(t,exp(-(t/inv.b0.parm[1,1])^a.parm[2,1]),lty='dotted')
# lines(t,exp(-(t/inv.b0.parm[3,1])^a.parm[2,1]),lty='dotted')
# 
# 
# # sample from the same rows of the posterior to keep correlations between parameters
# # par(mfrow=c(1,3))
# # plot(a.wb,b.wb[,1])
# # plot(a.wb,b.wb[,2])
# # plot(a.wb,b.wb[,3])
# 
# par(mfrow=c(1,3))
## AGE 1
  for(j in 1:20){
plot(t,exp(-(t/inv.b.parm[2,j])^a.parm[2,j]),
     xlim=c(0,max(data$months)),ylim=c(0,1),type='l')
for(i in 1:100){
  tmp=sample(1:dim(a.wb)[1],1)
  lines(t,exp(-(t/inv.b[tmp,j])^a.wb[tmp,j]),
        lwd=.5,col=ifelse(a.wb[tmp,j]>1,"purple","orange"))
}
lines(t,exp(-(t/inv.b.parm[2,j])^a.parm[2,j]),lwd=2)

# df<-survivalData[survivalData$yearSurvival==2006,]
# points(df$months,df$y/df$seedStart,pch=16,cex=.75)
}
## AGE 2
plot(t,exp(-(t/inv.b.parm[2,2])^a.parm[2,1]),
     xlim=c(0,max(data$months)),ylim=c(0,1),type='l')
for(i in 1:100){
  tmp=sample(1:dim(a.wb)[1],1)
  lines(t,exp(-(t/inv.b[tmp,2])^a.wb[tmp,1]),
        lwd=.5,col=ifelse(a.wb[tmp,1]>1,"purple","orange"))
}
lines(t,exp(-(t/inv.b.parm[2,2])^a.parm[2,1]),lwd=2)


df<-survivalData[survivalData$yearSurvival==2007,]
points(df$months,df$y/df$seedStart,pch=16,cex=.75)


## AGE 3
plot(t,exp(-(t/inv.b.parm[2,3])^a.parm[2,1]),
     xlim=c(0,max(data$months)),ylim=c(0,1),type='l')
for(i in 1:100){
  tmp=sample(1:dim(a.wb)[1],1)
  lines(t,exp(-(t/inv.b[tmp,3])^a.wb[tmp,1]),
        lwd=.5,col=ifelse(a.wb[tmp,1]>1,"purple","orange"))
}
lines(t,exp(-(t/inv.b.parm[2,3])^a.parm[2,1]),lwd=2)

df<-survivalData[survivalData$yearSurvival==2008,]
points(df$months,df$y/df$seedStart,pch=16,cex=.75)



# -------------------------------------------------------------------
# Binned residual plots
# see check 6: https://moodle2.units.it/pluginfile.php/290133/mod_resource/content/1/model_check_script.R
# https://discourse.mc-stan.org/t/posterior-prediction-from-logit-regression/12217
# -------------------------------------------------------------------
# 
# sims <- MCMCchains(mcmcSamples, params = "y_sim")
# sims.subset=sims[sample(dim(sims)[1],10000),]
# 
# sims <- MCMCchains(mcmcSamples, params = "seedlingJan_sim")
# 
# for(i in 1:6){
#   
#   par(mfrow=c(2,3))
#   bayesplot::ppc_error_binned(data$seedlingJan[data$gIndex==3], 
#                               sims[1:6,data$gIndex==3])
#   
# }
