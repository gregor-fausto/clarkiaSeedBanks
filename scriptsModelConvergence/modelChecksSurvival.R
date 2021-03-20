rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)

library(MCMCvis)
library(tidybayes)
library(tidyverse)
library(magrittr)
library(bayesplot)
library(rethinking)

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedlingSurvival/"
modelFittingFiles <- paste0(directory,list.files(directory))

mcmcSamples <- readRDS(modelFittingFiles[[grep("seedSurvivalSamplesChecks.rds",modelFittingFiles)]])
data <- readRDS(modelFittingFiles[[grep("data.rds",modelFittingFiles)]])


censusSeedlingsFruitingPlants <- readRDS("~/Dropbox/dataLibrary/postProcessingData/censusSeedlingsFruitingPlants.RDS")

siteNames = unique(censusSeedlingsFruitingPlants$site)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Convergence diagnostics
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# MCMCsummary(mcmcSamples, params = c("mu0"))
# MCMCsummary(mcmcSamples, params = c("sigma0"))
# MCMCsummary(mcmcSamples, params = c("sigma"))

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
year=2006:2019


# par(mfrow = c(1,1),
#     oma = c(2,2,0,0) + 0.1,
#     mar = c(1,0,1,1) + 0.1)
# 
# for(h in 1:20){
# site.index = data$site==h
# total.n=length(y.obs[site.index])
# 
# counter=0
# 
# plot(NA,NA,xlim=c(0,140),ylim=c(0,1),
#      xlab='',ylab='',yaxt='n',xaxt='n',
#      axes=FALSE,frame=FALSE,main=siteNames[h])
# 
# for(j in 1:14){
#   index=data$site==h&data$year==j
#   p.obs=y.obs[index]/n.obs[index]
#   p.sim=sweep(y.sim[,index], 2, n.obs[index], FUN = '/')
#   p.sim.sum=apply(p.sim,2,quantile,probs=c(.025,.5,.975),na.rm=TRUE)
#   
#   subset=sample(1:length(p.obs),10)
#   
#   order.index=order(p.obs[subset])
#   year.n = length(p.obs[order.index])-1
#   for(i in 1:dim(p.sim.sum)[2]){
#     segments(x0=counter:(counter+year.n), 
#              y0=p.sim.sum[1,order.index],
#              y1=p.sim.sum[3,order.index])
#     points(x=counter:(counter+year.n), 
#            y=p.sim.sum[2,order.index],
#            pch=21,bg='white',cex=1)
#     points(x=counter:(counter+year.n), 
#            y=p.obs[order.index],
#            pch=16,cex=.75,col='purple')
#   }
#  # abline(v=counter+year.n+.5,lty='dotted',col='gray')
#   if(j %% 2 != 0) {
#     polygon(x=c(counter+year.n+.5,counter+year.n+.5,
#                 counter+10+year.n+.5,counter+10+year.n+.5),
#             y=c(-1,1,1,-1),col='gray90',border=0) }
#   axis(1, counter+year.n/2,
#        labels = year[j], las = 1, 
#        col = NA, col.ticks = 1, cex.axis = 1)
#   counter=counter+year.n+1
# }
# axis(2,  seq(0,1,by=.2), col.ticks = 1)
# }
# dev.off()
# #   text(x=3,y=-.9,siteNames[i])
# #   ifelse(i%in%c(16:20),axis(1L),NA)
# #   ifelse(i%in%c(1,6,11,16),axis(2L),NA)





for(h in 1:20){
  site.index = data$site==h
  total.n=length(y.obs[site.index])
  
 # counter=0

  #for(j in 1:14){
    index=data$site==h#&data$year==j
    p.obs=y.obs[index]#/n.obs[index]
   # p.sim=sweep(y.sim[,index], 2, n.obs[index], FUN = '/')
    p.sim=y.sim[,index]
    p.sim.sum=apply(p.sim,2,quantile,probs=c(.025,.5,.975),na.rm=TRUE)
    
    plot(NA,NA,xlim=c(0,max(p.sim)),ylim=c(0,max(p.sim)),
         xlab='',ylab='',yaxt='n',xaxt='n')
   # subset=sample(1:length(p.obs),10)
    
    order.index=order(p.obs)
 #   year.n = length(p.obs[order.index])-1
      segments(x0=p.obs[order.index], 
               y0=p.sim.sum[1,order.index],
               y1=p.sim.sum[3,order.index])
      points(x=p.obs[order.index], 
             y=p.sim.sum[2,order.index],pch=16)
      abline(a=0,b=1)
    # 
    # abline(v=counter+year.n+.5,lty='dotted',col='gray')
    # ifelse(h%in%c(8:10),axis(1, counter+year.n/2,
    #                          labels = year[j], las = 1, 
    #                          col = NA, col.ticks = 1, cex.axis = 1),NA)
    # counter=counter+year.n+1
  }
  #ifelse(h%in%c(1,6),axis(2,  seq(0,1,by=.2), col.ticks = 1),NA)

pdf(file="~/Dropbox/clarkiaSeedBanks/products/figures/modelChecks/seedlingSurvivalFruiting-ppc-population.pdf",height=6,width=6)

for(i in 1:14){
  par(mfrow = c(4,5),
      oma = c(4,5,0,0) + 0.1,
      mar = c(1,0,1,1) + 0.1)
  for(j in 1:20){
    n.samples = 100
    iter.ind = sample(1:n.iter,n.samples)
    
    index=data$site==j&data$year==i
    p.obs=y.obs[index]#/n.obs[index]
    # p.sim=sweep(y.sim[,index], 2, n.obs[index], FUN = '/')
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
    
    p = data$fruitplNumber[index]#/data$seedlingNumber[index]
    m.max=max(p,na.rm=TRUE)
    m.min=min(p,na.rm=TRUE)
    
    dens.x = density(p,from=m.min,to=m.max,na.rm=TRUE)
    
    lines(x=dens.x$x,y=dens.x$y)
    
    ifelse(i%in%c(1),axis(2L),NA)
    
    axis(2,cex=.25,tick=FALSE,line=-1)
    axis(1,cex=.25,tick=FALSE,line=-1)
    
    text(x=all.max.x*.7,y=all.max.y*.9,siteNames[j],pos=4)
    #  ifelse(j%in%c(16:20),axis(1L),NA)
    #ifelse(j%in%c(1,6,11,16),axis(2L),NA)
  }
  
  mtext(paste0("Probability of seedling survival (",year[i],")"), side = 1, outer = TRUE, line = 1.5)
  mtext("Density", side = 2, outer = TRUE, line = 2.2)
}
dev.off()



# -------------------------------------------------------------------
# Fruiting plants
# -------------------------------------------------------------------
# par(mfrow=c(1,3))
# 
# for(i in 1:20){
# 
#   hist(data$fruitplNumber[data$site==1], breaks = 100,
#        freq = FALSE, main = "Simulated and real data for germination",
#        xlab = expression(paste("germinant count")), cex.lab = 1.2)
#   fruitplNumber_sim=MCMCchains(mcmcSamples, params = "fruitplNumber_sim")[,data$site==1]
#   lines(density(fruitplNumber_sim,adjust=10), col = "red")
# }
# 
# par(mfrow = c(4,5),
#     oma = c(5,4,0,0) + 0.1,
#     mar = c(0,0,1,1) + 0.1)
# for(i in 1:20){
#   
#   hist(data$fruitplNumber[data$site==i], breaks = 100,
#        freq = FALSE, main = "Simulated and real data for germination",
#        xlab = expression(paste("germinant count")), cex.lab = 1.2)
#   fruitplNumber_sim=MCMCchains(mcmcSamples, params = "fruitplNumber_sim")[,data$site==i]
#   lines(density(fruitplNumber_sim,adjust=10), col = "red")
# }

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Posterior Predictive Checks
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Binned residual plots
# see check 6: https://moodle2.units.it/pluginfile.php/290133/mod_resource/content/1/model_check_script.R
# https://discourse.mc-stan.org/t/posterior-prediction-from-logit-regression/12217
# -------------------------------------------------------------------
# 
# sims <- MCMCchains(mcmcSamples, params = "fruitplNumber_sim")
# #sims.subset=sims[sample(dim(sims)[1],10000),]
# 
# #sims <- MCMCchains(mcmcSamples, params = "seedlingJan_sim")
# 
# 
# par(mfrow=c(2,3))
# bayesplot::ppc_error_binned(data$fruitplNumber[data$site==1],
#                             sims[1:6,data$site==1])
# bayesplot::ppc_error_binned(data$fruitplNumber[data$site==2],
#                             sims[1:6,data$site==2])
# 
# par(mfrow=c(3,5))
# 
# for(i in 1:14){
# hist(apply(sims[,data$site==1&data$year==i],2,mean)/data$seedlingNumber[data$site==1&data$year==i],breaks=100,main='')
# abline(v=mean(data$fruitplNumber[data$site==1&data$year==i]/data$seedlingNumber[data$site==1&data$year==i],na.rm=TRUE),col='red',lwd=2)
# }

# ---
# Test statistics ----
# ---

chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.obs"))
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.sim"))

n.iter = dim(chi2.obs)[1]

p.chi.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=14)
for(j in 1:20){
  
  for(i in 1:14){
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
  p.pop = matrix(NA,nrow=n.iter,ncol=14)
  for(j in 1:20){
    
    for(i in 1:14){
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
col.vec=colfunc(14)

pdf(file="~/Dropbox/clarkiaSeedBanks/products/figures/modelChecks/seedlingSurvivalFruiting-pvals.pdf",height=6,width=6)

par(mfrow=c(1,1))
time.sample = 1:14
plot(NA,NA,
     ylim=c(0,1),pch=16,xlim=c(-.5,24.5),
     xlab="",ylab="",
     main="Seedling survival to fruiting",type='n',
     xaxt='n',yaxt='n',frame=FALSE)

start = c(0,5,10,15,20)
offset = seq(0,4,length.out=14)#c(-.125,.125)
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

for(j in 1:14){

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
# 
# f2=function(y.sim=chains,y.obs=data,n.obs=data2,class=1,model.fun=mean){
#   
#   n.iter=dim(y.sim)[1]
#   
#   test.list = list()
#   for(j in 1:20){
#     
#     index = data$sitePlot==j&data$yearPlot==class
#     
#     tmp.obs=y.obs[index]/n.obs[index]
#     tmp.sim=sweep(y.sim[,index], 2, n.obs[index], FUN = '/')
#     test.obs=model.fun(tmp.obs)
#     test.sim=apply(tmp.sim,1,model.fun)
#     test.list[[j]] = list(test.obs,test.sim)
#     
#   }
#   return(test.list)
# }
# 
# sims=MCMCchains(mcmcSamples,params=c("plotSeedlings_sim"))
# df=data$plotSeedlings
# df2=data$fec
# 
# seeds.min=f2(y.sim=sims,y.obs=df,n.obs=df2,class=2,model.fun=min)
# 
# par(mfrow = c(4,5),
#     oma = c(5,4,0,0) + 0.1,
#     mar = c(1,0,1,1) + 0.1)
# 
# for(i in 1:20){
#   tmp = seeds.min[[i]]
#   mat=ifelse(is.finite(tmp[[2]]),tmp[[2]],NA)
#   if(is.na(mat)[1]){
#     plot(0,type='n',axes=FALSE,ann=FALSE,xlim=c(0,1),ylim=c(0,1))
#     text(x=0.05,y=.95,
#          siteNames[i],pos=4)
#   } else {
#     hist( mat, breaks=100,main='',freq=FALSE,ylab='',yaxt='n')
#     abline(v=tmp[[1]],col='red',lwd=2,lty='dotted')
#     y.max = max(hist( mat, breaks=100,plot=FALSE)$density)*.75
#     x.max = min(hist( mat, breaks=100,plot=FALSE)$mids)*1
#     
#     text(x=x.max,y=y.max,
#          siteNames[i],pos=4)
#   }
#   
# }
# 
# mtext("Minimum (2008)", side = 1, outer = TRUE, line = 2.2)
# mtext("Density", side = 2, outer = TRUE, line = 2.2)
# #dev.off()


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


pdf(file="~/Dropbox/clarkiaSeedBanks/products/figures/modelChecks/seedlingSurvivalFruiting-population.pdf",height=6,width=8)

par(mfrow=c(1,2))
pop.sample = 1:20
plot(y=(pop.sample),x=rev(apply(p.pop,2,mean)),
     xlim=c(0,1),pch=16,
     ylab="Population",xlab="p-Value",
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
p.pop = matrix(NA,nrow=n.iter,ncol=14)
for(j in 1:20){
  for(i in 1:14){
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

time.sample = 1:14+2005
plot(time.sample,p.chi.mat[1,],
     ylim=c(0,1),pch=16,xlim=c(2006,2019),
     xlab="Year",ylab="p-Value",
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

time.sample = 1:14+2005
for(i in 1:20){
  
  tot.seeds=c()
  for(j in 1:14){
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
#mtext("Population*year-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)


par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1,
    mgp=c(3,0,0))

for(i in 1:20){

  tot.fruitpl=c()
  for(j in 1:14){
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
mtext("Population*year-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)


par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1,
    mgp=c(3,0,0))

for(i in 1:20){
  
  var.fruitpl=c()
  for(j in 1:14){
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
mtext("Population*year-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)

# Years in which almost no plants survived, this model does not capture that


pdf(file="~/Dropbox/clarkiaSeedBanks/products/figures/modelChecks/seedlingSurvivalFruiting-populationyear.pdf",height=6,width=6)


par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1,
    mgp=c(3,0,0))

for(i in 1:20){
  
  tot.seeds=c()
  for(j in 1:14){
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
  for(j in 1:14){
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

p.pop = matrix(NA,nrow=dim(chi2.obs)[1],ncol=20)
for(i in 1:20){
  tmp.chi2.obs=chi2.obs[,data2$site==i]
  tmp.chi2.sim=chi2.sim[,data2$site==i]
  fit.obs=apply(tmp.chi2.obs,1,sum)
  fit.sim=apply(tmp.chi2.sim,1,sum)
  p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
  p.pop[,i] = p.chi2.calc
}
apply(p.pop,2,mean)


#pdf(file="~/Dropbox/clarkiaSeedBanks/products/figures/modelChecks/seedlingSurvivalFruiting-population.pdf",height=6,width=8)

par(mfrow=c(1,2))
pop.sample = 1:20
plot(y=(pop.sample),x=rev(apply(p.pop,2,mean)),
     xlim=c(0,1),pch=16,
     ylab="Population",xlab="p-Value",
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
p.pop = matrix(NA,nrow=n.iter,ncol=14)
for(j in 1:20){
  for(i in 1:14){
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

time.sample = 1:14+2005
plot(time.sample,p.chi.mat[1,],
     ylim=c(0,1),pch=16,xlim=c(2006,2019),
     xlab="Year",ylab="p-Value",
     main="Fruiting plant counts",type='n')
for(i in 1:20){
  points(time.sample+rnorm(1,0,sd=.05),
         p.chi.mat[i,],pch=19,cex=.5,col='black')
}
abline(h=c(.1,.9),lty='dotted')
abline(h=c(.2,.8),lty='dotted',col='gray')
#dev.off()

cbind(time.sample,apply(ifelse(p.chi.mat>.9,1,0),2,sum)/20)

# per year*population shows some values that are not well modeled this way

pdf(file="~/Dropbox/clarkiaSeedBanks/products/figures/modelChecks/seedlingSurvivalFruiting-zerotrials.pdf",height=6,width=6)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

time.sample = 1:14+2005
for(i in 1:20){
  
  # tot.seeds=c()
  # for(j in 1:14){
  #   tot.seeds[j]=sum(data$seedlingNumber[data$site==i&data$year==j],na.rm=TRUE)
  # }
  # 
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
  ifelse(i%in%c(5), legend(x = 15, y = 1,
                           col = c('gray','orange'),
                           lty = c(1,1),
                           legend = c("Persistence only","Persistence & viability"),
                           cex=.55,
                           box.lty=0), NA)
}
mtext("Year", side = 1, outer = TRUE, line = 2.2)
mtext("Bayesian p-value (chi-2)", side = 2, outer = TRUE, line = 2.2)
#mtext("Population*year-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)
dev.off()

# -------------------------------------------------------------------
# Bayesian p-values: mean
# entire dataset, age-specific
# -------------------------------------------------------------------

y.sim=MCMCchains(mcmcSamples,params=c("fruitplNumber_sim"))

p.test.list = list()
p.test.year = matrix(NA,nrow=n.iter,ncol=20)
for(i in 1:20){
  # here is the test function
  test.data=mean(data$fruitplNumber[data$site==i])
  y_samples=y.sim[,data$site==i]
  test.sim=apply(y_samples,1,mean)
  p.test = ifelse(test.sim-test.data>=0,1,0)
  p.test.year[,i] = p.test
}

dev.off()
par(mfrow=c(1,2))
pop.sample = 1:20
plot(y=(pop.sample),x=rev(apply(p.test.year,2,mean)),
     xlim=c(0,1),pch=16,
     ylab="Population",xlab="p-Value",
     main="Fruiting plant counts",    
     axes=FALSE,frame=FALSE,xaxt='n',yaxt='n')
abline(v=c(.1,.9),lty='dotted')
abline(v=c(.2,.8),lty='dotted',col='gray')

axis(2, (1:20),
     labels = rev(siteNames), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(1,  seq(0,1,by=.2), col.ticks = 1)


n.iter = dim(y.sim)[1]

p.test.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=14)
for(j in 1:20){
  for(i in 1:14){
    tmp.obs=data$seedlingNumber[data$site==j& data$year==i]
    tmp.sim=y.sim[,data$site==j& data$year==i]
    #test.obs=mean(tmp.obs)
    test.sim=apply(tmp.sim,2,mean)
    p.test.calc=ifelse(test.sim-test.obs>=0,1,0)
    p.pop[,i] = p.test.calc
  }
  p.test.list[[j]] = apply(p.pop,2,mean,na.rm=TRUE)
}
p.test.mat=do.call(rbind,p.test.list)

time.sample = 1:14+2005
plot(time.sample,p.test.mat[1,],
     ylim=c(0,1),pch=16,xlim=c(2006,2019),
     xlab="Year",ylab="p-Value",
     main="Fruiting plant counts",type='n')
for(i in 1:20){
  points(time.sample+rnorm(1,0,sd=.05),
         p.test.mat[i,],pch=19,cex=.5,col='black')
}
abline(h=c(.1,.9),lty='dotted')
abline(h=c(.2,.8),lty='dotted',col='gray')

# -------------------------------------------------------------------
# Bayesian p-values: additional
# entire dataset, age-specific
# -------------------------------------------------------------------

y.sim=MCMCchains(mcmcSamples,params=c("fruitplNumber_sim"))

p.sd.list = list()
p.sd.year = matrix(NA,nrow=n.iter,ncol=20)
for(i in 1:20){
    sd.data=sd(data$fruitplNumber[data$site==i])
    y_samples=y.sim[,data$site==i]
    sd.sim=apply(y_samples,1,sd)
    p.sd = ifelse(sd.sim-sd.data>=0,1,0)
    p.sd.year[,i] = p.sd
}

dev.off()
par(mfrow=c(1,2))
pop.sample = 1:20
plot(y=(pop.sample),x=rev(apply(p.sd.year,2,mean)),
     xlim=c(0,1),pch=16,
     ylab="Population",xlab="p-Value",
     main="Fruiting plant counts",    
     axes=FALSE,frame=FALSE,xaxt='n',yaxt='n')
abline(v=c(.1,.9),lty='dotted')
abline(v=c(.2,.8),lty='dotted',col='gray')

axis(2, (1:20),
     labels = rev(siteNames), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(1,  seq(0,1,by=.2), col.ticks = 1)


n.iter = dim(chi2.obs)[1]

p.sd.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=14)
for(j in 1:20){
  for(i in 1:14){
    tmp.obs=data$seedlingNumber[data$site==j& data$year==i]
    tmp.sim=y.sim[,data$site==j& data$year==i]
    sd.obs=sd(tmp.obs)
    sd.sim=apply(tmp.sim,1,sd)
    # fit.obs=apply(tmp.chi2.obs,1,sum)
    # fit.sim=apply(tmp.chi2.sim,1,sum)
    #   exclude=fit.obs!=fit.sim
    #   p.chi2.calc=ifelse(fit.sim[exclude]-fit.obs[exclude]>=0,1,0)
    p.sd.calc=ifelse(sd.sim-sd.obs>=0,1,0)
    p.pop[,i] = p.sd.calc
  }
  p.sd.list[[j]] = apply(p.pop,2,mean,na.rm=TRUE)
}
p.sd.mat=do.call(rbind,p.sd.list)

time.sample = 1:14+2005
plot(time.sample,p.sd.mat[1,],
     ylim=c(0,1),pch=16,xlim=c(2006,2019),
     xlab="Year",ylab="p-Value",
     main="Fruiting plant counts",type='n')
for(i in 1:20){
  points(time.sample+rnorm(1,0,sd=.05),
         p.sd.mat[i,],pch=19,cex=.5,col='black')
}
abline(h=c(.1,.9),lty='dotted')
abline(h=c(.2,.8),lty='dotted',col='gray')

# -------------------------------------------------------------------
# Bayesian p-values: additional mean
# entire dataset, age-specific
# -------------------------------------------------------------------
mean.data=mean(data$seedlingJan[data$siteGermination==1 & data$gIndex==1])
#mean.sim=MCMCchains(mcmcSamples,params=c("mean.sim"))
# plot(apply(MCMCchains(mcmcSamples,params=c("seedlingJan_sim")),1,mean),
#      mean.sim);abline(a=0,b=1)
mean.sim=apply(MCMCchains(mcmcSamples,params=c("seedlingJan_sim"))[,data$siteGermination==1 & data$gIndex==1],1,mean)
p.mean = ifelse(mean.sim-mean.data>=0,1,0)
mean(p.mean)

n.iter=dim(MCMCchains(mcmcSamples,params=c("seedlingJan_sim")))[1]

p.mean.list = list()
p.mean.age = matrix(NA,nrow=n.iter,ncol=3)
for(h in 1:20){
for(i in 1:3){
  mean.data=mean(data$seedlingJan[data$siteGermination==h & data$gIndex==i])
  y_samples=MCMCchains(mcmcSamples,params=c("seedlingJan_sim"))[,data$siteGermination==h & data$gIndex==i]
  mean.sim=apply(y_samples,1,mean)
  p.mean = ifelse(mean.sim-mean.data>=0,1,0)
  p.mean.age[,i] = p.mean
}
  p.mean.list[[h]] = apply(p.mean.age,2,mean)
}
p.mean.mat=do.call(rbind,p.mean.list)


par(mfrow=c(1,1))
time.sample = 1:3
plot(time.sample,p.mean.mat[1,],
     ylim=c(0,1),pch=16,xlim=c(0,4),
     xlab="Months",ylab="p-Value",
     main="Germinant counts",type='n')
for(i in 1:20){
  points(time.sample+rnorm(1,0,sd=.05),
         p.mean.mat[i,],pch=1,cex=.5)
}
abline(h=c(.1,.9),lty='dotted')

# -------------------------------------------------------------------
# Bayesian p-values: additional median
# entire dataset, age-specific
# -------------------------------------------------------------------
median.data=median(data$seedlingJan)
#median.sim=MCMCchains(mcmcSamples,params=c("median.sim"))
# plot(apply(MCMCchains(mcmcSamples,params=c("seedlingJan_sim")),1,median),
#      median.sim);abline(a=0,b=1)
median.sim=apply(MCMCchains(mcmcSamples,params=c("seedlingJan_sim")),1,median)
p.median = ifelse(median.sim-median.data>=0,1,0)
mean(p.median)

p.median.age = matrix(NA,nrow=dim(n.iter)[1],ncol=3)
for(i in 1:3){
  median.data=median(data$seedlingJan[data$gIndex==i])
  y_samples=MCMCchains(mcmcSamples,params=c("seedlingJan_sim"))[,data$gIndex==i]
  median.sim=apply(y_samples,1,median)
  p.median = ifelse(median.sim-median.data>=0,1,0)
  p.median.age[,i] = p.median
}
apply(p.median.age,2,mean)

par(mfrow=c(1,1))
time.sample = 1:3
plot(time.sample,apply(p.median.age,2,mean),
     ylim=c(0,1),pch=16,
     xlab="Months",ylab="p-Value",
     main="Germinant counts")
abline(h=c(.1,.9),lty='dotted')

# -------------------------------------------------------------------
# Bayesian p-values: additional CV
# entire dataset, age-specific
# -------------------------------------------------------------------
cv.data=sd(data$seedlingJan)/mean(data$seedlingJan)
#median.sim=MCMCchains(mcmcSamples,params=c("median.sim"))
# plot(apply(MCMCchains(mcmcSamples,params=c("seedlingJan_sim")),1,median),
#      median.sim);abline(a=0,b=1)
cv.sim=sd.sim/mean.sim
p.cv = ifelse(cv.sim-cv.data>=0,1,0)
mean(p.cv)

p.cv.age = matrix(NA,nrow=dim(n.iter)[1],ncol=3)
for(i in 1:3){
  cv.data=sd(data$seedlingJan[data$gIndex==i])/mean(data$seedlingJan[data$gIndex==i])
  y_samples=MCMCchains(mcmcSamples,params=c("seedlingJan_sim"))[,data$gIndex==i]
  cv.sim = apply(y_samples,1,sd)/apply(y_samples,1,mean)
  p.cv = ifelse(cv.sim-cv.data>=0,1,0)
  p.cv.age[,i] = p.cv
}
apply(p.cv.age,2,mean)

par(mfrow=c(1,1))
time.sample = 1:3
plot(time.sample,apply(p.cv.age,2,mean),
     ylim=c(0,1),pch=16,
     xlab="Months",ylab="p-Value",
     main="Germinant counts")
abline(h=c(.1,.9),lty='dotted')



# -------------------------------------------------------------------
# Bayesian p-values: mean for seed counts
# entire dataset, age-specific
# -------------------------------------------------------------------

par(mfrow=c(1,1))
time.sample = 1:3



p.mean.list = list()
p.mean.age = matrix(NA,nrow=n.iter,ncol=12)
for(h in 1:20){
for(i in 1:12){
  mean.data=mean(data$y[data$siteSurvival==h & data$compIndex==i])
  y_samples=MCMCchains(mcmcSamples,params=c("y_sim"))[,data$siteSurvival==h & data$compIndex==i]
  mean.sim=apply(y_samples,1,mean)
  p.mean = ifelse(mean.sim-mean.data>=0,1,0)
  p.mean.age[,i] = p.mean
}
  p.mean.list[[h]] = apply(p.mean.age,2,mean)
}
p.mean.mat=do.call(rbind,p.mean.list)

par(mfrow=c(1,1))
time.sample = c(3,12,15,24,27,36,
                3,12,15,24,
                3,12)
year.sample = c(1,1,1,1,1,1,
                2,2,2,2,
                3,3)
plot(time.sample,p.mean.mat[1,],
     ylim=c(0,1),pch=16,
     xlab="Months",ylab="p-Value",
     main="Seed counts",type='n')
for(i in 1:20){
  points(time.sample+rnorm(1,0,sd=.05),
         p.mean.mat[i,],pch=1,cex=.5,col=year.sample[i])
}
abline(h=c(.1,.9),lty='dotted')

# -------------------------------------------------------------------
# Bayesian p-values: sd for seed counts
# entire dataset, age-specific
# -------------------------------------------------------------------

p.sd.list = list()
p.sd.age = matrix(NA,nrow=n.iter,ncol=12)
for(h in 1:20){
  for(i in 1:12){
    sd.data=sd(data$y[data$siteSurvival==h & data$compIndex==i])
    y_samples=MCMCchains(mcmcSamples,params=c("y_sim"))[,data$siteSurvival==h & data$compIndex==i]
    sd.sim=apply(y_samples,1,sd)
    p.sd = ifelse(sd.sim-sd.data>=0,1,0)
    p.sd.age[,i] = p.sd
  }
  p.sd.list[[h]] = apply(p.sd.age,2,sd)
}
p.sd.mat=do.call(rbind,p.sd.list)

par(mfrow=c(1,1))
time.sample = c(3,12,15,24,27,36,
                3,12,15,24,
                3,12)
year.sample = c(1,1,1,1,1,1,
                2,2,2,2,
                3,3)
plot(time.sample,p.sd.mat[1,],
     ylim=c(0,1),pch=16,
     xlab="Months",ylab="p-Value",
     main="Seed counts",type='n')
for(i in 1:20){
  points(time.sample+rnorm(1,0,sd=.05),
         p.sd.mat[i,],pch=1,cex=.5,col=year.sample[i])
}
abline(h=c(.1,.9),lty='dotted')

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


