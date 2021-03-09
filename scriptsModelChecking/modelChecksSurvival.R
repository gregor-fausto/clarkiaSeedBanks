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

# -------------------------------------------------------------------
# Bayesian p-values: omnibus Chi-squared for counts of fruiting plants
# entire dataset, population-level
# -------------------------------------------------------------------
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

par(mfrow=c(1,1))
pop.sample = 1:20
plot(pop.sample,apply(p.pop,2,mean),
     ylim=c(0,1),pch=16,
     xlab="Population",ylab="p-Value",
     main="Fruiting plant counts",    
     axes=FALSE,frame=FALSE,xaxt='n',yaxt='n')
abline(h=c(.1,.9),lty='dotted')

axis(1, (1:20),
     labels = siteNames, las = 1, 
     col = NA, col.ticks = 1, cex.axis = .5)
axis(2,  seq(0,1,by=.2), col.ticks = 1)

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
    p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
    p.pop[,i] = p.chi2.calc
  }
  p.chi.list[[j]] = apply(p.pop,2,mean)
}
p.chi.mat=do.call(rbind,p.chi.list)


par(mfrow=c(1,1))
time.sample = 1:14+2005
plot(time.sample,p.chi.mat[1,],
     ylim=c(0,1),pch=16,xlim=c(2006,2019),
     xlab="Months",ylab="p-Value",
     main="Germinant counts",type='n')
for(i in 1:20){
  points(time.sample+rnorm(1,0,sd=.05),
         p.chi.mat[i,],pch=1,cex=.5)
}
abline(h=c(.1,.9),lty='dotted')

# per year*population shows some values that are not well modeled this way

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

time.sample = 1:14+2005
for(i in 1:20){
  plot(NA,NA,
       ylim=c(0,1),pch=16,xlim=c(2006,2019),
       ylab='',xlab='',xaxt='n',yaxt='n')
  points(time.sample+rnorm(1,0,sd=.05),
         p.chi.mat[i,],pch=19,cex=1,
         col=ifelse(p.chi.mat[i,]>.9|p.chi.mat[i,]<.1,"purple","black"))
  
  abline(h=c(.1,.9),lty='dotted')
  
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
mtext("Population*year-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)


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
       ylim=c(0,1),pch=16,xlim=c(0,max(tot.seeds)),
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
mtext("Population*year-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)

# years with low number of seedlings (trials) and fruiting plants (successes)
# pull the mean towards the population-level average
# due to pooling

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

# -------------------------------------------------------------------
# Bayesian p-values: additional
# entire dataset, age-specific
# -------------------------------------------------------------------

sd.data=sd(data$fruitplNumber)
sd.sim=apply(MCMCchains(mcmcSamples,params=c("fruitplNumber_sim")),1,sd)
p.sd = ifelse(sd.sim-sd.data>=0,1,0)
mean(p.sd)

p.sd.list = list()
p.sd.year = matrix(NA,nrow=n.iter,ncol=20)
for(h in 1:20){
    sd.data=sd(data$fruitplNumber[data$site==h])
    y_samples=MCMCchains(mcmcSamples,params=c("fruitplNumber_sim"))[,data$site==h]
    sd.sim=apply(y_samples,1,sd)
    p.sd = ifelse(sd.sim-sd.data>=0,1,0)
    p.sd.year[,i] = p.sd
}
p.sd.mat=do.call(rbind,p.sd.list)


par(mfrow=c(1,1))
pop.sample = 1:20
plot(pop.sample,p.sd.mat[1,],
     ylim=c(0,1),pch=16,xlim=c(0,4),
     xlab="Months",ylab="p-Value",
     main="Germinant counts",type='n')
for(i in 1:20){
  points(pop.sample+rnorm(1,0,sd=.05),
         p.sd.mat[i,],pch=1,cex=.5)
}
abline(h=c(.1,.9),lty='dotted')

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


