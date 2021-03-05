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

mcmcSamples <- readRDS(modelFittingFiles[[2]])
data <- readRDS(modelFittingFiles[[1]])


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Convergence diagnostics
# -------------------------------------------------------------------
# -------------------------------------------------------------------
MCMCsummary(mcmcSamples, params = c("mu0_g"))
MCMCsummary(mcmcSamples, params = c("sigma0_g"))
MCMCsummary(mcmcSamples, params = c("sigma_g"))

MCMCsummary(mcmcSamples, params = c("mu0_s"))
MCMCsummary(mcmcSamples, params = c("sigma0_s"))
MCMCsummary(mcmcSamples, params = c("sigma_s"))

MCMCsummary(mcmcSamples, params = c("a"))

alpha<-MCMCchains(mcmcSamples, params = c("a"))
alpha.sum<-apply(alpha,2,quantile,c(.025,.5,.975))

par(mfrow=c(1,1))

plot(NA,NA,type='n',xlim=c(0,2),ylim=c(0,20))
for(i in 1:20){
  tmp<-alpha.sum[,i]
  segments(x0=tmp[1],x1=tmp[3],y0=i)
  points(x=tmp[2],y=i,pch=19)
}

for(i in 1:20){
  hist(MCMCchains(mcmcSamples,params="a")[,i],
       breaks=100,freq=FALSE,xlim=c(0,2));
  abline(v=1,col='red')
}

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Graphical Posterior Predictive Checks
# -------------------------------------------------------------------
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# Germination
# -------------------------------------------------------------------
par(mfrow=c(1,3))

for(i in 1:3){
  hist(data$seedlingJan[data$siteGermination==1&data$gIndex==i], breaks = 10, 
       freq = FALSE, main = "Simulated and real data for germination", 
       xlab = expression(paste("germinant count")), cex.lab = 1.2) 
  seedlingJan_sim=MCMCchains(mcmcSamples, params = "seedlingJan_sim")[,data$siteGermination==1&data$gIndex==i]
  lines(density(seedlingJan_sim,adjust=5), col = "red")
}


# -------------------------------------------------------------------
# Intact seed counts
# -------------------------------------------------------------------
par(mfrow=c(2,3))

for(i in 1:6){
  hist(data$y[data$siteSurvival==1&data$compIndex==i], breaks = 10, 
       freq = FALSE, main = "Simulated and real data for germination", 
       xlab = expression(paste("germinant count")), cex.lab = 1.2) 
  y_sim=MCMCchains(mcmcSamples, params = "y_sim")[,data$siteSurvival==1&data$compIndex==i]
  lines(density(y_sim,adjust=5), col = "red")
}

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

# -------------------------------------------------------------------
# Bayesian p-values: omnibus Chi-squared for germination counts
# entire dataset, age-specific
# -------------------------------------------------------------------
chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.obs"))
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.sim"))
# calculations are rowwise
fit.obs=apply(chi2.obs,1,sum)
fit.sim=apply(chi2.sim,1,sum)
p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
mean(p.chi2.calc)

#MCMCsummary(mcmcSamples, params = c("p.chi2"))
# hist(MCMCchains(mcmcSamples, params = c("p.chi2")))
# hist(p.chi2.calc)

p.age = matrix(NA,nrow=dim(chi2.obs)[1],ncol=3)
for(i in 1:3){
  chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.obs"))[,data$gIndex==i]
  chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.sim"))[,data$gIndex==i]
  fit.obs=apply(chi2.obs,1,sum)
  fit.sim=apply(chi2.sim,1,sum)
  p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
  p.age[,i] = p.chi2.calc
}
apply(p.age,2,mean)

par(mfrow=c(1,1))
time.sample = 1:3
plot(time.sample,apply(p.age,2,mean),
     ylim=c(0,1),pch=16,
     xlab="Months",ylab="p-Value",
     main="Germinant counts")
abline(h=c(.1,.9),lty='dotted')

# note that for particular parameter combinations the counts
# are very low in year 3, making it challenging to obtain reasonable
# estimates of conditional germination rates
# see this by reducing lambda and comparing
# but I think this might trade off with ability to estimate survival

# -------------------------------------------------------------------
# Bayesian p-values: omnibus Chi-squared for seed counts
# entire dataset, age-specific
# -------------------------------------------------------------------
chi2.yobs=MCMCchains(mcmcSamples,params=c("chi2.yobs"))
chi2.ysim=MCMCchains(mcmcSamples,params=c("chi2.ysim"))
# calculations are rowwise
fit.obs=apply(chi2.yobs,1,sum)
fit.sim=apply(chi2.ysim,1,sum)
p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
mean(p.chi2.calc)

# hist(MCMCchains(mcmcSamples, params = c("p.chi2")))
# hist(p.chi2.calc)

p.index = matrix(NA,nrow=dim(chi2.obs)[1],ncol=6)
for(i in 1:6){
  chi2.yobs.index=chi2.yobs[,data$compIndex==i]
  chi2.ysim.index=chi2.ysim[,data$compIndex==i]
  fit.obs=apply(chi2.yobs.index,1,sum)
  fit.sim=apply(chi2.ysim.index,1,sum)
  p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
  p.index[,i] = p.chi2.calc
}
apply(p.index,2,mean)

par(mfrow=c(1,1))
time.sample = c(3,12,15,24,27,36)
plot(time.sample,apply(p.index,2,mean),
     ylim=c(0,1),pch=16,
     xlab="Months",ylab="p-Value",
     main="Seed counts")
abline(h=c(.1,.9),lty='dotted')

# -------------------------------------------------------------------
# Bayesian p-values: additional
# entire dataset, age-specific
# -------------------------------------------------------------------


sd.data=sd(data[data$siteGermination==1]$seedlingJan)
#sd.sim=MCMCchains(mcmcSamples,params=c("sd.sim"))
# plot(apply(MCMCchains(mcmcSamples,params=c("seedlingJan_sim")),1,sd),
#      sd.sim);abline(a=0,b=1)
sd.sim=apply(MCMCchains(mcmcSamples,params=c("seedlingJan_sim"))[,data$siteGermination==1],1,sd)
p.sd = ifelse(sd.sim-sd.data>=0,1,0)
mean(p.sd)

n.iter=dim(MCMCchains(mcmcSamples,params=c("seedlingJan_sim")))[1]

p.sd.list = list()
p.sd.age = matrix(NA,nrow=n.iter,ncol=3)
for(h in 1:20){
  for(i in 1:3){
    sd.data=sd(data$seedlingJan[data$siteGermination==h & data$gIndex==i])
    y_samples=MCMCchains(mcmcSamples,params=c("seedlingJan_sim"))[,data$siteGermination==h & data$gIndex==i]
    sd.sim=apply(y_samples,1,sd)
    p.sd = ifelse(sd.sim-sd.data>=0,1,0)
    p.sd.age[,i] = p.sd
  }
  p.sd.list[[h]] = apply(p.sd.age,2,sd)
}
p.sd.mat=do.call(rbind,p.sd.list)


par(mfrow=c(1,1))
time.sample = 1:3
plot(time.sample,p.sd.mat[1,],
     ylim=c(0,1),pch=16,xlim=c(0,4),
     xlab="Months",ylab="p-Value",
     main="Germinant counts",type='n')
for(i in 1:20){
  points(time.sample+rnorm(1,0,sd=.05),
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


