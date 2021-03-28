# rm(list=ls(all=TRUE)) # clear R environment
rm(list=setdiff(ls(all=TRUE),c("fileDirectory","outputDirectory"))) # if using in source(script)
options(stringsAsFactors = FALSE)

library(MCMCvis)
library(tidybayes)
library(tidyverse)
library(magrittr)
library(bayesplot)
library(rethinking)

mcmcSampleDirectory <- paste0(fileDirectory,list.files(fileDirectory))
mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("seedSamples.rds",mcmcSampleDirectory)]])
data <- readRDS(mcmcSampleDirectory[[grep("seedData.rds",mcmcSampleDirectory)]])

# -------------------------------------------------------------------
# Get site names and position
# -------------------------------------------------------------------

position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,dominant.surface.rock.type) %>%
  dplyr::mutate(easting=easting/1000)

siteNames <- unique(position$site)

# -------------------------------------------------------------------
# Incorporating loss of viability 
# -------------------------------------------------------------------

mcmcSamplesViability <- readRDS(mcmcSampleDirectory[[grep("viabilityTrialSamples.rds",mcmcSampleDirectory)]])
dataViability <- readRDS(mcmcSampleDirectory[[grep("viabilityData.rds",mcmcSampleDirectory)]])

# -------------------------------------------------------------------
# Germination posterior
# -------------------------------------------------------------------

mu0_g=MCMCchains(mcmcSamples,params="mu0_g")

index.g1=grep(",1\\]",colnames(mu0_g))
index.g2=grep(",2\\]",colnames(mu0_g))
index.g3=grep(",3\\]",colnames(mu0_g))

gamma1 = boot::inv.logit(mu0_g[,index.g1])
gamma1.sum=apply(gamma1,2,quantile,probs=c(0.025,.25,.5,.75,.975))

gamma2 = boot::inv.logit(mu0_g[,index.g2])
gamma2.sum=apply(gamma2,2,quantile,probs=c(0.025,.25,.5,.75,.975))

gamma3 = boot::inv.logit(mu0_g[,index.g3])
gamma3.sum=apply(gamma3,2,quantile,probs=c(0.025,.25,.5,.75,.975))

# -------------------------------------------------------------------
# Discrete survival component
# -------------------------------------------------------------------

theta_c2 = (1-gamma1)
theta_c4 = (1-gamma2)
theta_c6 = (1-gamma3)

# -------------------------------------------------------------------
# Discretize product integral of survival function
# -------------------------------------------------------------------
a.wb=MCMCchains(mcmcSamples,params="a")
b0.wb=MCMCchains(mcmcSamples,params="mu0_s")

inv.b0.wb=(exp(-b0.wb/a.wb))

x = c(1,2,2,3,4,4,5,6,6,7)
t = c(0,4,12,16,24,28,36)/36

theta_0 = 1
theta_1 = function(t,inv.b0,alpha){
  exp(-(t/inv.b0)^alpha)
}
th_1=theta_1(t[x[2]],inv.b0=inv.b0.wb,alpha=a.wb)

theta_2 = function(t,inv.b0,alpha,theta_c2){
  theta_c2*exp(-(t/inv.b0)^alpha)
}
th_2=theta_2(t[x[3]],inv.b0=inv.b0.wb,alpha=a.wb,theta_c2=theta_c2)

theta_3 = theta_2
th_3=theta_3(t[x[4]],inv.b0=inv.b0.wb,alpha=a.wb,theta_c2=theta_c2)

theta_4 = theta_2
th_4=theta_4(t[x[5]],inv.b0=inv.b0.wb,alpha=a.wb,theta_c2=theta_c2)

theta_5 = function(t,inv.b0,alpha,theta_c2=theta_c2,theta_c4=theta_c4){
  theta_c4*theta_c2*exp(-(t/inv.b0)^alpha)
}
th_5=theta_5(t[x[6]],inv.b0=inv.b0.wb,alpha=a.wb,theta_c2=theta_c2,theta_c4=theta_c4)

theta_6 = theta_5
th_6=theta_6(t[x[7]],inv.b0=inv.b0.wb,alpha=a.wb,theta_c2=theta_c2,theta_c4=theta_c4)

theta_7 = theta_5
th_7=theta_7(t[x[8]],inv.b0=inv.b0.wb,alpha=a.wb,theta_c2=theta_c2,theta_c4=theta_c4)

theta_8 = function(t,inv.b0,alpha,theta_c2=theta_c2,theta_c4=theta_c4,theta_c6=theta_c6){
  theta_c6*theta_c4*theta_c2*exp(-(t/inv.b0)^alpha)
}
th_8=theta_8(t[x[9]],inv.b0=inv.b0.wb,alpha=a.wb,theta_c2=theta_c2,theta_c4=theta_c4,theta_c6=theta_c6)

theta_9 = function(t,inv.b0,alpha,theta_c2=theta_c2,theta_c4=theta_c4,theta_c6=theta_c6){
  theta_c6*theta_c4*theta_c2*exp(-(t/inv.b0)^alpha)
}
th_9=theta_9(t[x[10]],inv.b0=inv.b0.wb,alpha=a.wb,theta_c2=theta_c2,theta_c4=theta_c4,theta_c6=theta_c6)

# fill in theta_0
th_0=matrix(1,nrow=dim(th_9)[1],ncol=dim(th_9)[2])

# -------------------------------------------------------------------
# Viability
# -------------------------------------------------------------------

# get population estimates for year 1-3
mu0_g<-MCMCchains(mcmcSamplesViability, params = "mu0_g")
mu0_v<-MCMCchains(mcmcSamplesViability, params = "mu0_v")
# calculate total viability; calculation on latent scale
nu0=exp(mu0_g)/(1+exp(mu0_g))+(exp(-mu0_g+mu0_v)/(1+exp(-mu0_g)+exp(mu0_v)+exp(-mu0_g+mu0_v)))

# interpolate for January estimates
nu0_ratio1=(nu0[,index.g1])^(1/3)
nu0_ratio2=nu0[,index.g1]*(nu0[,index.g2]/nu0[,index.g1])^(1/3)
nu0_ratio3=nu0[,index.g2]*(nu0[,index.g3]/nu0[,index.g2])^(1/3)
nu0_ratio2=nu0_ratio2
nu0_ratio3=nu0_ratio3

nu_01.inter = nu0_ratio1
nu_1 = nu0[,index.g1]
nu_12.inter = nu0_ratio2
nu_2 = nu0[,index.g2]
nu_23.inter = nu0_ratio3
nu_3 = nu0[,index.g3]

# -------------------------------------------------------------------
# Combine survival function + viability
# tried 2 ways: first as survival function at each time step
# then by combining age-specific survival function
# the second is what I use
# -------------------------------------------------------------------

## conditional germination 

gamma1.v = gamma1/(1-(1-nu_01.inter)*(1-gamma1))
gamma2.v = gamma2/(1-(1-nu_12.inter)*(1-gamma2))
gamma3.v = gamma3/(1-(1-nu_23.inter)*(1-gamma3))

f.wb = function(t,inv.b0,alpha){
  exp(-(t/inv.b0)^alpha)
}

## age-specific probability formulation function formulation

# Calculate age-specific survival probability
Oct_0.v = f.wb(t=t[x[1]],inv.b0=inv.b0.wb,alpha=a.wb)

Janpre_1.v = f.wb(t=t[x[2]],inv.b0=inv.b0.wb,alpha=a.wb)*(gamma1+(1-gamma1)*nu_01.inter)
Jangerm_1.v = gamma1.v
Janpost_1.v = (1-gamma1.v)
Oct_1.v = (f.wb(t=t[x[4]],inv.b0=inv.b0.wb,alpha=a.wb)*nu_1)/(f.wb(t=t[x[3]],inv.b0=inv.b0.wb,alpha=a.wb)*nu_01.inter)

Janpre_2.v =  (f.wb(t=t[x[5]],inv.b0=inv.b0.wb,alpha=a.wb)*(gamma2+(1-gamma2)*nu_12.inter))/(f.wb(t=t[x[4]],inv.b0=inv.b0.wb,alpha=a.wb)*nu_1)
Jangerm_2.v = gamma2.v
Janpost_2.v = (1-gamma2.v)
Oct_2.v = (f.wb(t=t[x[7]],inv.b0=inv.b0.wb,alpha=a.wb)*nu_2)/(f.wb(t=t[x[6]],inv.b0=inv.b0.wb,alpha=a.wb)*nu_12.inter)

Janpre_3.v = (f.wb(t=t[x[8]],inv.b0=inv.b0.wb,alpha=a.wb)*(gamma3+(1-gamma3)*nu_23.inter))/(f.wb(t=t[x[7]],inv.b0=inv.b0.wb,alpha=a.wb)*nu_2)
Jangerm_3.v = gamma3.v
Janpost_3.v = (1-gamma3.v)
Oct_3.v = (f.wb(t=t[x[10]],inv.b0=inv.b0.wb,alpha=a.wb)*nu_3)/(f.wb(t=t[x[9]],inv.b0=inv.b0.wb,alpha=a.wb)*nu_23.inter)

# Calculate survivorship schedule
phi_0=Oct_0.v
phi_1=Janpre_1.v
phi_2=Janpost_1.v*Janpre_1.v
phi_3=Oct_1.v*Janpost_1.v*Janpre_1.v
phi_4=Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v
phi_5=Janpost_2.v*Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v
phi_6=Oct_2.v*Janpost_2.v*Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v
phi_7=Janpre_3.v*Oct_2.v*Janpost_2.v*Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v
phi_8=Janpost_3.v*Janpre_3.v*Oct_2.v*Janpost_2.v*Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v
phi_9=Oct_3.v*Janpost_3.v*Janpre_3.v*Oct_2.v*Janpost_2.v*Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v

# -------------------------------------------------------------------
# Check for properties of the survival function
# that no longer make it a proper survival function;
# specifically the pmf of the survival function has values <0
# commented out nwo but the code below shows how this is the case for the initial
# survival function
# -------------------------------------------------------------------
## if calculating survival by incorporating viability, end up with survival > 1
## by including viability, the trajectory is no longer a proper survival function
# plot(rep(1,20),apply(phi_0-phi_1,2,min),col=ifelse(apply(phi_0-phi_1,2,min) <0,'orange','black'),ylim=c(-.25,1.25),xlim=c(.5,9.5))
# points(rep(1,20),apply(phi_0-phi_1,2,max),col=ifelse(apply(phi_0-phi_1,2,max) >1,'orange','black'),ylim=c(-.25,1.25))
# points(rep(2,20),apply(phi_1-phi_2,2,min),col=ifelse(apply(phi_1-phi_2,2,min) <0,'orange','black'))
# points(rep(2,20),apply(phi_1-phi_2,2,max),col=ifelse(apply(phi_1-phi_2,2,max) >1,'orange','black'))
# del = phi_2-phi_3
# points(rep(3,20),apply(del,2,min),col=ifelse(apply(del,2,min) <0,'orange','black'))
# points(rep(3,20),apply(del,2,max),col=ifelse(apply(del,2,max) >1,'orange','black'))
# del = phi_3-phi_4
# points(rep(4,20),apply(del,2,min),col=ifelse(apply(del,2,min) <0,'orange','black'))
# points(rep(4,20),apply(del,2,max),col=ifelse(apply(del,2,max) >1,'orange','black'))
# del = phi_4-phi_5
# points(rep(5,20),apply(del,2,min),col=ifelse(apply(del,2,min) <0,'orange','black'))
# points(rep(5,20),apply(del,2,max),col=ifelse(apply(del,2,max) >1,'orange','black'))
# del = phi_5-phi_6
# points(rep(6,20),apply(del,2,min),col=ifelse(apply(del,2,min) <0,'orange','black'))
# points(rep(6,20),apply(del,2,max),col=ifelse(apply(del,2,max) >1,'orange','black'))
# del = phi_6-phi_7
# points(rep(7,20),apply(del,2,min),col=ifelse(apply(del,2,min) <0,'orange','black'))
# points(rep(7,20),apply(del,2,max),col=ifelse(apply(del,2,max) >1,'orange','black'))
# del = phi_7-phi_8
# points(rep(8,20),apply(del,2,min),col=ifelse(apply(del,2,min) <0,'orange','black'))
# points(rep(8,20),apply(del,2,max),col=ifelse(apply(del,2,max) >1,'orange','black'))
# del = phi_8-phi_9
# points(rep(9,20),apply(del,2,min),col=ifelse(apply(del,2,min) <0,'orange','black'))
# points(rep(9,20),apply(del,2,max),col=ifelse(apply(del,2,max) >1,'orange','black'))

## check distribution of difference of curves
# par(mfrow = c(4,5),
#     oma = c(5,4,0,0) + 0.1,
#     mar = c(0,0,1,1) + 0.1)
# # for(i in 1:20){hist(phi_0[,i]-phi_1[,i],main='');abline(v=0,col='red')}
# # for(i in 1:20){hist(phi_1[,i]-phi_2[,i],main='');abline(v=0,col='red')}
# # for(i in 1:20){hist(phi_2[,i]-phi_3[,i],main='');abline(v=c(0,1),col='red')}
# for(i in 1:20){hist(phi_3[,i]-phi_4[,i],main='');abline(v=c(0,1),col='red')}
# # for(i in 1:20){hist(phi_4[,i]-phi_5[,i],main='');abline(v=0,col='red')}
# for(i in 1:20){hist(phi_5[,i]-phi_6[,i],main='');abline(v=0,col='red')}
# for(i in 1:20){hist(phi_6[,i]-phi_7[,i],main='');abline(v=0,col='red')}
# # for(i in 1:20){hist(phi_7[,i]-phi_8[,i],main='');abline(v=0,col='red')}
# for(i in 1:20){hist(phi_8[,i]-phi_9[,i],main='');abline(v=0,col='red')}

## check distribution of components
# # for(i in 1:20){hist(Oct_0.v[,i],main='');abline(v=c(0,1),col='red')}
# # for(i in 1:20){hist(Janpre_1.v[,i],main='');abline(v=c(0,1),col='red')}
# # for(i in 1:20){hist(Janpost_1.v[,i],main='');abline(v=c(0,1),col='red')}
# #for(i in 1:20){hist(Oct_1.v[,i],main='');abline(v=c(0,1),col='red')}
# for(i in 1:20){hist(Janpre_2.v[,i],main='');abline(v=c(0,1),col='red')}
# #for(i in 1:20){hist(Janpost_2.v[,i],main='');abline(v=c(0,1),col='red')}
# for(i in 1:20){hist(Oct_2.v[,i],main='');abline(v=c(0,1),col='red')}
# for(i in 1:20){hist(Janpre_3.v[,i],main='');abline(v=c(0,1),col='red')}
# #for(i in 1:20){hist(Janpost_3.v[,i],main='');abline(v=c(0,1),col='red')}
# for(i in 1:20){hist(Oct_3.v[,i],main='');abline(v=c(0,1),col='red')}

## examine Janpre_2.v
# par(mfrow = c(4,5),
#     oma = c(5,4,0,0) + 0.1,
#     mar = c(0,0,1,1) + 0.1)
# test = (f.wb(t=t[x[5]],inv.b0=inv.b0.wb,alpha=a.wb)*(gamma2+(1-gamma2)*nu_12.inter))/(f.wb(t=t[x[4]],inv.b0=inv.b0.wb,alpha=a.wb)*nu_1)
# for(i in 1:20){hist(test[,i],main='');abline(v=c(0,1),col='red')}
# test = (f.wb(t=t[x[5]],inv.b0=inv.b0.wb,alpha=a.wb)*(gamma2+(1-gamma2)*nu_12.inter))/(f.wb(t=t[x[4]],inv.b0=inv.b0.wb,alpha=a.wb)*nu_1)
# for(i in 1:20){hist(test[,i],main='');abline(v=c(0,1),col='red')}
# 
# numerator=(f.wb(t=t[x[5]],inv.b0=inv.b0.wb,alpha=a.wb)*(gamma2+(1-gamma2)*nu_12.inter))
# denominator = (f.wb(t=t[x[4]],inv.b0=inv.b0.wb,alpha=a.wb)*nu_1)
# plot(numerator[,18],denominator[,18]);abline(a=0,b=1,lwd=2,col='red')
# 
# hist((numerator/denominator)[,18],breaks=100)
# hist(ifelse((numerator/denominator)[,18]>1,1,(numerator/denominator)[,18]),add=TRUE,col='red',breaks=100)
# 
# numerator[,18]=ifelse(numerator[,18]<denominator[,18],numerator[,18],denominator[,18])
# hist((numerator/denominator)[,18],add=TRUE,col='blue',breaks=100)
# 
# plot(numerator[,18],denominator[,18]);abline(a=0,b=1,lwd=2,col='red')
# 
# numerator18=numerator[,18];denominator18=denominator[,18]
# plot(numerator18[order(numerator18)],denominator18[order(denominator18)]);abline(a=0,b=1,lwd=2,col='red')
# hist(denominator18-numerator18,breaks=100)
# hist(denominator18[order(denominator18)]-numerator18[order(numerator18)],breaks=100)

# end message:
# if calculating survival by incorporating viability, end up with survival > 1
# by including viability, the trajectory is no longer a proper survival function

# -------------------------------------------------------------------
# Truncate distributions and recalculate survival function
# here I considered 3 options
# 1. truncate distribution: this piled a lot of probability mass up at 1
# 2. truncate distribution and resample from same distribution [chose this option]
# 3. resample but from the distribution of the transition without viability
# -------------------------------------------------------------------
# truncate distributions
# hist(ifelse(x>1,resample(x),x),col='red',breaks=100,add=TRUE)
# Janpre_2.v=ifelse(Janpre_2.v>1,1,Janpre_2.v)
# Oct_2.v=ifelse(Oct_2.v>1,1,Oct_2.v)
# Janpre_3.v=ifelse(Janpre_3.v>1,1,Janpre_3.v)
# Oct_3.v=ifelse(Oct_3.v>1,1,Oct_3.v)

# instead of truncating try redistributing by sampling from same distribution?
# resample with replacement

# x = Janpre_2.v[,1]
resample=function(x){
  tmp = ifelse(x<=1,x,NA)
  len=sum(is.na(tmp))
  tmp2=x[x<=1]
  out = ifelse(is.na(tmp),sample(tmp2,len,replace=TRUE),tmp)
  return(out)
}

Janpre_2.v=apply(Janpre_2.v,2,resample)
Oct_2.v=apply(Oct_2.v,2,resample)
Janpre_3.v=apply(Janpre_3.v,2,resample)
Oct_3.v=apply(Oct_3.v,2,resample)

# -------------------------------------------------------------------
# Diagnostic
# -------------------------------------------------------------------
## Commenting this out but this was a useful diagnostic to complement above
## check if there are any PMF values <0
# par(mfrow = c(4,5),
#     oma = c(5,4,0,0) + 0.1,
#     mar = c(0,0,1,1) + 0.1)
# for(i in 1:20){
#   plot(NA,NA,xlim=c(0,10),ylim=c(-1,1),ylab='',xlab='',xaxt='n',yaxt='n')
#   abline(h=0)
#   surv.fun_persistence=apply(cbind(th_0[,i],th_1[,i],th_2[,i],th_3[,i],th_4[,i],th_5[,i],th_6[,i],th_7[,i],th_8[,i],th_9[,i]),2,median)
#   # lines(t[x]*36,surv.fun_persistence,lty='solid',col='lightgray')
#   # points(t[x]*36,surv.fun_persistence,pch=16,col='black')
#   
#   surv.fun_persistence.viability=apply(cbind(phi_0[,i],phi_1[,i],phi_2[,i],phi_3[,i],phi_4[,i],phi_5[,i],phi_6[,i],phi_7[,i],phi_8[,i],phi_9[,i]),2,median)
#   # lines(t[x]*36,surv.fun_persistence.viability,lty='solid',col=
#   #         ifelse(all(diff(surv.fun_persistence.viability) <= 0),'orange','blue'))
#   # points(t[x]*36,surv.fun_persistence.viability,pch=16,
#   #        col=ifelse(all(diff(surv.fun_persistence.viability) <= 0),'orange','blue'))
#   
#   pmf = surv.fun_persistence.viability[1:9]-surv.fun_persistence.viability[2:10]
#   h = pmf/surv.fun_persistence.viability[1:9]
#   points(1:9,h,col=ifelse(h<0,"orange","blue"))
#   
#   text(x=3,y=-.9,siteNames[i])
#   ifelse(i%in%c(16:20),axis(1L),NA)
#   ifelse(i%in%c(1,6,11,16),axis(2L),NA)
# }


# -------------------------------------------------------------------
# Recalculate survivorship schedule using truncated transitions
# -------------------------------------------------------------------
phi_0=Oct_0.v
phi_1=Janpre_1.v
phi_2=Janpost_1.v*Janpre_1.v
phi_3=Oct_1.v*Janpost_1.v*Janpre_1.v
phi_4=Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v
phi_5=Janpost_2.v*Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v
phi_6=Oct_2.v*Janpost_2.v*Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v
phi_7=Janpre_3.v*Oct_2.v*Janpost_2.v*Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v
phi_8=Janpost_3.v*Janpre_3.v*Oct_2.v*Janpost_2.v*Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v
phi_9=Oct_3.v*Janpost_3.v*Janpre_3.v*Oct_2.v*Janpost_2.v*Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v


## Check to make sure no values <0
# if calculating survival by incorporating viability, end up with survival > 1
# by including viability, the trajectory is no longer a proper survival function
# code written so that values <0 show up in orange; should all be black to pass check
dev.off()
plot(rep(1,20),apply(phi_0-phi_1,2,min),
     col=ifelse(apply(phi_0-phi_1,2,min) <0,'orange','black'),
     ylim=c(-.25,1.25),xlim=c(.5,9.5),
     xlab="Probability Mass Function",ylab="Min (lower) or Max (upper)")
points(rep(1,20),apply(phi_0-phi_1,2,max),col=ifelse(apply(phi_0-phi_1,2,max) >1,'orange','black'),ylim=c(-.25,1.25))
points(rep(2,20),apply(phi_1-phi_2,2,min),col=ifelse(apply(phi_1-phi_2,2,min) <0,'orange','black'))
points(rep(2,20),apply(phi_1-phi_2,2,max),col=ifelse(apply(phi_1-phi_2,2,max) >1,'orange','black'))
del = phi_2-phi_3
points(rep(3,20),apply(del,2,min),col=ifelse(apply(del,2,min) <0,'orange','black'))
points(rep(3,20),apply(del,2,max),col=ifelse(apply(del,2,max) >1,'orange','black'))
del = phi_3-phi_4
points(rep(4,20),apply(del,2,min),col=ifelse(apply(del,2,min) <0,'orange','black'))
points(rep(4,20),apply(del,2,max),col=ifelse(apply(del,2,max) >1,'orange','black'))
del = phi_4-phi_5
points(rep(5,20),apply(del,2,min),col=ifelse(apply(del,2,min) <0,'orange','black'))
points(rep(5,20),apply(del,2,max),col=ifelse(apply(del,2,max) >1,'orange','black'))
del = phi_5-phi_6
points(rep(6,20),apply(del,2,min),col=ifelse(apply(del,2,min) <0,'orange','black'))
points(rep(6,20),apply(del,2,max),col=ifelse(apply(del,2,max) >1,'orange','black'))
del = phi_6-phi_7
points(rep(7,20),apply(del,2,min),col=ifelse(apply(del,2,min) <0,'orange','black'))
points(rep(7,20),apply(del,2,max),col=ifelse(apply(del,2,max) >1,'orange','black'))
del = phi_7-phi_8
points(rep(8,20),apply(del,2,min),col=ifelse(apply(del,2,min) <0,'orange','black'))
points(rep(8,20),apply(del,2,max),col=ifelse(apply(del,2,max) >1,'orange','black'))
del = phi_8-phi_9
points(rep(9,20),apply(del,2,min),col=ifelse(apply(del,2,min) <0,'orange','black'))
points(rep(9,20),apply(del,2,max),col=ifelse(apply(del,2,max) >1,'orange','black'))


# -------------------------------------------------------------------
# Plot survival function: persistence only and persistence+viability
# -------------------------------------------------------------------
pdf("~/Dropbox/clarkiaSeedBanks/products/figures/survival-function-persistence-viability.pdf",width=8,height=6)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  plot(NA,NA,xlim=c(0,36),ylim=c(0,1),ylab='',xlab='',xaxt='n',yaxt='n')
  
  surv.fun_persistence=apply(cbind(th_0[,i],th_1[,i],th_2[,i],th_3[,i],th_4[,i],th_5[,i],th_6[,i],th_7[,i],th_8[,i],th_9[,i]),2,median)
  lines((t[x]*36)[1:10],surv.fun_persistence[1:10],lty='solid',col='black')
  #points((t[x]*36)[1:10],surv.fun_persistence[1:10],pch=16,col='black')
  
  surv.fun_persistence.viability=apply(cbind(phi_0[,i],phi_1[,i],phi_2[,i],phi_3[,i],phi_4[,i],phi_5[,i],phi_6[,i],phi_7[,i],phi_8[,i],phi_9[,i]),2,median)
  lines((t[x]*36)[1:10],surv.fun_persistence.viability[1:10],lty='solid',col=
          ifelse(all(diff(surv.fun_persistence.viability[1:10]) <= 0),'orange','blue'))
  # points((t[x]*36)[1:10],surv.fun_persistence.viability[1:10],pch=16,
  #        col=ifelse(all(diff(surv.fun_persistence.viability[1:10]) <= 0),'orange','blue'))
  text(x=-2,y=.05,siteNames[i],pos=4)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
  ifelse(i%in%c(5), legend(x = 15, y = 1,
                           col = c('gray','orange'),
                           lty = c(1,1),
                           legend = c("Persistence only","Persistence & viability"),
                           cex=.55,
                           box.lty=0), NA)
}
mtext("Time (months)", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of persistence and viability", side = 2, outer = TRUE, line = 2.2)
mtext("Population-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)

# -------------------------------------------------------------------
# Plot discrete KM-style survival function plot
# -------------------------------------------------------------------

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  plot(NA,NA,xlim=c(0,36),ylim=c(0,1),ylab='',xlab='',xaxt='n',yaxt='n')
  surv.fun=apply(cbind(th_0[,i],th_1[,i],th_2[,i],th_3[,i],th_4[,i],th_5[,i],th_6[,i],th_7[,i],th_8[,i],th_9[,i]),2,quantile,c(0.025,.5,.975))
   y.vals=surv.fun
  
  lines(t[x[c(1,2,2)]]*36,y.vals[1,c(1,1,2)],col='gray')
  lines(t[x[c(2,3,3)]]*36,y.vals[1,c(2,2,3)],lty='dotted',col='gray')
  lines(t[x[c(3,4,4)]]*36,y.vals[1,c(3,3,4)],col='gray')
  lines(t[x[c(4,5,5)]]*36,y.vals[1,c(4,4,5)],col='gray')
  lines(t[x[c(5,6,6)]]*36,y.vals[1,c(5,5,6)],lty='dotted',col='gray')
  lines(t[x[c(6,7,7)]]*36,y.vals[1,c(6,6,7)],col='gray')
  lines(t[x[c(7,8,8)]]*36,y.vals[1,c(7,7,8)],col='gray')
  lines(t[x[c(8,9,9)]]*36,y.vals[1,c(8,8,9)],lty='dotted',col='gray')
  lines(t[x[c(9,10,10)]]*36,y.vals[1,c(9,9,10)],col='gray')
  
  lines(t[x[c(1,2,2)]]*36,y.vals[3,c(1,1,2)],col='gray')
  lines(t[x[c(2,3,3)]]*36,y.vals[3,c(2,2,3)],lty='dotted',col='gray')
  lines(t[x[c(3,4,4)]]*36,y.vals[3,c(3,3,4)],col='gray')
  lines(t[x[c(4,5,5)]]*36,y.vals[3,c(4,4,5)],col='gray')
  lines(t[x[c(5,6,6)]]*36,y.vals[3,c(5,5,6)],lty='dotted',col='gray')
  lines(t[x[c(6,7,7)]]*36,y.vals[3,c(6,6,7)],col='gray')
  lines(t[x[c(7,8,8)]]*36,y.vals[3,c(7,7,8)],col='gray')
  lines(t[x[c(8,9,9)]]*36,y.vals[3,c(8,8,9)],lty='dotted',col='gray')
  lines(t[x[c(9,10,10)]]*36,y.vals[3,c(9,9,10)],col='gray')
  
  lines(t[x[c(1,2,2)]]*36,y.vals[2,c(1,1,2)])
  lines(t[x[c(2,3,3)]]*36,y.vals[2,c(2,2,3)],lty='dotted')
  lines(t[x[c(3,4,4)]]*36,y.vals[2,c(3,3,4)])
  lines(t[x[c(4,5,5)]]*36,y.vals[2,c(4,4,5)])
  lines(t[x[c(5,6,6)]]*36,y.vals[2,c(5,5,6)],lty='dotted')
  lines(t[x[c(6,7,7)]]*36,y.vals[2,c(6,6,7)])
  lines(t[x[c(7,8,8)]]*36,y.vals[2,c(7,7,8)])
  lines(t[x[c(8,9,9)]]*36,y.vals[2,c(8,8,9)],lty='dotted')
  lines(t[x[c(9,10,10)]]*36,y.vals[2,c(9,9,10)])
  
  surv.fun_persistence.viability=apply(cbind(phi_0[,i],phi_1[,i],phi_2[,i],phi_3[,i],phi_4[,i],phi_5[,i],phi_6[,i],phi_7[,i],phi_8[,i],phi_9[,i]),2,quantile,c(0.025,.5,.975))
  y.vals=surv.fun_persistence.viability
  
  lines(t[x[c(1,2,2)]]*36,y.vals[1,c(1,1,2)],col='orange')
  lines(t[x[c(2,3,3)]]*36,y.vals[1,c(2,2,3)],lty='dotted',col='orange')
  lines(t[x[c(3,4,4)]]*36,y.vals[1,c(3,3,4)],col='orange')
  lines(t[x[c(4,5,5)]]*36,y.vals[1,c(4,4,5)],col='orange')
  lines(t[x[c(5,6,6)]]*36,y.vals[1,c(5,5,6)],lty='dotted',col='orange')
  lines(t[x[c(6,7,7)]]*36,y.vals[1,c(6,6,7)],col='orange')
  lines(t[x[c(7,8,8)]]*36,y.vals[1,c(7,7,8)],col='orange')
  lines(t[x[c(8,9,9)]]*36,y.vals[1,c(8,8,9)],lty='dotted',col='orange')
  lines(t[x[c(9,10,10)]]*36,y.vals[1,c(9,9,10)],col='orange')
  
  lines(t[x[c(1,2,2)]]*36,y.vals[3,c(1,1,2)],col='orange')
  lines(t[x[c(2,3,3)]]*36,y.vals[3,c(2,2,3)],lty='dotted',col='orange')
  lines(t[x[c(3,4,4)]]*36,y.vals[3,c(3,3,4)],col='orange')
  lines(t[x[c(4,5,5)]]*36,y.vals[3,c(4,4,5)],col='orange')
  lines(t[x[c(5,6,6)]]*36,y.vals[3,c(5,5,6)],lty='dotted',col='orange')
  lines(t[x[c(6,7,7)]]*36,y.vals[3,c(6,6,7)],col='orange')
  lines(t[x[c(7,8,8)]]*36,y.vals[3,c(7,7,8)],col='orange')
  lines(t[x[c(8,9,9)]]*36,y.vals[3,c(8,8,9)],lty='dotted',col='orange')
  lines(t[x[c(9,10,10)]]*36,y.vals[3,c(9,9,10)],col='orange')
  
  lines(t[x[c(1,2,2)]]*36,y.vals[2,c(1,1,2)],col='orange')
  lines(t[x[c(2,3,3)]]*36,y.vals[2,c(2,2,3)],lty='dotted',col='orange')
  lines(t[x[c(3,4,4)]]*36,y.vals[2,c(3,3,4)],col='orange')
  lines(t[x[c(4,5,5)]]*36,y.vals[2,c(4,4,5)],col='orange')
  lines(t[x[c(5,6,6)]]*36,y.vals[2,c(5,5,6)],lty='dotted',col='orange')
  lines(t[x[c(6,7,7)]]*36,y.vals[2,c(6,6,7)],col='orange')
  lines(t[x[c(7,8,8)]]*36,y.vals[2,c(7,7,8)],col='orange')
  lines(t[x[c(8,9,9)]]*36,y.vals[2,c(8,8,9)],lty='dotted',col='orange')
  lines(t[x[c(9,10,10)]]*36,y.vals[2,c(9,9,10)],col='orange')
 
  text(x=3,y=.05,siteNames[i])
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}



par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  plot(NA,NA,xlim=c(0,36),ylim=c(0,1),ylab='',xlab='',xaxt='n',yaxt='n')
  surv.fun=apply(cbind(th_0[,i],th_1[,i],th_2[,i],th_3[,i],th_4[,i],th_5[,i],th_6[,i],th_7[,i],th_8[,i],th_9[,i]),2,quantile,c(0.025,.5,.975))
  y.vals=surv.fun
  
  # build credible interval region
  x.ci = t[x[c(2,3)]]*36+c(0,.5)
  polygon(c(x.ci,rev(x.ci)), c(rep(y.vals[1,2],2),rev(rep(y.vals[3,2],2))) ,
          col="gray90",border=NA)
  
  x.ci = t[x[c(3,4)]]*36+c(.5,0)
  polygon(c(x.ci,rev(x.ci)), c(rep(y.vals[1,3],2),rev(rep(y.vals[3,3],2))) ,
          col="gray90",border=NA)
  
  x.ci = t[x[c(4,5)]]*36+c(0,0)
  polygon(c(x.ci,rev(x.ci)), c(rep(y.vals[1,4],2),rev(rep(y.vals[3,4],2))) ,
          col="gray90",border=NA)
  
  x.ci = t[x[c(5,6)]]*36+c(.5,0)
  polygon(c(x.ci,rev(x.ci)), c(rep(y.vals[1,5],2),rev(rep(y.vals[3,5],2))) ,
          col="gray90",border=NA)
  
  x.ci = t[x[c(6,7)]]*36+c(.5,0)
  polygon(c(x.ci,rev(x.ci)), c(rep(y.vals[1,6],2),rev(rep(y.vals[3,6],2))) ,
          col="gray90",border=NA)
  
  x.ci = t[x[c(7,8)]]*36+c(0,0)
  polygon(c(x.ci,rev(x.ci)), c(rep(y.vals[1,7],2),rev(rep(y.vals[3,7],2))) ,
          col="gray90",border=NA)
  
  x.ci = t[x[c(8,9)]]*36+c(.5,0)
  polygon(c(x.ci,rev(x.ci)), c(rep(y.vals[1,8],2),rev(rep(y.vals[3,8],2))) ,
          col="gray90",border=NA)
  
  x.ci = t[x[c(9,10)]]*36+c(.5,0)
  polygon(c(x.ci,rev(x.ci)), c(rep(y.vals[1,9],2),rev(rep(y.vals[3,9],2))) ,
          col="gray90",border=NA)
  
  x.ci = c(t[x[c(10)]]*36,37)
  polygon(c(x.ci,rev(x.ci)), c(rep(y.vals[1,10],2),rev(rep(y.vals[3,10],2))) ,
          col="gray90",border=NA)
  
  lines(t[x[c(1,2,2)]]*36,y.vals[2,c(1,1,2)])
  lines(t[x[c(2,3)]]*36+c(0,.5),y.vals[2,c(2,2)])

  lines(t[x[c(3,3)]]*36+.5,y.vals[2,c(2,3)],lty='dotted')
  lines(t[x[c(3,4,4)]]*36+c(.5,0,0),y.vals[2,c(3,3,4)])

  lines(t[x[c(4,5,5)]]*36,y.vals[2,c(4,4,5)])
  
  lines(t[x[c(5,6)]]*36+c(0,.5),y.vals[2,c(5,5)])
  lines(t[x[c(5,5)]]*36+.5,y.vals[2,c(5,6)],lty='dotted')
  lines(t[x[c(6,7,7)]]*36+c(.5,0,0),y.vals[2,c(6,6,7)])
  
  lines(t[x[c(7,8,8)]]*36,y.vals[2,c(7,7,8)])

  lines(t[x[c(8,9)]]*36+c(0,.5),y.vals[2,c(8,8)])
  lines(t[x[c(8,8)]]*36+.5,y.vals[2,c(8,9)],lty='dotted')
  lines(t[x[c(9,10,10)]]*36+c(.5,0,0),y.vals[2,c(9,9,10)])

  lines(x.ci,y.vals[2,c(10,10)])
  
  
  
  surv.fun_persistence.viability=apply(cbind(phi_0[,i],phi_1[,i],phi_2[,i],phi_3[,i],phi_4[,i],phi_5[,i],phi_6[,i],phi_7[,i],phi_8[,i],phi_9[,i]),2,quantile,c(0.025,.5,.975))
  y.vals=surv.fun_persistence.viability
 
  # build credible interval region
  x.ci = t[x[c(2,3)]]*36+c(0,.5)
  polygon(c(x.ci,rev(x.ci)), c(rep(y.vals[1,2],2),rev(rep(y.vals[3,2],2))) ,
          col=rgb(1,.647,0,.5),border=NA)
  
  x.ci = t[x[c(3,4)]]*36+c(.5,0)
  polygon(c(x.ci,rev(x.ci)), c(rep(y.vals[1,3],2),rev(rep(y.vals[3,3],2))) ,
          col=rgb(1,.647,0,.5),border=NA)
  
  x.ci = t[x[c(4,5)]]*36+c(0,0)
  polygon(c(x.ci,rev(x.ci)), c(rep(y.vals[1,4],2),rev(rep(y.vals[3,4],2))) ,
          col=rgb(1,.647,0,.5),border=NA)
  
  x.ci = t[x[c(5,6)]]*36+c(.5,0)
  polygon(c(x.ci,rev(x.ci)), c(rep(y.vals[1,5],2),rev(rep(y.vals[3,5],2))) ,
          col=rgb(1,.647,0,.5),border=NA)
  
  x.ci = t[x[c(6,7)]]*36+c(.5,0)
  polygon(c(x.ci,rev(x.ci)), c(rep(y.vals[1,6],2),rev(rep(y.vals[3,6],2))) ,
          col=rgb(1,.647,0,.5),border=NA)
  
  x.ci = t[x[c(7,8)]]*36+c(0,0)
  polygon(c(x.ci,rev(x.ci)), c(rep(y.vals[1,7],2),rev(rep(y.vals[3,7],2))) ,
          col=rgb(1,.647,0,.5),border=NA)
  
  x.ci = t[x[c(8,9)]]*36+c(.5,0)
  polygon(c(x.ci,rev(x.ci)), c(rep(y.vals[1,8],2),rev(rep(y.vals[3,8],2))) ,
          col=rgb(1,.647,0,.5),border=NA)
  
  x.ci = t[x[c(9,10)]]*36+c(.5,0)
  polygon(c(x.ci,rev(x.ci)), c(rep(y.vals[1,9],2),rev(rep(y.vals[3,9],2))) ,
          col=rgb(1,.647,0,.5),border=NA)
  
  x.ci = c(t[x[c(10)]]*36,37)
  polygon(c(x.ci,rev(x.ci)), c(rep(y.vals[1,10],2),rev(rep(y.vals[3,10],2))) ,
          col=rgb(1,.647,0,.5),border=NA)
  
  lines(t[x[c(1,2,2)]]*36,y.vals[2,c(1,1,2)],col='orange')
  lines(t[x[c(2,3)]]*36+c(0,.5),y.vals[2,c(2,2)],col='orange')
  
  lines(t[x[c(3,3)]]*36+.5,y.vals[2,c(2,3)],lty='dotted',col='orange')
  lines(t[x[c(3,4,4)]]*36+c(.5,0,0),y.vals[2,c(3,3,4)],col='orange')
  
  lines(t[x[c(4,5,5)]]*36,y.vals[2,c(4,4,5)],col='orange')
  
  lines(t[x[c(5,6)]]*36+c(0,.5),y.vals[2,c(5,5)],col='orange')
  lines(t[x[c(5,5)]]*36+.5,y.vals[2,c(5,6)],lty='dotted',col='orange')
  lines(t[x[c(6,7,7)]]*36+c(.5,0,0),y.vals[2,c(6,6,7)],col='orange')
  
  lines(t[x[c(7,8,8)]]*36,y.vals[2,c(7,7,8)],col='orange')
  
  lines(t[x[c(8,9)]]*36+c(0,.5),y.vals[2,c(8,8)],col='orange')
  lines(t[x[c(8,8)]]*36+.5,y.vals[2,c(8,9)],lty='dotted',col='orange')
  lines(t[x[c(9,10,10)]]*36+c(.5,0,0),y.vals[2,c(9,9,10)],col='orange')
  
  lines(x.ci,y.vals[2,c(10,10)],col='orange')
  
  
  text(x=-2,y=.05,siteNames[i],pos=4)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
  ifelse(i%in%c(5), legend(x = 15, y = 1,
                           col = c('gray','orange'),
                           lty = c(1,1),
                           legend = c("Persistence only","Persistence & viability"),
                           cex=.55,
                           box.lty=0), NA)
}
mtext("Time (months)", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of persistence and viability", side = 2, outer = TRUE, line = 2.2)
mtext("Population-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)

dev.off()


# -------------------------------------------------------------------
# Calculate structured model parameters
# -------------------------------------------------------------------

# look at germination, conditional on persistence vs. conditional on persistence and viability

# conditional on persistence only
g1=gamma1
g2=gamma2
g3=gamma3

# conditional on persistence and viability
g1.v=gamma1.v
g2.v=gamma2.v
g3.v=gamma3.v

saveRDS(g1.v,"/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/g1-pop.RDS")
saveRDS(g2.v,"/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/g2-pop.RDS")
saveRDS(g3.v,"/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/g3-pop.RDS")

# compare differences
# unconditional calculation from both persistence and viability
hist(th_1*gamma1-phi_1*gamma1.v)
hist(th_4*gamma2-phi_4*gamma2.v)
hist(th_7*gamma3-phi_7*gamma3.v)

dev.off()

pdf("~/Dropbox/clarkiaSeedBanks/products/figures/compare-structured-germination.pdf",width=8,height=4)

# age 2 and 3 unconditional germination estimates are lower than expected for viability 
# when looking at full distribution
par(mfrow=c(1,3))

# plot comparing unconditional probabilities of germination
plot(apply(th_1*gamma1,2,median),apply(phi_1*gamma1.v,2,median),
     frame=FALSE,xlab="",ylab="",pch=21,xlim=c(0,.2),ylim=c(0,.2));
abline(a=0,b=1);
mtext("Unconditional germination probability; \n calculated from persistence & viability",
      side=2,line=2,adj=.5,col='black',cex=.75)
mtext("Age 1",
      side=3,line=0,adj=0,col='black',cex=1)

plot(apply(th_4*gamma2,2,median),apply(phi_4*gamma2.v,2,median),
     frame=FALSE,xlab="",ylab="",pch=21,xlim=c(0,.2),ylim=c(0,.2));
abline(a=0,b=1)
mtext("Unconditional germination probability; calculated from persistence only",
      side=1,line=2.5,adj=.5,col='black',cex=1)
mtext("Age 2",
      side=3,line=0,adj=0,col='black',cex=1)
plot(apply(th_7*gamma3,2,median),apply(phi_7*gamma3.v,2,median),
     frame=FALSE,xlab="",ylab="",pch=21,xlim=c(0,.2),ylim=c(0,.2));
abline(a=0,b=1)
mtext("Age 3",
      side=3,line=0,adj=0,col='black',cex=1)


par(mfrow=c(1,3))

# plot comparing conditional probabilities of germination
plot(apply(g1,2,median),apply(g1.v,2,median),
     frame=FALSE,xlab="",ylab="",pch=21,xlim=c(0,.52),ylim=c(0,.52));
abline(a=0,b=1);
mtext("Germination conditional on persistence & viability",
      side=2,line=2,adj=.5,col='black',cex=.75)
mtext("Age 1",
      side=3,line=0,adj=0,col='black',cex=1)

plot(apply(g2,2,median),apply(g2.v,2,median),
     frame=FALSE,xlab="",ylab="",pch=21,xlim=c(0,.52),ylim=c(0,.52));
abline(a=0,b=1)
mtext("Germination conditional on persistence only",
      side=1,line=2.5,adj=.5,col='black',cex=1)
mtext("Age 2",
      side=3,line=0,adj=0,col='black',cex=1)

plot(apply(g3,2,median),apply(g3.v,2,median),
     frame=FALSE,xlab="",ylab="",pch=21,xlim=c(0,.52),ylim=c(0,.52));
abline(a=0,b=1)
mtext("Age 3",
      side=3,line=0,adj=0,col='black',cex=1)

dev.off()


pdf("~/Dropbox/clarkiaSeedBanks/products/figures/compare-structured-germination2.pdf",width=8,height=6)

g.sum = rbind(apply(g1,2,quantile,c(.025,.5,.975)),apply(g2,2,quantile,c(.025,.5,.975)),apply(g3,2,quantile,c(.025,.5,.975)))
g.v.sum = rbind(apply(g1.v,2,quantile,c(.025,.5,.975)),apply(g2.v,2,quantile,c(.025,.5,.975)),apply(g3.v,2,quantile,c(.025,.5,.975)))

# medians only
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  plot(NA,NA,xlim=c(.5,3.5),ylim=c(0,.8),ylab='',xlab='',xaxt='n',yaxt='n')
  lines(c(1,2,3),g.v.sum[c(2,5,8),i],type='b',col='orange',pch=19)

  lines(c(1,2,3),g.sum[c(2,5,8),i],type='b',col='black',pch=21)

  text(x=.25,y=.75,siteNames[i],pos=4)
  ifelse(i%in%c(16:20),axis(1,c(1,2,3), col.ticks = 1),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Age (years)", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of germination", side = 2, outer = TRUE, line = 2.2)
mtext("Population-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)


# medians+ci
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  plot(NA,NA,xlim=c(.5,3.5),ylim=c(0,.8),ylab='',xlab='',xaxt='n',yaxt='n')
  lines(c(1,2,3)-.1,g.v.sum[c(2,5,8),i],type='b',col='orange',pch=19)
  segments(x0=c(1,2,3)-.1,y0=g.v.sum[c(2,5,8)-1,i],y1=g.v.sum[c(2,5,8)+1,i],col='orange',pch=19)
  
  segments(x0=c(1,2,3)+.1,y0=g.sum[c(2,5,8)-1,i],y1=g.sum[c(2,5,8)+1,i],col='black',lty='dotted')
  lines(c(1,2,3)+.1,g.sum[c(2,5,8),i],type='b',col='black',pch=21,bg='white')
  
  text(x=.25,y=.75,siteNames[i],pos=4)
  ifelse(i%in%c(16:20),axis(1,c(1,2,3), col.ticks = 1),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Age (years)", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of germination", side = 2, outer = TRUE, line = 2.2)
mtext("Population-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)

dev.off()

# survival parameters
s1=th_1
s2=th_3/th_2
s3=th_4/th_3
s4=th_6/th_5
s5=th_7/th_6
s6=th_9/th_8

s1.v=phi_1
s2.v=phi_3/phi_2
s3.v=phi_4/phi_3
s4.v=phi_6/phi_5
s5.v=phi_7/phi_6
s6.v=phi_9/phi_8

saveRDS(s1.v,"/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/s1-pop.RDS")
saveRDS(s2.v,"/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/s2-pop.RDS")
saveRDS(s3.v,"/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/s3-pop.RDS")
saveRDS(s4.v,"/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/s4-pop.RDS")
saveRDS(s5.v,"/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/s5-pop.RDS")
saveRDS(s6.v,"/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/s6-pop.RDS")

## GET S0

mu0_s0=MCMCchains(mcmcSamples,params="mu0_s0")
s0 = boot::inv.logit(mu0_s0)
saveRDS(s0,"/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/s0-pop.RDS")


yr1=s0*s1.v*(1-g1.v)*s2.v
yr2=s0*s1.v*(1-g1.v)*s2.v*s3.v*(1-g2.v)*s4.v
yr3=s0*s1.v*(1-g1.v)*s2.v*s3.v*(1-g2.v)*s4.v*s5.v*(1-g3.v)*s6.v

par(mfrow=c(1,2))
plot(apply(yr1,2,median),apply(yr2,2,median),xlim=c(0,.5),ylim=c(0,.5));abline(a=0,b=1)
text(.05,.45,signif(cor(apply(yr1,2,median),apply(yr2,2,median)),2))
plot(apply(yr1,2,median),apply(yr3,2,median),xlim=c(0,.5),ylim=c(0,.5));abline(a=0,b=1)
text(.05,.45,signif(cor(apply(yr1,2,median),apply(yr3,2,median)),2))

pdf("~/Dropbox/clarkiaSeedBanks/products/figures/compare-structured-survival.pdf",width=8,height=6)

par(mfrow=c(2,3))
# plot comparing conditional probabilities of survival
plot(apply(s1,2,median),apply(s1.v,2,median),
     frame=FALSE,xlab="",ylab="",pch=21,xlim=c(.4,1),ylim=c(.4,1));
abline(a=0,b=1);
mtext("Survival parameter conditional on persistence & viability",
      side=2,line=2,adj=1,col='black',cex=1)
mtext("s1",
      side=3,line=0,adj=0,col='black',cex=1)

plot(apply(s2,2,median),apply(s2.v,2,median),
     frame=FALSE,xlab="",ylab="",pch=21,xlim=c(.4,1),ylim=c(.4,1));
abline(a=0,b=1)
mtext("s2",
      side=3,line=0,adj=0,col='black',cex=1)

plot(apply(s3,2,median),apply(s3.v,2,median),
     frame=FALSE,xlab="",ylab="",pch=21,xlim=c(.4,1),ylim=c(.4,1));
abline(a=0,b=1)
mtext("s3",
      side=3,line=0,adj=0,col='black',cex=1)

plot(apply(s4,2,median),apply(s4.v,2,median),
     frame=FALSE,xlab="",ylab="",pch=21,xlim=c(.4,1),ylim=c(.4,1));
abline(a=0,b=1)
mtext("s4",
      side=3,line=0,adj=0,col='black',cex=1)

plot(apply(s5,2,median),apply(s5.v,2,median),
     frame=FALSE,xlab="",ylab="",pch=21,xlim=c(.4,1),ylim=c(.4,1));
abline(a=0,b=1)
mtext("s5",
      side=3,line=0,adj=0,col='black',cex=1)

mtext("Survival parameter conditional on persistence only",
      side=1,line=2.5,adj=.5,col='black',cex=1)

plot(apply(s6,2,median),apply(s6.v,2,median),
     frame=FALSE,xlab="",ylab="",pch=21,xlim=c(.4,1),ylim=c(.4,1));
abline(a=0,b=1)
mtext("s6",
      side=3,line=0,adj=0,col='black',cex=1)

dev.off()

# compare patterns
s.sum = rbind(apply(s1,2,quantile,c(.025,.5,.975)),
              apply(s2,2,quantile,c(.025,.5,.975)),
              apply(s3,2,quantile,c(.025,.5,.975)),
              apply(s4,2,quantile,c(.025,.5,.975)),
              apply(s5,2,quantile,c(.025,.5,.975)),
              apply(s6,2,quantile,c(.025,.5,.975)))
s.v.sum = rbind(apply(s1.v,2,quantile,c(.025,.5,.975)),
              apply(s2.v,2,quantile,c(.025,.5,.975)),
              apply(s3.v,2,quantile,c(.025,.5,.975)),
              apply(s4.v,2,quantile,c(.025,.5,.975)),
              apply(s5.v,2,quantile,c(.025,.5,.975)),
              apply(s6.v,2,quantile,c(.025,.5,.975)))

pdf("~/Dropbox/clarkiaSeedBanks/products/figures/compare-structured-survival2.pdf",width=8,height=6)

# medians only
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  plot(NA,NA,xlim=c(.5,6.5),ylim=c(.4,1),ylab='',xlab='',xaxt='n',yaxt='n')
  lines(c(1,2,3,4,5,6),s.v.sum[c(2,2+c(3,6,9,12,15)),i],type='b',col='orange',pch=19)
  
  lines(c(1,2,3,4,5,6),s.sum[c(2,2+c(3,6,9,12,15)),i],type='b',col='black',pch=21)
  
  text(x=.25,y=.95,siteNames[i],pos=4)
  ifelse(i%in%c(16:20),axis(1,c(1,2,3,4,5,6), col.ticks = 1),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Survival parameter", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of survival", side = 2, outer = TRUE, line = 2.2)
mtext("Population-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)


# medians+ci
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  plot(NA,NA,xlim=c(.5,6.5),ylim=c(0,1),ylab='',xlab='',xaxt='n',yaxt='n')
  lines(c(1,2,3,4,5,6)-.1,s.v.sum[c(2,2+c(3,6,9,12,15)),i],type='b',col='orange',pch=19)
  segments(x0=c(1,2,3,4,5,6)-.1,y0=s.v.sum[c(2,2+c(3,6,9,12,15))-1,i],y1=s.v.sum[c(2,2+c(3,6,9,12,15))+1,i],col='orange',pch=19)
  
  segments(x0=c(1,2,3,4,5,6)+.1,y0=s.sum[c(2,2+c(3,6,9,12,15))-1,i],y1=s.sum[c(2,2+c(3,6,9,12,15))+1,i],col='black',lty='dotted')
  lines(c(1,2,3,4,5,6)+.1,s.sum[c(2,2+c(3,6,9,12,15)),i],type='b',col='black',pch=21,bg='white')
  
  text(x=.25,y=.95,siteNames[i],pos=4)
  ifelse(i%in%c(16:20),axis(1,c(1,2,3,4,5,6), col.ticks = 1),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Survival parameter", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of survival", side = 2, outer = TRUE, line = 2.2)
mtext("Population-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)

dev.off()

# site 1 is a good example where viability increases in the last
# time step leading to increased survival relative to observed from
# persistence alone
# histograms show this with distribution but hard to summarise, try something else
# par(mfrow=c(2,3))
# 
# hist(s1[,1],col='gray',breaks=100)
# hist(s1.v[,1],col='orange',breaks=100,add=TRUE)
# 
# hist(s2[,1],col='gray',breaks=100)
# hist(s2.v[,1],col='orange',breaks=100,add=TRUE)
# 
# hist(s3[,1],col='gray',breaks=100)
# hist(s3.v[,1],col='orange',breaks=100,add=TRUE)
# 
# hist(s4[,1],col='gray',breaks=100)
# hist(s4.v[,1],col='orange',breaks=100,add=TRUE)
# 
# hist(s5[,1],col='gray',breaks=100)
# hist(s5.v[,1],col='orange',breaks=100,add=TRUE)
# 
# hist(s6[,1],col='gray',breaks=100)
# hist(s6.v[,1],col='orange',breaks=100,add=TRUE)

# could show full distribution of all transitions?
# par(mfrow = c(4,5),
#     oma = c(5,4,0,0) + 0.1,
#     mar = c(0,0,1,1) + 0.1)
# for(i in 1:20){
#   hist(s3[,i],col='gray',breaks=50,ylab='',xlab='',xaxt='n',yaxt='n',main='',freq=FALSE)
#   hist(s3.v[,i],col='orange',breaks=50,add=TRUE,freq=FALSE)
#   text(x=.65,y=6,siteNames[i])
#   ifelse(i%in%c(16:20),axis(1L),NA)
#   ifelse(i%in%c(1,6,11,16),axis(2L),NA)
# }
# mtext("Age (years)", side = 1, outer = TRUE, line = 2.2)
# mtext("Probability of germination", side = 2, outer = TRUE, line = 2.2)
# mtext("Population-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)



# PLOTS TO MAKE
# population level probability of germination conditional on viability and persistence (Figure 7)

# function to plot estimates for both by population
f = function(x,x.v){
  
  x.sum=apply(x,2,quantile,probs=c(0.025,.25,.5,.75,.975))
  x.v.xum=apply(x.v,2,quantile,probs=c(0.025,.25,.5,.75,.975))
  
  plot(NA,NA,type='n',xlim=c(0,1),ylim=c(0,20),
       axes=FALSE,frame=FALSE,
       xlab="",ylab="")
  y.pt = 20:1
  for(i in 20:1){
    tmp<-x.sum[,i]
    segments(x0=tmp[1],x1=tmp[5],y0=y.pt[i]-.2)
    segments(x0=tmp[2],x1=tmp[4],y0=y.pt[i]-.2,lwd=3)
    points(x=tmp[3],y=y.pt[i]-.2,pch=21,bg='white')
    
    tmp2<-x.v.xum[,i]
    segments(x0=tmp2[1],x1=tmp2[5],y0=y.pt[i]+.2,col='orange')
    segments(x0=tmp2[2],x1=tmp2[4],y0=y.pt[i]+.2,lwd=3,col='orange')
    points(x=tmp2[3],y=y.pt[i]+.2,pch=21,bg='white',col='orange')
  }
  axis(1,  seq(0,1,by=.2), col.ticks = 1)
  axis(2, (1:20),
       labels = rev(siteNames), las = 1, 
       col = NA, col.ticks = 1, cex.axis = 1)
}

pdf("~/Dropbox/clarkiaSeedBanks/products/figures/compare-structured-parameters.pdf",width=8,height=4)

# plot germination
par(mfrow=c(1,3))
f(gamma1,gamma1.v)
mtext("Germination probability",
      side=1,line=2.5,adj=.5,col='black',cex=1)
mtext("Age 1",
      side=3,line=0,adj=0,col='black',cex=1)

f(gamma2,gamma2.v)
mtext("Germination probability",
      side=1,line=2.5,adj=.5,col='black',cex=1)
mtext("Age 2",
      side=3,line=0,adj=0,col='black',cex=1)

f(gamma3,gamma3.v)
mtext("Germination probability",
      side=1,line=2.5,adj=.5,col='black',cex=1)
mtext("Age 2",
      side=3,line=0,adj=0,col='black',cex=1)


# plot survival 

par(mfrow=c(1,3))
f(s1,s1.v)
mtext("Survival probability",
      side=1,line=2.5,adj=.5,col='black',cex=1)
mtext("s1: age-1 winter",
      side=3,line=0,adj=0,col='black',cex=1)

f(s3,s3.v)
mtext("Survival probability",
      side=1,line=2.5,adj=.5,col='black',cex=1)
mtext("s3: age-2 winter",
      side=3,line=0,adj=0,col='black',cex=1)

f(s5,s5.v)
mtext("Survival probability",
      side=1,line=2.5,adj=.5,col='black',cex=1)
mtext("s5: age-3 winter",
      side=3,line=0,adj=0,col='black',cex=1)

# plot survival 

par(mfrow=c(1,3))
f(s2,s2.v)
mtext("Survival probability",
      side=1,line=2.5,adj=.5,col='black',cex=1)
mtext("s2: age-1 summer",
      side=3,line=0,adj=0,col='black',cex=1)

f(s4,s4.v)
mtext("Survival probability",
      side=1,line=2.5,adj=.5,col='black',cex=1)
mtext("s4: age-2 summer",
      side=3,line=0,adj=0,col='black',cex=1)

f(s6,s6.v)
mtext("Survival probability",
      side=1,line=2.5,adj=.5,col='black',cex=1)
mtext("s6: age-3 summer",
      side=3,line=0,adj=0,col='black',cex=1)
dev.off()

# function to plot 1:1 comparison
f.compare = function(x,x.v,parm="g1"){
  
  x.sum=apply(x,2,quantile,probs=c(0.025,.25,.5,.75,.975))
  x.v.sum=apply(x.v,2,quantile,probs=c(0.025,.25,.5,.75,.975))
  
  
  plot(NA,NA,type='n',xlim=c(0,1),ylim=c(0,1),
       axes=FALSE,frame=FALSE,
       xlab="",ylab="")
  abline(a=0,b=1,col='gray',lty='dotted')
  for(i in 1:20){
    tmp<-x.sum[,i]
    tmp2<-x.v.sum[,i]
    
    segments(x0=tmp[1],x1=tmp[5],y0=tmp2[3])
    segments(x0=tmp[2],x1=tmp[4],y0=tmp2[3],lwd=2)
    points(x=tmp[3],y=tmp2[3],pch=21,bg='white')
    
    segments(y0=tmp2[1],y1=tmp2[5],x0=tmp[3],col='black')
    segments(y0=tmp2[2],y1=tmp2[4],x0=tmp[3],lwd=2,col='black')
    points(y=tmp2[3],x=tmp[3],pch=21,bg='white',col='black')
  }
  axis(1,  seq(0,1,by=.2), col.ticks = 1)
  axis(2,  seq(0,1,by=.2), col.ticks = 1)
  mtext(paste("Parameter:",parm),
        side=3,line=0,adj=0,col='black',cex=1)
}

pdf("~/Dropbox/clarkiaSeedBanks/products/figures/compare-structured-parameters-1to1.pdf",width=8,height=4)


par(mfrow=c(1,3))
f.compare(g1,g1.v,parm="g1")
mtext("Probability conditional on persistence & viability",
      side=2,line=2.5,adj=.5,col='black',cex=1)
f.compare(g2,g2.v,parm="g2")
mtext("Probability conditional on persistence only",
      side=1,line=2.5,adj=.5,col='black',cex=1)
f.compare(g3,g3.v,parm="g3")

par(mfrow=c(1,3))
f.compare(s1,s1.v,parm="s1")
mtext("Probability conditional on persistence & viability",
      side=2,line=2.5,adj=.5,col='black',cex=1)
f.compare(s3,s3.v,parm="s3")
mtext("Probability conditional on persistence only",
      side=1,line=2.5,adj=.5,col='black',cex=1)
f.compare(s5,s5.v,parm="s5")

f.compare(s2,s2.v,parm="s2")
mtext("Probability conditional on persistence & viability",
      side=2,line=2.5,adj=.5,col='black',cex=1)
f.compare(s4,s4.v,parm="s4")
mtext("Probability conditional on persistence only",
      side=1,line=2.5,adj=.5,col='black',cex=1)
f.compare(s6,s6.v,parm="s6")

dev.off()

# visualize difference as distribution of differences


# function to calculate and plot differences
f=function(x,x.v,parm = "g1"){
  
  x.sum=apply(x.v-x,2,quantile,probs=c(0.025,.25,.5,.75,.975))
  
  plot(NA,NA,type='n',xlim=c(0,max(x.sum)+.05),ylim=c(0,20),
       axes=FALSE,frame=FALSE,
       xlab="",ylab="")
  y.pt = 20:1
  for(i in 20:1){
    tmp<-x.sum[,i]
    segments(x0=tmp[1],x1=tmp[5],y0=y.pt[i])
    segments(x0=tmp[2],x1=tmp[4],y0=y.pt[i],lwd=3)
    points(x=tmp[3],y=y.pt[i],pch=21,bg='white')
    
  }
  axis(1,  seq(0,max(x.sum)+.05,by=.02), col.ticks = 1)
  axis(2, (1:20),
       labels = rev(siteNames), las = 1, 
       col = NA, col.ticks = 1, cex.axis = 1)
  mtext("Delta",
        side=1,line=2.5,adj=.5,col='black',cex=1)
  mtext(paste("Parameter:",parm),
        side=3,line=0,adj=0,col='black',cex=1)
  
  x.sum.df<-data.frame(t(x.sum),position)
  names(x.sum.df)[1:5] = c("ci.lolo","ci.lo","ci.med","ci.hi","ci.hihi")
  
  plot(NA,NA,type='n',ylim=c(0,max(x.sum)+.05),xlim=c(340,375),
       axes=FALSE,frame=FALSE,
       xlab="",ylab="")
  
  segments(x0=x.sum.df$easting,y0=x.sum.df$ci.lolo,y1=x.sum.df$ci.hihi,lwd=1)
  segments(x0=x.sum.df$easting,y0=x.sum.df$ci.lo,y1=x.sum.df$ci.hi,lwd=3)
  points(x=x.sum.df$easting,y=x.sum.df$ci.med,pch=21,bg='white')
  
  axis(2, seq(0,max(x.sum)+.05,by=.02), col.ticks = 1)
  axis(1, seq(340,375,by=5),
       labels = seq(340,375,by=5), las = 1, 
       col.ticks = 1, cex.axis = 1)
  mtext("Delta",
        side=2,line=2.5,adj=.5,col='black',cex=1)
  mtext("Easting (km)",
        side=1,line=2.5,adj=.5,col='black',cex=1)
  # mtext("Age 1",
  #       side=3,line=0,adj=0,col='black',cex=1)
}

f2=function(x,x.v,parm = "s1"){
  
  x.sum=apply(x.v-x,2,quantile,probs=c(0.025,.25,.5,.75,.975))
  
  plot(NA,NA,type='n',xlim=c(min(signif(x.sum,1))-.05,max(signif(x.sum,1))+.05),ylim=c(0,20),
       axes=FALSE,frame=FALSE,
       xlab="",ylab="")
  abline(v=0,col='gray')
  y.pt = 20:1
  for(i in 20:1){
    tmp<-x.sum[,i]
    segments(x0=tmp[1],x1=tmp[5],y0=y.pt[i])
    segments(x0=tmp[2],x1=tmp[4],y0=y.pt[i],lwd=3)
    points(x=tmp[3],y=y.pt[i],pch=21,bg='white')
  }
  
  axis(1,  seq(min(signif(x.sum,1))-.05,max(signif(x.sum,1))+.05,by=.02), col.ticks = 1)
  axis(2, (1:20),
       labels = rev(siteNames), las = 1, 
       col = NA, col.ticks = 1, cex.axis = 1)
  mtext("Delta",
        side=1,line=2.5,adj=.5,col='black',cex=1)
  mtext(paste("Parameter:",parm),
        side=3,line=0,adj=0,col='black',cex=1)
  
  x.sum.df<-data.frame(t(x.sum),position)
  names(x.sum.df)[1:5] = c("ci.lolo","ci.lo","ci.med","ci.hi","ci.hihi")
  
  plot(NA,NA,type='n',ylim=c(min(signif(x.sum,1))-.05,max(signif(x.sum,1))+.05),xlim=c(340,375),
       axes=FALSE,frame=FALSE,
       xlab="",ylab="")
  abline(h=0,col='gray')
  
  segments(x0=x.sum.df$easting,y0=x.sum.df$ci.lolo,y1=x.sum.df$ci.hihi,lwd=1)
  segments(x0=x.sum.df$easting,y0=x.sum.df$ci.lo,y1=x.sum.df$ci.hi,lwd=3)
  points(x=x.sum.df$easting,y=x.sum.df$ci.med,pch=21,bg='white')
  
  axis(2, seq(min(signif(x.sum,1))-.05,max(signif(x.sum,1))+.05,by=.02), col.ticks = 1)
  axis(1, seq(340,375,by=5),
       labels = seq(340,375,by=5), las = 1, 
       col.ticks = 1, cex.axis = 1)
  mtext("Delta",
        side=2,line=2.5,adj=.5,col='black',cex=1)
  mtext("Easting (km)",
        side=1,line=2.5,adj=.5,col='black',cex=1)
}

pdf("~/Dropbox/clarkiaSeedBanks/products/figures/structured-parameters-delta.pdf",width=8,height=4)

par(mfcol=c(2,3),
    oma = c(5,4,0,0) + 0.1,
    mar = c(2,2,1,1) + 0.1)
f(g1,g1.v,parm="g1")
f(g2,g2.v,parm="g2")
f(g3,g3.v,parm="g3")

f2(s1,s1.v,parm="s1")
f2(s2,s2.v,parm="s2")
f2(s3,s3.v,parm="s3")
f2(s4,s4.v,parm="s4")
f2(s5,s5.v,parm="s5")
f2(s6,s6.v,parm="s6")

dev.off()

# population-level estimates of all parameters, by pop and easting (Figure 8)

f.param = function(x.v,parm="g1"){
  
  x.sum=apply(x.v,2,quantile,probs=c(0.025,.25,.5,.75,.975))
  
  plot(NA,NA,type='n',xlim=c(0,1),ylim=c(0,20),
       axes=FALSE,frame=FALSE,
       xlab="",ylab="")
  #abline(v=0,col='gray')
  y.pt = 20:1
  for(i in 20:1){
    tmp<-x.sum[,i]
    segments(x0=tmp[1],x1=tmp[5],y0=y.pt[i])
    segments(x0=tmp[2],x1=tmp[4],y0=y.pt[i],lwd=3)
    points(x=tmp[3],y=y.pt[i],pch=21,bg='white')
  }
  
  axis(1,  seq(0,1,by=.2), col.ticks = 1)
  axis(2, (1:20),
       labels = rev(siteNames), las = 1, 
       col = NA, col.ticks = 1, cex.axis = 1)
  mtext("Probability",
        side=1,line=2.5,adj=.5,col='black',cex=1)
  mtext(paste("Parameter:",parm),
        side=3,line=0,adj=0,col='black',cex=.65)
  
  x.sum.df<-data.frame(t(x.sum),position)
  names(x.sum.df)[1:5] = c("ci.lolo","ci.lo","ci.med","ci.hi","ci.hihi")
  
  plot(NA,NA,type='n',ylim=seq(0,1),xlim=c(340,375),
       axes=FALSE,frame=FALSE,
       xlab="",ylab="")
  #abline(h=0,col='gray')
  
  segments(x0=x.sum.df$easting,y0=x.sum.df$ci.lolo,y1=x.sum.df$ci.hihi,lwd=1)
  segments(x0=x.sum.df$easting,y0=x.sum.df$ci.lo,y1=x.sum.df$ci.hi,lwd=3)
  points(x=x.sum.df$easting,y=x.sum.df$ci.med,pch=21,bg='white')
  
  axis(2, seq(0,1,by=.2), col.ticks = 1)
  axis(1, seq(340,375,by=5),
       labels = seq(340,375,by=5), las = 1, 
       col.ticks = 1, cex.axis = 1)
  mtext("Probability",
        side=2,line=2.5,adj=.5,col='black',cex=1)
  mtext("Easting (km)",
        side=1,line=2.5,adj=.5,col='black',cex=1)
}

pdf("~/Dropbox/clarkiaSeedBanks/products/figures/structured-parameters-all.pdf",width=8,height=6)

par(mfcol=c(2,3),
    oma = c(5,4,0,0) + 0.1,
    mar = c(2,2,1,1) + 0.1)
f.param(g1.v,parm="g1")
f.param(g2.v,parm="g2")
f.param(g3.v,parm="g3")

par(mfcol=c(2,3),
    oma = c(5,4,0,0) + 0.1,
    mar = c(2,2,1,1) + 0.1)
f.param(s1.v,parm="s1")
f.param(s2.v,parm="s2")
f.param(s3.v,parm="s3")

par(mfcol=c(2,3),
    oma = c(5,4,0,0) + 0.1,
    mar = c(2,2,1,1) + 0.1)
f.param(s4.v,parm="s4")
f.param(s5.v,parm="s5")
f.param(s6.v,parm="s6")

dev.off()

# probability that seed persists and is viable after 1, 2, 3 years

pdf("~/Dropbox/clarkiaSeedBanks/products/figures/cumulative-survival.pdf",width=8,height=6)

par(mfcol=c(2,3),
    oma = c(5,4,0,0) + 0.1,
    mar = c(2,2,1,1) + 0.1)
f.param(s1.v*(1-g1.v)*s2.v,parm="Persistent and viable after 1 year")
f.param(s1.v*(1-g1.v)*s2.v*s3.v*(1-g2.v)*s4.v,parm="Persistent and viable after 2 years")
f.param(s1.v*(1-g1.v)*s2.v*s3.v*(1-g2.v)*s4.v*s5.v*(1-g3.v)*s6.v,parm="Persistent and viable after 3 year")

f.param(s1*(1-g1)*s2,parm="Persistent after 1 year")
f.param(s1*(1-g1)*s2*s3*(1-g2)*s4,parm="Persistent after 2 years")
f.param(s1*(1-g1)*s2*s3*(1-g2)*s4*s5*(1-g3)*s6,parm="Persistent after 3 years")

par(mfrow=c(1,1))
f.compare(s1*(1-g1)*s2,s1.v*(1-g1.v)*s2.v,parm="After 1 year")
mtext("Probability conditional on persistence & viability",
      side=2,line=2.5,adj=.5,col='black',cex=1)
mtext("Probability conditional on persistence only",
      side=1,line=2.5,adj=.5,col='black',cex=1)

f.compare(s1*(1-g1)*s2*s3*(1-g2)*s4,s1.v*(1-g1.v)*s2.v*s3.v*(1-g2.v)*s4.v ,parm="After 2 years")
mtext("Probability conditional on persistence & viability",
      side=2,line=2.5,adj=.5,col='black',cex=1)
mtext("Probability conditional on persistence only",
      side=1,line=2.5,adj=.5,col='black',cex=1)

f.compare(s1*(1-g1)*s2*s3*(1-g2)*s4*s5*(1-g3)*s6,s1.v*(1-g1.v)*s2.v*s3.v*(1-g2.v)*s4.v*s5.v*(1-g3.v)*s6.v,parm="After 3 years")
mtext("Probability conditional on persistence & viability",
      side=2,line=2.5,adj=.5,col='black',cex=1)
mtext("Probability conditional on persistence only",
      side=1,line=2.5,adj=.5,col='black',cex=1)

dev.off()
# 
# f2(x=s1*(1-g1)*s2,x.v=s1.v*(1-g1.v)*s2.v,parm="Persistent and viable after 1 year (delta)")
# f2(x=s1*(1-g1)*s2*s3*(1-g2)*s4,x.v=s1.v*(1-g1.v)*s2.v*s3.v*(1-g2.v)*s4.v,parm="Persistent and viable after 2 years (delta)")
# f2(x=s1*(1-g1)*s2*s3*(1-g2)*s4*s5*(1-g3)*s6,x.v=s1.v*(1-g1.v)*s2.v*s3.v*(1-g2.v)*s4.v*s5.v*(1-g3.v)*s6.v,parm="Persistent and viable after 3 years (delta)")

# estimated proportion of seeds that are viable at each time

