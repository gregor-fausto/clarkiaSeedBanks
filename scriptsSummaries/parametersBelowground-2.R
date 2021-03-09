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
# Get site names and position
# -------------------------------------------------------------------

directory2 = "/Users/Gregor/Dropbox/dataLibrary/workflow/tidyData/"
dataFiles <- paste0(directory2,list.files(directory2))

data2 <- readRDS(dataFiles[[1]])
siteNames = unique(data2$siteBags)

position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,dominant.surface.rock.type) %>%
  dplyr::mutate(easting=easting/1000)

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

# par(mfrow=c(1,1))
# 
# plot(NA,NA,type='n',xlim=c(0,2),ylim=c(0,20))
# for(i in 1:20){
#   tmp<-alpha.sum[,i]
#   segments(x0=tmp[1],x1=tmp[3],y0=i)
#   points(x=tmp[2],y=i,pch=19)
# }

# for(i in 1:20){
#   hist(MCMCchains(mcmcSamples,params="a")[,i],
#        breaks=100,freq=FALSE,xlim=c(0,2));
#   abline(v=1,col='red')
# }

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Characterizing survival curves for populations
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Shape parameter
# -------------------------------------------------------------------
alpha<-MCMCchains(mcmcSamples, params = c("a"))
alpha.sum<-apply(alpha,2,quantile,c(.025,.5,.975))

#pdf("~/Dropbox/clarkiaSeedBanks/products/figures/survival-function-parms-population.pdf",width=8,height=6)

par(mfrow=c(1,2))

plot(NA,NA,type='n',xlim=c(0,1.6),ylim=c(0,20),
     axes=FALSE,frame=FALSE,
     xlab="",ylab="")
abline(v=1,lty='dotted')
y.pt = 20:1
for(i in 20:1){
  tmp<-alpha.sum[,i]
  segments(x0=tmp[1],x1=tmp[3],y0=y.pt[i])
  points(x=tmp[2],y=y.pt[i],pch=19)
}
axis(1, c(0,.5,1,1.5), col.ticks = 1)
axis(2, (1:20),
     labels = rev(siteNames), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
mtext(expression(paste("Shape parameter (",alpha,")")),
      side=1,line=2.5,adj=.5,col='black',cex=1)


a.df<-data.frame(t(alpha.sum),position)
names(a.df)[1:3] = c("ci.lo","ci.med","ci.hi")

plot(NA,NA,type='n',ylim=c(0,1.6),xlim=c(340,375),
     axes=FALSE,frame=FALSE,
     xlab="",ylab="")
abline(h=1,lty='dotted')
#y.pt = 20:1

points(x=a.df$easting,y=a.df$ci.med,pch=19)
segments(x0=a.df$easting,y0=a.df$ci.lo,y1=a.df$ci.hi,lwd=1.5)

axis(2, c(0,.5,1,1.5), col.ticks = 1)
axis(1, seq(340,375,by=5),
     labels = seq(340,375,by=5), las = 1, 
     col.ticks = 1, cex.axis = 1)
mtext(expression(paste("Shape parameter (",alpha,")")),
      side=2,line=2.5,adj=.5,col='black',cex=1)
mtext("Easting (km)",
      side=1,line=2.5,adj=.5,col='black',cex=1)

# -------------------------------------------------------------------
# Inverse scale parameter
# -------------------------------------------------------------------
b0=MCMCchains(mcmcSamples,params="mu0_s")
b0.sum<-apply(b0,2,quantile,c(.025,.5,.975))

par(mfrow=c(1,2))

plot(NA,NA,type='n',xlim=c(-1.5,1.5),ylim=c(0,20),
     axes=FALSE,frame=FALSE,
     xlab="",ylab="")
abline(v=0,lty='dotted')
y.pt = 20:1
for(i in 20:1){
  tmp<-b0.sum[,i]
  segments(x0=tmp[1],x1=tmp[3],y0=y.pt[i])
  points(x=tmp[2],y=y.pt[i],pch=19)
}
axis(1, c(-1.5,-1,-.5,0,.5,1,1.5), col.ticks = 1)
axis(2, (1:20),
     labels = rev(siteNames), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
mtext(expression(paste("Inverse scale parameter (",beta,")")),
      side=1,line=2.5,adj=.5,col='black',cex=1)


b0.df<-data.frame(t(b0.sum),position)
names(b0.df)[1:3] = c("ci.lo","ci.med","ci.hi")

plot(NA,NA,type='n',ylim=c(-1.5,1.5),xlim=c(340,375),
     axes=FALSE,frame=FALSE,
     xlab="",ylab="")
abline(h=0,lty='dotted')
#y.pt = 20:1

points(x=b0.df$easting,y=b0.df$ci.med,pch=19)
segments(x0=b0.df$easting,y0=b0.df$ci.lo,y1=b0.df$ci.hi,lwd=1.5)

axis(2, c(-1.5,-1,-.5,0,.5,1,1.5), col.ticks = 1)
axis(1, seq(340,375,by=5),
     labels = seq(340,375,by=5), las = 1, 
     col.ticks = 1, cex.axis = 1)
mtext(expression(paste("Inverse scale parameter (",beta,")")),
      side=2,line=2.5,adj=.5,col='black',cex=1)
mtext("Easting (km)",
      side=1,line=2.5,adj=.5,col='black',cex=1)
#dev.off()

# -------------------------------------------------------------------
# Continuous component of survival function
# -------------------------------------------------------------------
a.wb=MCMCchains(mcmcSamples,params="a")
b0.wb=MCMCchains(mcmcSamples,params="mu0_s")

inv.b0.wb=(exp(-b0.wb/a.wb))
inv.b0.parm=apply(inv.b0.wb,2,quantile,probs=c(0.025,.5,.975))

a.parm=apply(a.wb,2,quantile,probs=c(0.025,.5,.975))



t=seq(0,max(data$months),.01)

weibullSurvival = function(t,inv.b0,alpha){
  exp(-(t/inv.b0)^alpha)
}
weibullSurvival(0,inv.b0=inv.b0.wb[,1],alpha=a.wb[,1])

#pdf("~/Dropbox/clarkiaSeedBanks/products/figures/survival-function-population.pdf",width=8,height=6)


par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  # calculate survival at each time point using 
  # row samples from posterior
  out=lapply(t,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i])
  
  # create matrix 
  out.list=do.call(cbind,out)
  
  plot(t*36,out.list[1,],type='n',ylim=c(0,1),ylab='',xlab='',xaxt='n',yaxt='n')
  
  index=sample(1:dim(out.list)[1],50)
  for(j in index){
    lines(t*36,out.list[j,], lwd = .75,
          col=ifelse(a.wb[j,i]>1,"purple","orange"))
  }
  text(x=3,y=.05,siteNames[i])
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Time (months)", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of persistence", side = 2, outer = TRUE, line = 2.2)
mtext("Population-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)


par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  # calculate survival at each time point using 
  # row samples from posterior
  out=lapply(t,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i])

    # create matrix 
  out.list=do.call(cbind,out)
  out.list = apply(out.list,2,quantile,probs=c(0.025,.5,.975))
  plot(t*36,out.list[1,],type='n',ylim=c(0,1),ylab='',xlab='',xaxt='n',yaxt='n')
  
  polygon(x=c(t*36,rev(t*36)),y=c(out.list[1,],rev(out.list[3,])),
          col="cornsilk",border=NA)
  lines(t*36,out.list[2,], lwd = 2)

  text(x=3,y=.05,siteNames[i])
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Time (months)", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of persistence", side = 2, outer = TRUE, line = 2.2)
mtext("Population-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)

#dev.off()

# -------------------------------------------------------------------
# Germination probabilities
# -------------------------------------------------------------------
mu0_g=MCMCchains(mcmcSamples,params="mu0_g")
mu_g=MCMCchains(mcmcSamples,params="mu_g")

pdf("~/Dropbox/clarkiaSeedBanks/products/figures/germination-population.pdf",width=8,height=4)

par(mfrow=c(1,3),
    oma = c(5,4,0,0) + 0.1,
    mar = c(2,2,1,1) + 0.1)

gamma1 = boot::inv.logit(mu0_g[,1:20])
gamma1.sum=apply(gamma1,2,quantile,probs=c(0.025,.25,.5,.75,.975))

plot(NA,NA,type='n',xlim=c(0,1),ylim=c(0,20),
     axes=FALSE,frame=FALSE,
     xlab="",ylab="")
y.pt = 20:1
for(i in 20:1){
  tmp<-gamma1.sum[,i]
  segments(x0=tmp[1],x1=tmp[5],y0=y.pt[i])
  segments(x0=tmp[2],x1=tmp[4],y0=y.pt[i],lwd=3)
  
  points(x=tmp[3],y=y.pt[i],pch=21,bg='white')
}
axis(1,  seq(0,1,by=.2), col.ticks = 1)
axis(2, (1:20),
     labels = rev(siteNames), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
mtext("Germination probability",
      side=1,line=2.5,adj=.5,col='black',cex=1)
mtext("Age 1",
      side=3,line=0,adj=0,col='black',cex=1)

gamma2 = boot::inv.logit(mu0_g[,21:40])
gamma2.sum=apply(gamma2,2,quantile,probs=c(0.025,.25,.5,.75,.975))

plot(NA,NA,type='n',xlim=c(0,1),ylim=c(0,20),
     axes=FALSE,frame=FALSE,
     xlab="",ylab="")
y.pt = 20:1
for(i in 20:1){
  tmp<-gamma2.sum[,i]
  segments(x0=tmp[1],x1=tmp[5],y0=y.pt[i])
  segments(x0=tmp[2],x1=tmp[4],y0=y.pt[i],lwd=3)
  
  points(x=tmp[3],y=y.pt[i],pch=21,bg='white')
}
axis(1,  seq(0,1,by=.2), col.ticks = 1)
axis(2, (1:20),
     labels = rev(siteNames), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
mtext("Germination probability",
      side=1,line=2.5,adj=.5,col='black',cex=1)
mtext("Age 2",
      side=3,line=0,adj=0,col='black',cex=1)


# note here this is the population-year level 
# there is NO population-level parameter for age 3 germination
gamma3 = boot::inv.logit(mu_g[,101:120])
gamma3.sum=apply(gamma3,2,quantile,probs=c(0.025,.25,.5,.75,.975))

plot(NA,NA,type='n',xlim=c(0,1),ylim=c(0,20),
     axes=FALSE,frame=FALSE,
     xlab="",ylab="")
y.pt = 20:1
for(i in 20:1){
  tmp<-gamma3.sum[,i]
  segments(x0=tmp[1],x1=tmp[5],y0=y.pt[i])
  segments(x0=tmp[2],x1=tmp[4],y0=y.pt[i],lwd=3)
  
  points(x=tmp[3],y=y.pt[i],pch=21,bg='white')
}
axis(1,  seq(0,1,by=.2), col.ticks = 1)
axis(2, (1:20),
     labels = rev(siteNames), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
mtext("Germination probability",
      side=1,line=2.5,adj=.5,col='black',cex=1)
mtext("Age 3",
      side=3,line=0,adj=0,col='black',cex=1)



par(mfrow=c(1,3),
    oma = c(5,4,0,0) + 0.1,
    mar = c(2,2,1,1) + 0.1)

gamma1.sum<-data.frame(t(gamma1.sum),position)
names(gamma1.sum)[1:5] = c("ci.lolo","ci.lo","ci.med","ci.hi","ci.hihi")

plot(NA,NA,type='n',ylim=c(0,1),xlim=c(340,375),
     axes=FALSE,frame=FALSE,
     xlab="",ylab="")

segments(x0=gamma1.sum$easting,y0=gamma1.sum$ci.lolo,y1=gamma1.sum$ci.hihi,lwd=1)
segments(x0=gamma1.sum$easting,y0=gamma1.sum$ci.lo,y1=gamma1.sum$ci.hi,lwd=3)
points(x=gamma1.sum$easting,y=gamma1.sum$ci.med,pch=21,bg='white')

axis(2, seq(0,1,by=.2), col.ticks = 1)
axis(1, seq(340,375,by=5),
     labels = seq(340,375,by=5), las = 1, 
     col.ticks = 1, cex.axis = 1)
mtext("Germination probability",
      side=2,line=2.5,adj=.5,col='black',cex=1)
mtext("Easting (km)",
      side=1,line=2.5,adj=.5,col='black',cex=1)
mtext("Age 1",
      side=3,line=0,adj=0,col='black',cex=1)


gamma2.sum<-data.frame(t(gamma2.sum),position)
names(gamma2.sum)[1:5] = c("ci.lolo","ci.lo","ci.med","ci.hi","ci.hihi")

plot(NA,NA,type='n',ylim=c(0,1),xlim=c(340,375),
     axes=FALSE,frame=FALSE,
     xlab="",ylab="")

segments(x0=gamma2.sum$easting,y0=gamma2.sum$ci.lolo,y1=gamma2.sum$ci.hihi,lwd=1)
segments(x0=gamma2.sum$easting,y0=gamma2.sum$ci.lo,y1=gamma2.sum$ci.hi,lwd=3)
points(x=gamma2.sum$easting,y=gamma2.sum$ci.med,pch=21,bg='white')

axis(2, seq(0,1,by=.2), col.ticks = 1)
axis(1, seq(340,375,by=5),
     labels = seq(340,375,by=5), las = 1, 
     col.ticks = 1, cex.axis = 1)
# mtext("Germination probability",
#       side=2,line=2.5,adj=.5,col='black',cex=1)
mtext("Easting (km)",
      side=1,line=2.5,adj=.5,col='black',cex=1)
mtext("Age 2",
      side=3,line=0,adj=0,col='black',cex=1)

gamma3.sum<-data.frame(t(gamma3.sum),position)
names(gamma3.sum)[1:5] = c("ci.lolo","ci.lo","ci.med","ci.hi","ci.hihi")

plot(NA,NA,type='n',ylim=c(0,1),xlim=c(340,375),
     axes=FALSE,frame=FALSE,
     xlab="",ylab="")

segments(x0=gamma3.sum$easting,y0=gamma3.sum$ci.lolo,y1=gamma3.sum$ci.hihi,lwd=1)
segments(x0=gamma3.sum$easting,y0=gamma3.sum$ci.lo,y1=gamma3.sum$ci.hi,lwd=3)
points(x=gamma3.sum$easting,y=gamma3.sum$ci.med,pch=21,bg='white')

axis(2, seq(0,1,by=.2), col.ticks = 1)
axis(1, seq(340,375,by=5),
     labels = seq(340,375,by=5), las = 1, 
     col.ticks = 1, cex.axis = 1)
# mtext("Germination probability",
#       side=2,line=2.5,adj=.5,col='black',cex=1)
mtext("Easting (km)",
      side=1,line=2.5,adj=.5,col='black',cex=1)
mtext("Age 3",
      side=3,line=0,adj=0,col='black',cex=1)


dev.off()

########################################################
# Discrete component of survival function; age-specific
########################################################
theta_c2 = 1- gamma1
theta_c4 = (1-gamma2)
theta_c6 = (1-gamma3)



theta_c2.sum = apply(theta_c2,2,quantile,probs=c(0.025,.25,.5,.75,.975))
theta_c4.sum = apply(theta_c4,2,quantile,probs=c(0.025,.25,.5,.75,.975))
theta_c6.sum = apply(theta_c6,2,quantile,probs=c(0.025,.25,.5,.75,.975))

pdf("~/Dropbox/clarkiaSeedBanks/products/figures/survival-function-discrete.pdf",width=8,height=6)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)


for(i in 1:20){
  plot(NA,NA,type='n',ylim=c(0,1),
       xlim=c(0,36),ylab='',xlab='',xaxt='n',yaxt='n')
  
  segments(x0=c(4,16,28),
           y0=c(theta_c2.sum[1,i],theta_c4.sum[1,i],theta_c6.sum[1,i]),
           y1=c(theta_c2.sum[5,i],theta_c4.sum[1,i],theta_c6.sum[5,i]))
  segments(x0=c(4,16,28),
           y0=c(theta_c2.sum[2,i],theta_c4.sum[2,i],theta_c6.sum[2,i]),
           y1=c(theta_c2.sum[4,i],theta_c4.sum[4,i],theta_c6.sum[4,i]),lwd=3)
  points(c(4,16,28),c(theta_c2.sum[3,i],theta_c4.sum[3,i],theta_c6.sum[3,i]),
         pch=21,bg='white')
  
  text(x=3,y=.05,siteNames[i])
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Time (months)", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of persistence", side = 2, outer = TRUE, line = 2.2)
mtext("Population-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)

dev.off()

########################################################
# Combined
########################################################
# first, create compound discrete survival function components
theta_c2 = 1- gamma1
theta_c4 = (1- gamma1)*(1-gamma2)
theta_c6 = (1- gamma1)*(1-gamma2)*(1-gamma3)

a.wb=MCMCchains(mcmcSamples,params="a")
b0.wb=MCMCchains(mcmcSamples,params="mu0_s")

inv.b0.wb=(exp(-b0.wb/a.wb))

#inv.b0.parm=apply(inv.b0,2,quantile,probs=c(0.025,.5,.975))
#a.parm=apply(a.wb,2,quantile,probs=c(0.025,.5,.975))

pdf("~/Dropbox/clarkiaSeedBanks/products/figures/survival-function-compound.pdf",width=8,height=6)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  
  t=seq(0,max(data$months),.01)
  
  t1 = seq(0,4/36,length.out=25)
  
  weibullSurvival = function(t,inv.b0,alpha){
    exp(-(t/inv.b0)^alpha)
  }
  out=lapply(t1,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i])
  out.list=do.call(cbind,out)
  
  g1 = 4/36
  
  weibullSurvival = function(t,inv.b0,alpha,theta_c2){
    theta_c2*exp(-(t/inv.b0)^alpha)
  }
  out2=weibullSurvival(g1,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i])
  out.list=cbind(out.list,out2)
  
  t2 = seq(4/36,12/36,length.out=25)
  
  weibullSurvival = function(t,inv.b0,alpha,theta_c2){
    theta_c2*exp(-(t/inv.b0)^alpha)
  }
  out=lapply(t2,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i])
  out2=do.call(cbind,out)
  out.list=cbind(out.list,out2)
  
  t3 = seq(12/36,16/36,length.out=25)
  
  weibullSurvival = function(t,inv.b0,alpha,theta_c2){
    theta_c2*exp(-(t/inv.b0)^alpha)
  }
  out=lapply(t3,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i])
  out2=do.call(cbind,out)
  out.list=cbind(out.list,out2)
  
  g2 = 16/36
  
  weibullSurvival = function(t,inv.b0,alpha,theta_c2,theta_c4){
    theta_c4*theta_c2*exp(-(t/inv.b0)^alpha)
  }
  out2=weibullSurvival(g2,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i],theta_c4[,i])
  out.list=cbind(out.list,out2)
  
  t4 = seq(16/36,24/36,length.out=25)
  
  weibullSurvival = function(t,inv.b0,alpha,theta_c2,theta_c4){
    theta_c4*theta_c2*exp(-(t/inv.b0)^alpha)
  }
  out=lapply(t4,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i],theta_c4[,i])
  out2=do.call(cbind,out)
  out.list=cbind(out.list,out2)
  
  t5 = seq(24/36,28/36,length.out=25)
  
  weibullSurvival = function(t,inv.b0,alpha,theta_c2,theta_c4){
    theta_c4*theta_c2*exp(-(t/inv.b0)^alpha)
  }
  out=lapply(t5,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i],theta_c4[,i])
  out2=do.call(cbind,out)
  out.list=cbind(out.list,out2)
  
  g3 = 28/36
  
  weibullSurvival = function(t,inv.b0,alpha,theta_c2,theta_c4,theta_c6){
    theta_c6*theta_c4*theta_c2*exp(-(t/inv.b0)^alpha)
  }
  out2=weibullSurvival(g3,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i],theta_c4[,i],theta_c6[,i])
  out.list=cbind(out.list,out2)
  
  t6 = seq(28/36,36/36,length.out=25)
  
  weibullSurvival = function(t,inv.b0,alpha,theta_c2,theta_c4,theta_c6){
    theta_c6*theta_c4*theta_c2*exp(-(t/inv.b0)^alpha)
  }
  out=lapply(t6,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i],theta_c4[,i],theta_c6[,i])
  out2=do.call(cbind,out)
  out.list=cbind(out.list,out2)
  
  # plot(NA,NA,type='n',ylim=c(0,1),xlim=c(0,40),ylab='',xlab='',xaxt='n',yaxt='n')
  # 
  # index=sample(1:dim(out.list)[1],50)
  # for(j in index){
  # #for(j in c(1)){
  #   lines(t1*36,out.list[j,1:25], lwd = .75,
  #         col=ifelse(a.wb[j,i]>1,"purple","orange"))
  #   segments(x0=g1*36,y0=out.list[j,25],y1=out.list[j,26],lty='dotted',col='orange')
  #   lines(t2*36,out.list[j,27:51], lwd = .75,
  #         col=ifelse(a.wb[j,i]>1,"purple","orange"))
  #   lines(t3*36,out.list[j,52:76], lwd = .75,
  #         col=ifelse(a.wb[j,i]>1,"purple","orange"))
  #   segments(x0=g2*36,y0=out.list[j,76],y1=out.list[j,77],lty='dotted',col='orange')
  #   lines(t4*36,out.list[j,78:102], lwd = .75,
  #         col=ifelse(a.wb[j,i]>1,"purple","orange"))
  #   lines(t5*36,out.list[j,103:127], lwd = .75,
  #         col=ifelse(a.wb[j,i]>1,"purple","orange"))
  #   segments(x0=g3*36,y0=out.list[j,127],y1=out.list[j,128],lty='dotted',col='orange')
  #   lines(t6*36,out.list[j,129:153], lwd = .75,
  #         col=ifelse(a.wb[j,i]>1,"purple","orange")) 
  # }
  
  #pdf("~/Dropbox/clarkiaSeedBanks/products/figures/survival-function-population.pdf",width=8,height=6)
  
  out.list.sum = apply(out.list,2,quantile,probs=c(0.025,.25,.5,.75,.975))
  
  plot(NA,NA,type='n',ylim=c(0,1),xlim=c(0,36),ylab='',xlab='',xaxt='n',yaxt='n')
  
  tmp = out.list.sum
  
  polygon(x=c(t1*36,rev(t1*36)),y=c(tmp[1,1:25],rev(tmp[5,1:25])),col="cornsilk1",border=NA)
  polygon(x=c(t1*36,rev(t1*36)),y=c(tmp[2,1:25],rev(tmp[4,1:25])),col="cornsilk2",border=NA)
  lines(t1*36,tmp[3,1:25], lwd = .75)
  segments(x0=g1*36,y0=tmp[3,25],y1=tmp[3,26],lty='dotted',col='black')
  
  polygon(x=c(c(t2,t3)*36,rev(c(t2,t3)*36)),y=c(tmp[1,27:76],rev(tmp[5,27:76])),col="cornsilk1",border=NA)
  polygon(x=c(c(t2,t3)*36,rev(c(t2,t3)*36)),y=c(tmp[2,27:76],rev(tmp[4,27:76])),col="cornsilk2",border=NA)
  lines(t2*36,tmp[3,27:51], lwd = .75)
  lines(t3*36,tmp[3,52:76], lwd = .75)
  segments(x0=g2*36,y0=tmp[3,76],y1=tmp[3,77],lty='dotted',col='black')
  
  polygon(x=c(c(t4,t5)*36,rev(c(t4,t5)*36)),y=c(tmp[1,78:127],rev(tmp[5,78:127])),col="cornsilk1",border=NA)
  polygon(x=c(c(t4,t5)*36,rev(c(t4,t5)*36)),y=c(tmp[2,78:127],rev(tmp[4,78:127])),col="cornsilk2",border=NA)
  lines(t4*36,tmp[3,78:102], lwd = .75)
  lines(t5*36,tmp[3,103:127], lwd = .75)
  segments(x0=g3*36,y0=tmp[3,127],y1=tmp[3,128],lty='dotted',col='black')
  
  polygon(x=c(c(t6)*36,rev(c(t6)*36)),y=c(tmp[1,129:153],rev(tmp[5,129:153])),col="cornsilk1",border=NA)
  polygon(x=c(c(t6)*36,rev(c(t6)*36)),y=c(tmp[2,129:153],rev(tmp[4,129:153])),col="cornsilk2",border=NA)
  lines(t6*36,tmp[3,129:153], lwd = .75) 
  
  text(x=3,y=.05,siteNames[i])
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Time (months)", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of persistence", side = 2, outer = TRUE, line = 2.2)
mtext("Population-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)

dev.off()


pdf("~/Dropbox/clarkiaSeedBanks/products/figures/survival-function-summary.pdf",width=8,height=4)

par(mfrow=c(1,3),
    oma = c(5,4,0,0) + 0.1,
    mar = c(2,2,1,1) + 0.1)

plot(NA,NA,type='n',xlim=c(340,375),ylim=c(0,1),
     ylab='',xlab='',xaxt='n',yaxt='n',axes=FALSE,frame=FALSE)
for(i in 20:1){
  
  t=seq(0,max(data$months),.01)
  
  t1 = seq(0,4/36,length.out=25)
  
  weibullSurvival = function(t,inv.b0,alpha){
    exp(-(t/inv.b0)^alpha)
  }
  out=lapply(t1,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i])
  out.list=do.call(cbind,out)
  
  g1 = 4/36

  weibullSurvival = function(t,inv.b0,alpha,theta_c2){
    theta_c2*exp(-(t/inv.b0)^alpha)
  }
  out2=weibullSurvival(g1,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i])
  out.list=cbind(out.list,out2)
  
  t2 = seq(4/36,12/36,length.out=25)
  
  weibullSurvival = function(t,inv.b0,alpha,theta_c2){
    theta_c2*exp(-(t/inv.b0)^alpha)
  }
  out=lapply(t2,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i])
  out2=do.call(cbind,out)
  out.list=cbind(out.list,out2)
  # 
  # t3 = seq(12/36,16/36,length.out=25)
  # 
  # weibullSurvival = function(t,inv.b0,alpha,theta_c2){
  #   theta_c2*exp(-(t/inv.b0)^alpha)
  # }
  # out=lapply(t3,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i])
  # out2=do.call(cbind,out)
  # out.list=cbind(out.list,out2)
  # 
  # g2 = 16/36
  # 
  # weibullSurvival = function(t,inv.b0,alpha,theta_c2,theta_c4){
  #   theta_c4*theta_c2*exp(-(t/inv.b0)^alpha)
  # }
  # out2=weibullSurvival(g2,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i],theta_c4[,i])
  # out.list=cbind(out.list,out2)
  # 
  # t4 = seq(16/36,24/36,length.out=25)
  # 
  # weibullSurvival = function(t,inv.b0,alpha,theta_c2,theta_c4){
  #   theta_c4*theta_c2*exp(-(t/inv.b0)^alpha)
  # }
  # out=lapply(t4,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i],theta_c4[,i])
  # out2=do.call(cbind,out)
  # out.list=cbind(out.list,out2)
  # 
  # t5 = seq(24/36,28/36,length.out=25)
  # 
  # weibullSurvival = function(t,inv.b0,alpha,theta_c2,theta_c4){
  #   theta_c4*theta_c2*exp(-(t/inv.b0)^alpha)
  # }
  # out=lapply(t5,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i],theta_c4[,i])
  # out2=do.call(cbind,out)
  # out.list=cbind(out.list,out2)
  # 
  # g3 = 28/36
  # 
  # weibullSurvival = function(t,inv.b0,alpha,theta_c2,theta_c4,theta_c6){
  #   theta_c6*theta_c4*theta_c2*exp(-(t/inv.b0)^alpha)
  # }
  # out2=weibullSurvival(g3,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i],theta_c4[,i],theta_c6[,i])
  # out.list=cbind(out.list,out2)
  # 
  # t6 = seq(28/36,36/36,length.out=25)
  # 
  # weibullSurvival = function(t,inv.b0,alpha,theta_c2,theta_c4,theta_c6){
  #   theta_c6*theta_c4*theta_c2*exp(-(t/inv.b0)^alpha)
  # }
  # out=lapply(t6,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i],theta_c4[,i],theta_c6[,i])
  # out2=do.call(cbind,out)
  # out.list=cbind(out.list,out2)
  
  out.list.sum = apply(out.list,2,quantile,probs=c(0.025,.25,.5,.75,.975))
  

  segments(x0=position$easting[i],
           y0=out.list.sum[1,dim(out.list.sum)[2]],
           y1=out.list.sum[5,dim(out.list.sum)[2]])
  segments(x0=position$easting[i],
           y0=out.list.sum[2,dim(out.list.sum)[2]],
           y1=out.list.sum[4,dim(out.list.sum)[2]],lwd=3)
  points(position$easting[i],out.list.sum[3,dim(out.list.sum)[2]],pch=21,bg='white')#,
       # col=ifelse(position$dominant.surface.rock.type[i]=="igneous","black","orange"))
}

axis(2, seq(0,1,by=.1), col.ticks = 1)
axis(1, seq(340,375,by=5),
     labels = seq(340,375,by=5), las = 1, 
     col.ticks = 1, cex.axis = 1)
mtext("Probability of seed persistence",
      side=2,line=2.5,adj=.5,col='black',cex=1)
mtext("After 1 year",
      side=3,line=-.5,adj=0,col='black',cex=1)
mtext("Easting (km)",
      side=1,line=2.5,adj=.5,col='black',cex=1)



plot(NA,NA,type='n',xlim=c(340,375),ylim=c(0,1),
     ylab='',xlab='',xaxt='n',yaxt='n',axes=FALSE,frame=FALSE)
for(i in 20:1){
  
  t=seq(0,max(data$months),.01)
  
  t1 = seq(0,4/36,length.out=25)
  
  weibullSurvival = function(t,inv.b0,alpha){
    exp(-(t/inv.b0)^alpha)
  }
  out=lapply(t1,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i])
  out.list=do.call(cbind,out)
  
  g1 = 4/36
  
  weibullSurvival = function(t,inv.b0,alpha,theta_c2){
    theta_c2*exp(-(t/inv.b0)^alpha)
  }
  out2=weibullSurvival(g1,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i])
  out.list=cbind(out.list,out2)
  
  t2 = seq(4/36,12/36,length.out=25)
  
  weibullSurvival = function(t,inv.b0,alpha,theta_c2){
    theta_c2*exp(-(t/inv.b0)^alpha)
  }
  out=lapply(t2,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i])
  out2=do.call(cbind,out)
  out.list=cbind(out.list,out2)

  t3 = seq(12/36,16/36,length.out=25)

  weibullSurvival = function(t,inv.b0,alpha,theta_c2){
    theta_c2*exp(-(t/inv.b0)^alpha)
  }
  out=lapply(t3,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i])
  out2=do.call(cbind,out)
  out.list=cbind(out.list,out2)

  g2 = 16/36

  weibullSurvival = function(t,inv.b0,alpha,theta_c2,theta_c4){
    theta_c4*theta_c2*exp(-(t/inv.b0)^alpha)
  }
  out2=weibullSurvival(g2,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i],theta_c4[,i])
  out.list=cbind(out.list,out2)

  t4 = seq(16/36,24/36,length.out=25)

  weibullSurvival = function(t,inv.b0,alpha,theta_c2,theta_c4){
    theta_c4*theta_c2*exp(-(t/inv.b0)^alpha)
  }
  out=lapply(t4,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i],theta_c4[,i])
  out2=do.call(cbind,out)
  out.list=cbind(out.list,out2)
  # 
  # t5 = seq(24/36,28/36,length.out=25)
  # 
  # weibullSurvival = function(t,inv.b0,alpha,theta_c2,theta_c4){
  #   theta_c4*theta_c2*exp(-(t/inv.b0)^alpha)
  # }
  # out=lapply(t5,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i],theta_c4[,i])
  # out2=do.call(cbind,out)
  # out.list=cbind(out.list,out2)
  # 
  # g3 = 28/36
  # 
  # weibullSurvival = function(t,inv.b0,alpha,theta_c2,theta_c4,theta_c6){
  #   theta_c6*theta_c4*theta_c2*exp(-(t/inv.b0)^alpha)
  # }
  # out2=weibullSurvival(g3,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i],theta_c4[,i],theta_c6[,i])
  # out.list=cbind(out.list,out2)
  # 
  # t6 = seq(28/36,36/36,length.out=25)
  # 
  # weibullSurvival = function(t,inv.b0,alpha,theta_c2,theta_c4,theta_c6){
  #   theta_c6*theta_c4*theta_c2*exp(-(t/inv.b0)^alpha)
  # }
  # out=lapply(t6,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i],theta_c4[,i],theta_c6[,i])
  # out2=do.call(cbind,out)
  # out.list=cbind(out.list,out2)
  
  out.list.sum = apply(out.list,2,quantile,probs=c(0.025,.25,.5,.75,.975))
  

  segments(x0=position$easting[i],
           y0=out.list.sum[1,dim(out.list.sum)[2]],
           y1=out.list.sum[5,dim(out.list.sum)[2]])
  segments(x0=position$easting[i],
           y0=out.list.sum[2,dim(out.list.sum)[2]],
           y1=out.list.sum[4,dim(out.list.sum)[2]],lwd=3)
  points(position$easting[i],out.list.sum[3,dim(out.list.sum)[2]],pch=21,bg='white')#,
       # col=ifelse(position$dominant.surface.rock.type[i]=="igneous","black","orange"))
  
}

axis(2, seq(0,1,by=.1), col.ticks = 1)
axis(1, seq(340,375,by=5),
     labels = seq(340,375,by=5), las = 1, 
     col.ticks = 1, cex.axis = 1)
mtext("After 2 years",
      side=3,line=-.5,adj=0,col='black',cex=1)
mtext("Easting (km)",
      side=1,line=2.5,adj=.5,col='black',cex=1)



plot(NA,NA,type='n',xlim=c(340,375),ylim=c(0,1),
     ylab='',xlab='',xaxt='n',yaxt='n',axes=FALSE,frame=FALSE)
for(i in 20:1){
  
  t=seq(0,max(data$months),.01)
  
  t1 = seq(0,4/36,length.out=25)
  
  weibullSurvival = function(t,inv.b0,alpha){
    exp(-(t/inv.b0)^alpha)
  }
  out=lapply(t1,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i])
  out.list=do.call(cbind,out)
  
  g1 = 4/36
  
  weibullSurvival = function(t,inv.b0,alpha,theta_c2){
    theta_c2*exp(-(t/inv.b0)^alpha)
  }
  out2=weibullSurvival(g1,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i])
  out.list=cbind(out.list,out2)
  
  t2 = seq(4/36,12/36,length.out=25)
  
  weibullSurvival = function(t,inv.b0,alpha,theta_c2){
    theta_c2*exp(-(t/inv.b0)^alpha)
  }
  out=lapply(t2,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i])
  out2=do.call(cbind,out)
  out.list=cbind(out.list,out2)
  
  t3 = seq(12/36,16/36,length.out=25)
  
  weibullSurvival = function(t,inv.b0,alpha,theta_c2){
    theta_c2*exp(-(t/inv.b0)^alpha)
  }
  out=lapply(t3,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i])
  out2=do.call(cbind,out)
  out.list=cbind(out.list,out2)
  
  g2 = 16/36
  
  weibullSurvival = function(t,inv.b0,alpha,theta_c2,theta_c4){
    theta_c4*theta_c2*exp(-(t/inv.b0)^alpha)
  }
  out2=weibullSurvival(g2,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i],theta_c4[,i])
  out.list=cbind(out.list,out2)
  
  t4 = seq(16/36,24/36,length.out=25)
  
  weibullSurvival = function(t,inv.b0,alpha,theta_c2,theta_c4){
    theta_c4*theta_c2*exp(-(t/inv.b0)^alpha)
  }
  out=lapply(t4,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i],theta_c4[,i])
  out2=do.call(cbind,out)
  out.list=cbind(out.list,out2)

  t5 = seq(24/36,28/36,length.out=25)

  weibullSurvival = function(t,inv.b0,alpha,theta_c2,theta_c4){
    theta_c4*theta_c2*exp(-(t/inv.b0)^alpha)
  }
  out=lapply(t5,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i],theta_c4[,i])
  out2=do.call(cbind,out)
  out.list=cbind(out.list,out2)

  g3 = 28/36

  weibullSurvival = function(t,inv.b0,alpha,theta_c2,theta_c4,theta_c6){
    theta_c6*theta_c4*theta_c2*exp(-(t/inv.b0)^alpha)
  }
  out2=weibullSurvival(g3,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i],theta_c4[,i],theta_c6[,i])
  out.list=cbind(out.list,out2)

  t6 = seq(28/36,36/36,length.out=25)

  weibullSurvival = function(t,inv.b0,alpha,theta_c2,theta_c4,theta_c6){
    theta_c6*theta_c4*theta_c2*exp(-(t/inv.b0)^alpha)
  }
  out=lapply(t6,weibullSurvival,inv.b0=inv.b0.wb[,i],alpha=a.wb[,i],theta_c2[,i],theta_c4[,i],theta_c6[,i])
  out2=do.call(cbind,out)
  out.list=cbind(out.list,out2)
  
  out.list.sum = apply(out.list,2,quantile,probs=c(0.025,.25,.5,.75,.975))
  

  segments(x0=position$easting[i],
           y0=out.list.sum[1,dim(out.list.sum)[2]],
           y1=out.list.sum[5,dim(out.list.sum)[2]])
  segments(x0=position$easting[i],
           y0=out.list.sum[2,dim(out.list.sum)[2]],
           y1=out.list.sum[4,dim(out.list.sum)[2]],lwd=3)
  points(position$easting[i],out.list.sum[3,dim(out.list.sum)[2]],pch=21,bg='white')#,
        # col=ifelse(position$dominant.surface.rock.type[i]=="igneous","black","orange"))
  
}

axis(2, seq(0,1,by=.1), col.ticks = 1)
axis(1, seq(340,375,by=5),
     labels = seq(340,375,by=5), las = 1, 
     col.ticks = 1, cex.axis = 1)
mtext("After 3 years",
      side=3,line=-.5,adj=0,col='black',cex=1)
mtext("Easting (km)",
      side=1,line=2.5,adj=.5,col='black',cex=1)

dev.off()

########################################################
########################################################
# Incorporating loss of viability 
########################################################
########################################################

directoryViability = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBurial/"
modelFittingFilesViability <- paste0(directoryViability,list.files(directoryViability))

samples.rjagsViability <- readRDS(modelFittingFilesViability[[6]])
dataViability <- readRDS(modelFittingFilesViability[[4]])

## POPULATION PLOT

#par(mfrow=c(1,2))
mu_g<-MCMCchains(samples.rjagsViability, params = "mu_g")
mu_g.inv <- apply(mu_g,2,boot::inv.logit)

mu_v<-MCMCchains(samples.rjagsViability, params = "mu_v")
mu_v.inv <- apply(mu_v,2,boot::inv.logit)

nu=exp(mu_g)/(1+exp(mu_g))+(exp(-mu_g+mu_v)/(1+exp(-mu_g)+exp(mu_v)+exp(-mu_g+mu_v)))


mu0_g<-MCMCchains(samples.rjagsViability, params = "mu0_g")
mu0_g.inv <- apply(mu0_g,2,boot::inv.logit)

mu0_v<-MCMCchains(samples.rjagsViability, params = "mu0_v")
mu0_v.inv <- apply(mu0_v,2,boot::inv.logit)

nu0=exp(mu0_g)/(1+exp(mu0_g))+(exp(-mu0_g+mu0_v)/(1+exp(-mu0_g)+exp(mu0_v)+exp(-mu0_g+mu0_v)))

nu0_ratio1=(nu0[,c(1:20)])^(1/3)
nu0_ratio_sum1 = apply(nu0_ratio1,2,quantile,c(.025,.5,.975))

nu0_ratio2=nu0[,c(1:20)]*(nu0[,c(21:40)]/nu0[,c(1:20)])^(1/3)
nu0_ratio2.1=(((nu0[,c(1:20)])^(1+1/3))+((nu0[,c(21:40)])^(16/24)))/2
nu0_ratio2.2 = (2/3)*nu0[,c(1:20)]+(1/3)*nu0[,c(21:40)]
nu0_ratio_sum2 = apply(nu0_ratio2,2,quantile,c(.025,.5,.975))
nu0_ratio_sum2.1 = apply(nu0_ratio2.2,2,quantile,c(.025,.5,.975))

nu0_ratio3=nu0[,c(21:40)]*(nu[,c(101:120)]/nu0[,c(21:40)])^(1/3)
nu0_ratio_sum3 = apply(nu0_ratio3,2,quantile,c(.025,.5,.975))
nu0_ratio3.1=(((nu0[,c(21:40)])^(1+4/12))+((nu[,c(101:120)])^(28/36)))/2
nu0_ratio3.2 = (2/3)*nu0[,c(21:40)]+(1/3)*nu[,c(101:120)]

nu0_ratio_sum3.1 = apply(nu0_ratio3.2,2,quantile,c(.025,.5,.975))


g_popyear=apply(MCMCchains(samples.rjagsViability, params = "mu_g"),2,boot::inv.logit)
v_popyear=apply(MCMCchains(samples.rjagsViability, params = "mu_v"),2,boot::inv.logit)

nu_pop = g_popyear + v_popyear*(1-g_popyear)
nu.sum = apply(nu_pop,2,quantile,c(.025,.5,.975))

g_pop=apply(MCMCchains(samples.rjagsViability, params = "mu0_g"),2,boot::inv.logit)
v_pop=apply(MCMCchains(samples.rjagsViability, params = "mu0_v"),2,boot::inv.logit)

nu0_pop = g_pop + v_pop*(1-g_pop)
nu0.sum = apply(nu0_pop,2,quantile,c(.025,.5,.975))

dev.off()
op <- par('mfrow','oma','mar')

#pdf("~/Dropbox/clarkiaSeedBanks/products/figures/viability-estimates-population.pdf",width=8,height=6)
# 
# par(mfrow = c(4,5),
#     oma = c(5,4,0,0) + 0.1,
#     mar = c(0,0,1,1) + 0.1)
# for(i in 1:20){
#  hist(nu0_ratio1[,i],main='',
#       ylab='',xlab='',yaxt='n');
#   abline(v=1,col='red')
#   text(x=.5,y=.05,siteNames[i])
#   ifelse(i%in%c(16:20),axis(1L, at = c(0,1,2,3)),NA)
#   ifelse(i%in%c(1,6,11,16),axis(2L, at = c(0, .2, .4, .6, .8, 1)),NA)
# }
# 
# for(i in 1:20){
#   hist(nu0_ratio2[,i],main='',
#        ylab='',xlab='',yaxt='n');
#   abline(v=1,col='red')
#   text(x=.5,y=.05,siteNames[i])
#   ifelse(i%in%c(16:20),axis(1L, at = c(0,1,2,3)),NA)
#   ifelse(i%in%c(1,6,11,16),axis(2L, at = c(0, .2, .4, .6, .8, 1)),NA)
# }
# 
# for(i in 1:20){
#   hist(nu0_ratio3[,i],main='',
#        ylab='',xlab='',yaxt='n');
#   abline(v=1,col='red')
#   text(x=.5,y=.05,siteNames[i])
#   ifelse(i%in%c(16:20),axis(1L, at = c(0,1,2,3)),NA)
#   ifelse(i%in%c(1,6,11,16),axis(2L, at = c(0, .2, .4, .6, .8, 1)),NA)
# }

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  plot(NA,NA,type='n',xlim=c(0,3.5),ylim=c(0,1),
       ylab='',xlab='',xaxt='n',yaxt='n')
  
  points(x = 0, y = 1, pch = 21, col = 'lightgray',bg='white')
  
  tmp = nu0.sum[,c(i,i+20)]
  segments(y0=tmp[1,],y1=tmp[3,],x0=c(1,2))
  points(x=c(1,2),y=tmp[2,],pch=19)
  
  tmp2 = as.matrix(nu.sum[,c(i+100)])
  segments(y0=tmp2[1,],y1=tmp2[3,],x0=c(3))
  points(x=c(3),y=tmp2[2,],pch=19)
  
  segments(y0=nu0_ratio_sum1[1,i],y1=nu0_ratio_sum1[3,i],x0=c(1/3),col='black',lty='dotted')
  points(x=1/3,(nu0_ratio_sum1[2,i]),pch=21,bg='white')

  segments(y0=nu0_ratio_sum2[1,i],y1=nu0_ratio_sum2[3,i],x0=c(1+1/3),col='black',lty='dotted')
  points(x=1+1/3,(nu0_ratio_sum2[2,i]),pch=21,bg='white')
  
  segments(y0=nu0_ratio_sum2.1[1,i],y1=nu0_ratio_sum2.1[3,i],x0=c(1+1/3+.1),col='red',lty='dotted')
  points(x=1+1/3+.1,(nu0_ratio_sum2.1[2,i]),pch=21,bg='red')

    segments(y0=nu0_ratio_sum3[1,i],y1=nu0_ratio_sum3[3,i],x0=c(2+1/3),col='black',lty='dotted')
  points(x=2+1/3,(nu0_ratio_sum3[2,i]),pch=21,bg='white')

  segments(y0=nu0_ratio_sum3.1[1,i],y1=nu0_ratio_sum3.1[3,i],x0=c(2+1/3+.1),col='red',lty='dotted')
  points(x=2+1/3+.1,(nu0_ratio_sum3.1[2,i]),pch=21,bg='red')

  text(x=.5,y=.05,siteNames[i])
  ifelse(i%in%c(16:20),axis(1L, at = c(0,1,2,3)),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L, at = c(0, .2, .4, .6, .8, 1)),NA)
}
mtext("Age (years)", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of viability", side = 2, outer = TRUE, line = 2.2)
mtext("Population-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)


par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  plot(NA,NA,type='n',xlim=c(0,3.5),ylim=c(0,1),
       ylab='',xlab='',xaxt='n',yaxt='n')
  
  points(x = 0, y = 1, pch = 21, col = 'lightgray',bg='white')
  
  tmp = nu0.sum[,c(i,i+20)]
  segments(y0=tmp[1,],y1=tmp[3,],x0=c(1,2))
  points(x=c(1,2),y=tmp[2,],pch=19)
  
  tmp2 = as.matrix(nu.sum[,c(i+100)])
  segments(y0=tmp2[1,],y1=tmp2[3,],x0=c(3))
  points(x=c(3),y=tmp2[2,],pch=19)
  
  segments(y0=nu0_ratio_sum1[1,i],y1=nu0_ratio_sum1[3,i],x0=c(1/3),col='black',lty='dotted')
  points(x=1/3,(nu0_ratio_sum1[2,i]),pch=21,bg='white')
  
  for(j in 1:11){
    nu0_ratio.len=nu0[,i]^(seq(0,1,by=.1)[j])
    nu0_ratio_sum.len = quantile(nu0_ratio.len,c(.025,.5,.975))
    segments(y0=nu0_ratio_sum.len[1],y1=nu0_ratio_sum.len[3],x0=c(seq(0,1,by=.1)[j]),col='black',lty='dotted')
   # points(x=seq(0,1,by=.1)[j],(nu0_ratio_sum.len[2]),pch=21,bg='white')
  }
  
  segments(y0=nu0_ratio_sum2[1,i],y1=nu0_ratio_sum2[3,i],x0=c(1+1/3),col='black',lty='dotted')
  points(x=1+1/3,(nu0_ratio_sum2[2,i]),pch=21,bg='white')
  
  # segments(y0=nu0_ratio_sum2.1[1,i],y1=nu0_ratio_sum2.1[3,i],x0=c(1+1/3+.1),col='red',lty='dotted')
  # points(x=1+1/3+.1,(nu0_ratio_sum2.1[2,i]),pch=21,bg='red')
  
  for(j in 1:11){
    nu0_ratio.len=nu0[,i]*(nu0[,i+20]/nu0[,i])^(seq(0,1,by=.1)[j])
    nu0_ratio_sum.len = quantile(nu0_ratio.len,c(.025,.5,.975))
    segments(y0=nu0_ratio_sum.len[1],y1=nu0_ratio_sum.len[3],x0=c(1+seq(0,1,by=.1)[j]),col='black',lty='dotted')
   # points(x=1+seq(0,1,by=.1)[j],(nu0_ratio_sum.len[2]),pch=21,bg='white')
  }
  
  segments(y0=nu0_ratio_sum3[1,i],y1=nu0_ratio_sum3[3,i],x0=c(2+1/3),col='black',lty='dotted')
  points(x=2+1/3,(nu0_ratio_sum3[2,i]),pch=21,bg='white')
  
  for(j in 1:11){
    nu0_ratio.len=nu0[,i+20]*(nu[,i+100]/nu0[,i+20])^(seq(0,1,by=.1)[j])
    nu0_ratio_sum.len = quantile(nu0_ratio.len,c(.025,.5,.975))
    segments(y0=nu0_ratio_sum.len[1],y1=nu0_ratio_sum.len[3],x0=c(2+seq(0,1,by=.1)[j]),col='black',lty='dotted')
   # points(x=2+seq(0,1,by=.1)[j],(nu0_ratio_sum.len[2]),pch=21,bg='white')
  }
  
  # segments(y0=nu0_ratio_sum3.1[1,i],y1=nu0_ratio_sum3.1[3,i],x0=c(2+1/3+.1),col='red',lty='dotted')
  # points(x=2+1/3+.1,(nu0_ratio_sum3.1[2,i]),pch=21,bg='red')
  
  text(x=.5,y=.05,siteNames[i])
  ifelse(i%in%c(16:20),axis(1L, at = c(0,1,2,3)),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L, at = c(0, .2, .4, .6, .8, 1)),NA)
}
mtext("Age (years)", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of viability", side = 2, outer = TRUE, line = 2.2)
mtext("Population-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)



viability.df = cbind(rep(1,20),
                     nu0_ratio_sum1[2,],
                     nu0.sum[2,1:20],
                     nu0_ratio_sum2[2,],
                     nu0.sum[2,21:40],
                     nu0_ratio_sum3[2,],
                     nu.sum[2,101:120])


par(op)
plot(x = c(0,4,12,16,24,28,36),
     y = viability.df[1,], type = 'n',
     ylim=c(0,1),
     ylab='',xlab='',xaxt='n',yaxt='n')
for(i in 1:20){
  lines(x = c(0,4,12,16,24,28,36),
         y = viability.df[i,],
        type='b',
         col = c('lightgray',1,1,1,1,1,1),
         pch = c(21,21,19,21,19,21,19),
        cex = .75)
}
axis(1L, at = c(0,4,8,12,16,20,24,28,32,36))
axis(2L, at = c(0, .2, .4, .6, .8, 1))
mtext("Time (months)", side = 1, line = 2.2)
mtext("Probability of viability", side = 2, line = 2.2)

rect(c(0,4,12,16,24,28),0,c(4,12,16,24,28,36),-.025,
     col=c('gray50','gray90'),lwd=0)
text(c(0,4,12,16,24,28,36),0.025,
     c("O","J"),cex=.75)
dev.off()


########################################################
# Calculate survival function for persistence+viability
########################################################

## get germination probabilities

mu0_g=MCMCchains(mcmcSamples,params="mu0_g")
mu_g=MCMCchains(mcmcSamples,params="mu_g")

gamma1 = boot::inv.logit(mu0_g[,1:20])
gamma1.sum=apply(gamma1,2,quantile,probs=c(0.025,.25,.5,.75,.975))

gamma2 = boot::inv.logit(mu0_g[,21:40])
gamma2.sum=apply(gamma2,2,quantile,probs=c(0.025,.25,.5,.75,.975))

# note here this is the population-year level 
# there is NO population-level parameter for age 3 germination
gamma3 = boot::inv.logit(mu_g[,101:120])
gamma3.sum=apply(gamma3,2,quantile,probs=c(0.025,.25,.5,.75,.975))


theta_c2 = 1- gamma1
theta_c4 = (1-gamma2)
theta_c6 = (1-gamma3)

# discretize compound persistence function
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

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
plot(NA,NA,xlim=c(0,36),ylim=c(0,1),type='l')
surv.fun=apply(cbind(th_0[,i],th_1[,i],th_2[,i],th_3[,i],th_4[,i],th_5[,i],th_6[,i],th_7[,i],th_8[,i],th_9[,i]),2,median)
lines(t[x]*36,surv.fun,lty='solid',col='lightgray')
points(t[x]*36,surv.fun,pch=16,col='black')
}

for(i in 1:20){
  plot(NA,NA,xlim=c(0,36),ylim=c(0,1),ylab='',xlab='',xaxt='n',yaxt='n')
  surv.fun=apply(cbind(th_0[,i],th_1[,i],th_2[,i],th_3[,i],th_4[,i],th_5[,i],th_6[,i],th_7[,i],th_8[,i],th_9[,i]),2,quantile,c(0.025,.5,.975))
  #y.vals=rep(surv.fun,each=2)
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
  lines(t[x[c(2,3,3)]]*36,y.vals[2,c(2,2,3)],lty='dotted',col='red')
  lines(t[x[c(3,4,4)]]*36,y.vals[2,c(3,3,4)])
  lines(t[x[c(4,5,5)]]*36,y.vals[2,c(4,4,5)])
  lines(t[x[c(5,6,6)]]*36,y.vals[2,c(5,5,6)],lty='dotted',col='red')
  lines(t[x[c(6,7,7)]]*36,y.vals[2,c(6,6,7)])
  lines(t[x[c(7,8,8)]]*36,y.vals[2,c(7,7,8)])
  lines(t[x[c(8,9,9)]]*36,y.vals[2,c(8,8,9)],lty='dotted',col='red')
  lines(t[x[c(9,10,10)]]*36,y.vals[2,c(9,9,10)])
  
  text(x=3,y=.05,siteNames[i])
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}


## get viability probabilities

# get population estimates for year 1-2
mu0_g<-MCMCchains(samples.rjagsViability, params = "mu0_g")
mu0_v<-MCMCchains(samples.rjagsViability, params = "mu0_v")
# calculate total viability; calculation on latent scale
nu0=exp(mu0_g)/(1+exp(mu0_g))+(exp(-mu0_g+mu0_v)/(1+exp(-mu0_g)+exp(mu0_v)+exp(-mu0_g+mu0_v)))

# get population*year estimates for year 3
mu_g<-MCMCchains(samples.rjagsViability, params = "mu_g")
mu_v<-MCMCchains(samples.rjagsViability, params = "mu_v")
# calculate total viability; calculation on latent scale
nu=exp(mu_g)/(1+exp(mu_g))+(exp(-mu_g+mu_v)/(1+exp(-mu_g)+exp(mu_v)+exp(-mu_g+mu_v)))

# interpolate for January estimates
nu0_ratio1=(nu0[,c(1:20)])^(1/3)
#nu0_ratio2=nu0[,c(1:20)]*(nu0[,c(21:40)]/nu0[,c(1:20)])^(1/3)
#nu0_ratio3=nu0[,c(21:40)]*(nu[,c(101:120)]/nu0[,c(21:40)])^(1/3)
nu0_ratio2=nu0_ratio2.1
nu0_ratio3=nu0_ratio3.1

## to check, do calculations on probability scale as well
# g_popyear=apply(MCMCchains(samples.rjagsViability, params = "mu_g"),2,boot::inv.logit)
# v_popyear=apply(MCMCchains(samples.rjagsViability, params = "mu_v"),2,boot::inv.logit)
# nu_pop = g_popyear + v_popyear*(1-g_popyear)
# sum(nu_pop-nu)
# 
# g_pop=apply(MCMCchains(samples.rjagsViability, params = "mu0_g"),2,boot::inv.logit)
# v_pop=apply(MCMCchains(samples.rjagsViability, params = "mu0_v"),2,boot::inv.logit)
# nu0_pop = g_pop + v_pop*(1-g_pop)
# sum(nu0_pop-nu0)

nu_01.inter = nu0_ratio1
nu_1 = nu0[,1:20]
nu_12.inter = nu0_ratio2
nu_2 = nu0[,21:40]
nu_23.inter = nu0_ratio3
nu_3 = nu[,101:120]

## construct viability + persistence survival function

## conditional germination 

gamma1.v = gamma1/(1-(1-nu_01.inter)*(1-gamma1))
gamma2.v = gamma2/(1-(1-nu_12.inter)*(1-gamma2))
gamma3.v = gamma3/(1-(1-nu_23.inter)*(1-gamma3))

f.wb = function(t,inv.b0,alpha){
  exp(-(t/inv.b0)^alpha)
}

Oct_0.v = f.wb(t=t[x[1]],inv.b0=inv.b0.wb,alpha=a.wb)

Janpre_1.v = f.wb(t=t[x[2]],inv.b0=inv.b0.wb,alpha=a.wb)*(gamma1+(1-gamma1)*nu_01.inter)
Jangerm_1.v = gamma1.v
Janpost_1.v = f.wb(t=t[x[3]],inv.b0=inv.b0.wb,alpha=a.wb)*(1-gamma1)*nu_01.inter
Oct_1.v = (f.wb(t=t[x[4]],inv.b0=inv.b0.wb,alpha=a.wb)*nu_1*(1-gamma1))

Janpre_2.v =( f.wb(t=t[x[5]],inv.b0=inv.b0.wb,alpha=a.wb)*(1-gamma1)*(gamma2+(1-gamma2)*nu_12.inter))
Jangerm_2.v = gamma2.v
Janpost_2.v = f.wb(t=t[x[6]],inv.b0=inv.b0.wb,alpha=a.wb)*(1-gamma1)*(1-gamma2)*nu_12.inter
Oct_2.v = (f.wb(t=t[x[7]],inv.b0=inv.b0.wb,alpha=a.wb)*(1-gamma1)*(1-gamma2)*nu_2)

Janpre_3.v = ( f.wb(t=t[x[8]],inv.b0=inv.b0.wb,alpha=a.wb)*(1-gamma1)*(1-gamma2)*(gamma3+(1-gamma3)*nu_23.inter))
Jangerm_3.v = gamma3.v
Janpost_3.v = f.wb(t=t[x[9]],inv.b0=inv.b0.wb,alpha=a.wb)*(1-gamma1)*(1-gamma2)*(1-gamma3)*nu_23.inter
Oct_3.v =(f.wb(t=t[x[10]],inv.b0=inv.b0.wb,alpha=a.wb)*(1-gamma1)*(1-gamma2)*(1-gamma3)*nu_3)
# 
# for(i in 1:20){hist(Janpre_1.v[,i],main="")}
# for(i in 1:20){hist(Oct_1.v[,i],main="")}
# for(i in 1:20){hist(Janpre_2.v[,i],main="")}
# for(i in 1:20){hist(Oct_2.v[,i],main="")}
# for(i in 1:20){hist(Janpre_3.v[,i],main="")}
# for(i in 1:20){hist(Oct_3.v[,i],main="")}
# 
# for(i in 1:20){hist(Jangerm_1.v[,i],main="")}
# for(i in 1:20){hist(Jangerm_2.v[,i],main="")}
# for(i in 1:20){hist(Jangerm_3.v[,i],main="")}

phi_0 = Oct_0.v
phi_1 = Janpre_1.v
phi_2 = Janpost_1.v
phi_3 = Oct_1.v
phi_4 = Janpre_2.v
phi_5 = Janpost_2.v
phi_6 = Oct_2.v
phi_7 = Janpre_3.v
phi_8 = Janpost_3.v
phi_9 = Oct_3.v

# if calculating survival by incorporating viability, end up with survival > 1
# by including viability, the trajectory is no longer a proper survival function
#for(i in 1:20){hist(phi_0[,i]-phi_1[,i],main='');abline(v=0,col='red')}
#for(i in 1:20){hist(phi_1[,i]-phi_2[,i],main='');abline(v=0,col='red')}
#for(i in 1:20){hist(phi_2[,i]-phi_3[,i],main='');abline(v=c(0,1),col='red')}
for(i in 1:20){hist(phi_3[,i]-phi_4[,i],main='');abline(v=c(0,1),col='red')}
#for(i in 1:20){hist(phi_4[,i]-phi_5[,i],main='');abline(v=0,col='red')}
for(i in 1:20){hist(phi_5[,i]-phi_6[,i],main='');abline(v=0,col='red')}
for(i in 1:20){hist(phi_6[,i]-phi_7[,i],main='');abline(v=0,col='red')}
#for(i in 1:20){hist(phi_7[,i]-phi_8[,i],main='');abline(v=0,col='red')}
for(i in 1:20){hist(phi_8[,i]-phi_9[,i],main='');abline(v=0,col='red')}



# phi_0 = th_0
# phi_1 = th_1*(gamma1+(1-gamma1)*nu_01.inter)
# phi_2 = th_2*nu_01.inter
# phi_3 = th_3*nu_1
# phi_4 = th_4*(gamma2+(1-gamma2)*nu_12.inter)
# phi_5 = th_5*nu_12.inter
# phi_6 = th_6*nu_2
# phi_7 = th_7*(gamma3+(1-gamma3)*nu_23.inter)
# phi_8 = th_8*nu_23.inter
# phi_9 = th_9*nu_3


# plot
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  plot(NA,NA,xlim=c(0,36),ylim=c(0,1),ylab='',xlab='',xaxt='n',yaxt='n')
  
  surv.fun_persistence=apply(cbind(th_0[,i],th_1[,i],th_2[,i],th_3[,i],th_4[,i],th_5[,i],th_6[,i],th_7[,i],th_8[,i],th_9[,i]),2,median)
  lines((t[x]*36)[1:10],surv.fun_persistence[1:10],lty='solid',col='lightgray')
  points((t[x]*36)[1:10],surv.fun_persistence[1:10],pch=16,col='black')

  surv.fun_persistence.viability=apply(cbind(phi_0[,i],phi_1[,i],phi_2[,i],phi_3[,i],phi_4[,i],phi_5[,i],phi_6[,i],phi_7[,i],phi_8[,i],phi_9[,i]),2,median)
  lines((t[x]*36)[1:10],surv.fun_persistence.viability[1:10],lty='solid',col=
        ifelse(all(diff(surv.fun_persistence.viability[1:10]) <= 0),'orange','blue'))
  points((t[x]*36)[1:10],surv.fun_persistence.viability[1:10],pch=16,
         col=ifelse(all(diff(surv.fun_persistence.viability[1:10]) <= 0),'orange','blue'))
  text(x=3,y=.05,siteNames[i])
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Time (months)", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of persistence and viability", side = 2, outer = TRUE, line = 2.2)
mtext("Population-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  plot(NA,NA,xlim=c(0,36),ylim=c(0,1),ylab='',xlab='',xaxt='n',yaxt='n')
 # surv.fun=apply(cbind(th_0[,i],th_1[,i],th_2[,i],th_3[,i],th_4[,i],th_5[,i],th_6[,i],th_7[,i],th_8[,i],th_9[,i]),2,median)
#  y.vals=rep(surv.fun,each=2)
  
  surv.fun_persistence.viability=apply(cbind(phi_0[,i],phi_1[,i],phi_2[,i],phi_3[,i],phi_4[,i],phi_5[,i],phi_6[,i],phi_7[,i],phi_8[,i],phi_9[,i]),2,quantile,c(0.025,.5,.975))
  y.vals=surv.fun_persistence.viability
  
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
  lines(t[x[c(2,3,3)]]*36,y.vals[2,c(2,2,3)],lty='dotted',col='red')
  lines(t[x[c(3,4,4)]]*36,y.vals[2,c(3,3,4)])
  lines(t[x[c(4,5,5)]]*36,y.vals[2,c(4,4,5)])
  lines(t[x[c(5,6,6)]]*36,y.vals[2,c(5,5,6)],lty='dotted',col='red')
  lines(t[x[c(6,7,7)]]*36,y.vals[2,c(6,6,7)])
  lines(t[x[c(7,8,8)]]*36,y.vals[2,c(7,7,8)])
  lines(t[x[c(8,9,9)]]*36,y.vals[2,c(8,8,9)],lty='dotted',col='red')
  lines(t[x[c(9,10,10)]]*36,y.vals[2,c(9,9,10)])
  # 
  # lines(t[x[c(1,2,2)]]*36,y.vals[c(1,2,3)])
  # lines(t[x[c(2,3,3)]]*36,y.vals[c(3:5)],lty='dotted')
  # lines(t[x[c(2,3,3)]]*36,y.vals.v[c(3:5)],lty='dotted',col=ifelse(all(diff(y.vals.v[c(3:5)]) <= 0),'orange','blue'))
  # 
  # lines(t[x[c(3,4,4)]]*36,y.vals[c(5:7)])
  # lines(t[x[c(4,5,5)]]*36,y.vals[c(7:9)])
  # lines(t[x[c(5,6,6)]]*36,y.vals[c(9:11)],lty='dotted')
  # lines(t[x[c(6,7,7)]]*36,y.vals[c(11:13)])
  # lines(t[x[c(7,8,8)]]*36,y.vals[c(13:15)])
  # lines(t[x[c(8,9,9)]]*36,y.vals[c(15:17)],lty='dotted')
  # lines(t[x[c(9,10,10)]]*36,y.vals[c(17:19)])
  # 
  # 
  # lines(t[x[c(1,2,2)]]*36,y.vals.v[c(1,2,3)],col=ifelse(all(diff(y.vals.v[c(1,2,3)]) <= 0),'orange','blue'))
  # lines(t[x[c(2,3,3)]]*36,y.vals.v[c(3:5)],lty='dotted',col=ifelse(all(diff(y.vals.v[c(3:5)]) <= 0),'orange','blue'))
  # lines(t[x[c(3,4,4)]]*36,y.vals.v[c(5:7)],col=ifelse(all(diff(y.vals.v[c(5:7)]) <= 0),'orange','blue'))
  # lines(t[x[c(4,5,5)]]*36,y.vals.v[c(7:9)],col=ifelse(all(diff(y.vals.v[c(7:9)]) <= 0),'orange','blue'))
  # lines(t[x[c(5,6,6)]]*36,y.vals.v[c(9:11)],lty='dotted',col=ifelse(all(diff(y.vals.v[c(9:11)]) <= 0),'orange','blue'))
  # lines(t[x[c(6,7,7)]]*36,y.vals.v[c(11:13)],col=ifelse(all(diff(y.vals.v[c(11:13)]) <= 0),'orange','blue'))
  # lines(t[x[c(7,8,8)]]*36,y.vals.v[c(13:15)],col=ifelse(all(diff(y.vals.v[c(13:15)]) <= 0),'orange','blue'))
  # lines(t[x[c(8,9,9)]]*36,y.vals.v[c(15:17)],lty='dotted',col=ifelse(all(diff(y.vals.v[c(15:17)]) <= 0),'orange','blue'))
  # lines(t[x[c(9,10,10)]]*36,y.vals.v[c(17:19)],col=ifelse(all(diff(y.vals.v[c(17:19)]) <= 0),'orange','blue'))
  # 
  text(x=3,y=.05,siteNames[i])
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}


par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  plot(NA,NA,xlim=c(0,10),ylim=c(-1,1),ylab='',xlab='',xaxt='n',yaxt='n')
  abline(h=0)
  surv.fun_persistence=apply(cbind(th_0[,i],th_1[,i],th_2[,i],th_3[,i],th_4[,i],th_5[,i],th_6[,i],th_7[,i],th_8[,i],th_9[,i]),2,median)
  # lines(t[x]*36,surv.fun_persistence,lty='solid',col='lightgray')
  # points(t[x]*36,surv.fun_persistence,pch=16,col='black')
  
  surv.fun_persistence.viability=apply(cbind(phi_0[,i],phi_1[,i],phi_2[,i],phi_3[,i],phi_4[,i],phi_5[,i],phi_6[,i],phi_7[,i],phi_8[,i],phi_9[,i]),2,median)
  # lines(t[x]*36,surv.fun_persistence.viability,lty='solid',col=
  #         ifelse(all(diff(surv.fun_persistence.viability) <= 0),'orange','blue'))
  # points(t[x]*36,surv.fun_persistence.viability,pch=16,
  #        col=ifelse(all(diff(surv.fun_persistence.viability) <= 0),'orange','blue'))
  
  pmf = surv.fun_persistence.viability[1:9]-surv.fun_persistence.viability[2:10]
  h = pmf/surv.fun_persistence.viability[1:9]
  points(1:9,h,col=ifelse(h<0,"orange","blue"))
  
   text(x=3,y=-.9,siteNames[i])
    ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

# if calculating survival by incorporating viability, end up with survival > 1
# by including viability, the trajectory is no longer a proper survival function
#for(i in 1:20){hist(phi_0[,i]-phi_1[,i],main='');abline(v=0,col='red')}
#for(i in 1:20){hist(phi_1[,i]-phi_2[,i],main='');abline(v=0,col='red')}
#for(i in 1:20){hist(phi_2[,i]-phi_3[,i],main='');abline(v=0,col='red')}
for(i in 1:20){hist(phi_3[,i]-phi_4[,i],main='');abline(v=0,col='red')}
#for(i in 1:20){hist(phi_4[,i]-phi_5[,i],main='');abline(v=0,col='red')}
for(i in 1:20){hist(phi_5[,i]-phi_6[,i],main='');abline(v=0,col='red')}
for(i in 1:20){hist(phi_6[,i]-phi_7[,i],main='');abline(v=0,col='red')}
#for(i in 1:20){hist(phi_7[,i]-phi_8[,i],main='');abline(v=0,col='red')}
for(i in 1:20){hist(phi_8[,i]-phi_9[,i],main='');abline(v=0,col='red')}


# calculate structured model parameters
g1=gamma1/phi_1
g2=gamma2/phi_4
g3=gamma3/phi_7

g1.2 = gamma1/(1-(1-nu_01.inter)*(1-gamma1))
g2.2
g3.2

s1=phi_1
s2=phi_3/phi_2
s3=phi_4/phi_3
s4=phi_6/phi_5
s5=phi_7/phi_6
s6=phi_9/phi_8

s1.og=th_1
s2.og=th_3/th_2
s3.og=th_4/th_3
s4.og=th_6/th_5
s5.og=th_7/th_6
s6.og=th_9/th_8

dev.off()
plot(apply(gamma1,2,median),apply(g1,2,median),xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)
plot(apply(gamma2,2,median),apply(g2,2,median),xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)
plot(apply(gamma3,2,median),apply(g3,2,median),xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)

plot(apply(gamma1,2,median),apply(g1,2,median),xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)
plot(apply(g1.2,2,median),apply(g1,2,median),xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)
plot(apply(gamma1,2,median),apply(g1.2,2,median),xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)

par(mfrow=c(1,1))
plot(apply(s1.og,2,median),apply(s1,2,median),xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)
plot(apply(s2.og,2,median),apply(s2,2,median),xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)
plot(apply(s3.og,2,median),apply(s3,2,median),xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)
plot(apply(s4.og,2,median),apply(s4,2,median),xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)
plot(apply(s5.og,2,median),apply(s5,2,median),xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)
plot(apply(s6.og,2,median),apply(s6,2,median),xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)

hist(s5[,2])

