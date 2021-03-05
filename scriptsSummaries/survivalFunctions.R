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
# Incorporating loss of viability 
# -------------------------------------------------------------------


directoryViability = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBurial/"
modelFittingFilesViability <- paste0(directoryViability,list.files(directoryViability))

samples.rjagsViability <- readRDS(modelFittingFilesViability[[6]])
dataViability <- readRDS(modelFittingFilesViability[[4]])
# 
# # -------------------------------------------------------------------
# # Viability
# # -------------------------------------------------------------------
# 
# #par(mfrow=c(1,2))
# mu_g<-MCMCchains(samples.rjagsViability, params = "mu_g")
# mu_g.inv <- apply(mu_g,2,boot::inv.logit)
# 
# mu_v<-MCMCchains(samples.rjagsViability, params = "mu_v")
# mu_v.inv <- apply(mu_v,2,boot::inv.logit)
# 
# nu=exp(mu_g)/(1+exp(mu_g))+(exp(-mu_g+mu_v)/(1+exp(-mu_g)+exp(mu_v)+exp(-mu_g+mu_v)))
# 
# 
# mu0_g<-MCMCchains(samples.rjagsViability, params = "mu0_g")
# mu0_g.inv <- apply(mu0_g,2,boot::inv.logit)
# 
# mu0_v<-MCMCchains(samples.rjagsViability, params = "mu0_v")
# mu0_v.inv <- apply(mu0_v,2,boot::inv.logit)
# 
# nu0=exp(mu0_g)/(1+exp(mu0_g))+(exp(-mu0_g+mu0_v)/(1+exp(-mu0_g)+exp(mu0_v)+exp(-mu0_g+mu0_v)))
# 
# nu0_ratio1=(nu0[,c(1:20)])^(1/3)
# nu0_ratio_sum1 = apply(nu0_ratio1,2,quantile,c(.025,.5,.975))
# 
# nu0_ratio2=nu0[,c(1:20)]*(nu0[,c(21:40)]/nu0[,c(1:20)])^(1/3)
# nu0_ratio2.1=(((nu0[,c(1:20)])^(1+1/3))+((nu0[,c(21:40)])^(16/24)))/2
# nu0_ratio2.2 = (2/3)*nu0[,c(1:20)]+(1/3)*nu0[,c(21:40)]
# nu0_ratio_sum2 = apply(nu0_ratio2,2,quantile,c(.025,.5,.975))
# nu0_ratio_sum2.1 = apply(nu0_ratio2.2,2,quantile,c(.025,.5,.975))
# 
# nu0_ratio3=nu0[,c(21:40)]*(nu[,c(101:120)]/nu0[,c(21:40)])^(1/3)
# nu0_ratio_sum3 = apply(nu0_ratio3,2,quantile,c(.025,.5,.975))
# nu0_ratio3.1=(((nu0[,c(21:40)])^(1+4/12))+((nu[,c(101:120)])^(28/36)))/2
# nu0_ratio3.2 = (2/3)*nu0[,c(21:40)]+(1/3)*nu[,c(101:120)]
# 
# nu0_ratio_sum3.1 = apply(nu0_ratio3.2,2,quantile,c(.025,.5,.975))
# 
# 
# g_popyear=apply(MCMCchains(samples.rjagsViability, params = "mu_g"),2,boot::inv.logit)
# v_popyear=apply(MCMCchains(samples.rjagsViability, params = "mu_v"),2,boot::inv.logit)
# 
# nu_pop = g_popyear + v_popyear*(1-g_popyear)
# nu.sum = apply(nu_pop,2,quantile,c(.025,.5,.975))
# 
# g_pop=apply(MCMCchains(samples.rjagsViability, params = "mu0_g"),2,boot::inv.logit)
# v_pop=apply(MCMCchains(samples.rjagsViability, params = "mu0_v"),2,boot::inv.logit)
# 
# nu0_pop = g_pop + v_pop*(1-g_pop)
# nu0.sum = apply(nu0_pop,2,quantile,c(.025,.5,.975))
# 
# dev.off()
# op <- par('mfrow','oma','mar')
# 
# #pdf("~/Dropbox/clarkiaSeedBanks/products/figures/viability-estimates-population.pdf",width=8,height=6)
# # 
# # par(mfrow = c(4,5),
# #     oma = c(5,4,0,0) + 0.1,
# #     mar = c(0,0,1,1) + 0.1)
# # for(i in 1:20){
# #  hist(nu0_ratio1[,i],main='',
# #       ylab='',xlab='',yaxt='n');
# #   abline(v=1,col='red')
# #   text(x=.5,y=.05,siteNames[i])
# #   ifelse(i%in%c(16:20),axis(1L, at = c(0,1,2,3)),NA)
# #   ifelse(i%in%c(1,6,11,16),axis(2L, at = c(0, .2, .4, .6, .8, 1)),NA)
# # }
# # 
# # for(i in 1:20){
# #   hist(nu0_ratio2[,i],main='',
# #        ylab='',xlab='',yaxt='n');
# #   abline(v=1,col='red')
# #   text(x=.5,y=.05,siteNames[i])
# #   ifelse(i%in%c(16:20),axis(1L, at = c(0,1,2,3)),NA)
# #   ifelse(i%in%c(1,6,11,16),axis(2L, at = c(0, .2, .4, .6, .8, 1)),NA)
# # }
# # 
# # for(i in 1:20){
# #   hist(nu0_ratio3[,i],main='',
# #        ylab='',xlab='',yaxt='n');
# #   abline(v=1,col='red')
# #   text(x=.5,y=.05,siteNames[i])
# #   ifelse(i%in%c(16:20),axis(1L, at = c(0,1,2,3)),NA)
# #   ifelse(i%in%c(1,6,11,16),axis(2L, at = c(0, .2, .4, .6, .8, 1)),NA)
# # }
# 
# par(mfrow = c(4,5),
#     oma = c(5,4,0,0) + 0.1,
#     mar = c(0,0,1,1) + 0.1)
# for(i in 1:20){
#   plot(NA,NA,type='n',xlim=c(0,3.5),ylim=c(0,1),
#        ylab='',xlab='',xaxt='n',yaxt='n')
#   
#   points(x = 0, y = 1, pch = 21, col = 'lightgray',bg='white')
#   
#   tmp = nu0.sum[,c(i,i+20)]
#   segments(y0=tmp[1,],y1=tmp[3,],x0=c(1,2))
#   points(x=c(1,2),y=tmp[2,],pch=19)
#   
#   tmp2 = as.matrix(nu.sum[,c(i+100)])
#   segments(y0=tmp2[1,],y1=tmp2[3,],x0=c(3))
#   points(x=c(3),y=tmp2[2,],pch=19)
#   
#   segments(y0=nu0_ratio_sum1[1,i],y1=nu0_ratio_sum1[3,i],x0=c(1/3),col='black',lty='dotted')
#   points(x=1/3,(nu0_ratio_sum1[2,i]),pch=21,bg='white')
#   
#   segments(y0=nu0_ratio_sum2[1,i],y1=nu0_ratio_sum2[3,i],x0=c(1+1/3),col='black',lty='dotted')
#   points(x=1+1/3,(nu0_ratio_sum2[2,i]),pch=21,bg='white')
#   
#   segments(y0=nu0_ratio_sum2.1[1,i],y1=nu0_ratio_sum2.1[3,i],x0=c(1+1/3+.1),col='red',lty='dotted')
#   points(x=1+1/3+.1,(nu0_ratio_sum2.1[2,i]),pch=21,bg='red')
#   
#   segments(y0=nu0_ratio_sum3[1,i],y1=nu0_ratio_sum3[3,i],x0=c(2+1/3),col='black',lty='dotted')
#   points(x=2+1/3,(nu0_ratio_sum3[2,i]),pch=21,bg='white')
#   
#   segments(y0=nu0_ratio_sum3.1[1,i],y1=nu0_ratio_sum3.1[3,i],x0=c(2+1/3+.1),col='red',lty='dotted')
#   points(x=2+1/3+.1,(nu0_ratio_sum3.1[2,i]),pch=21,bg='red')
#   
#   text(x=.5,y=.05,siteNames[i])
#   ifelse(i%in%c(16:20),axis(1L, at = c(0,1,2,3)),NA)
#   ifelse(i%in%c(1,6,11,16),axis(2L, at = c(0, .2, .4, .6, .8, 1)),NA)
# }
# mtext("Age (years)", side = 1, outer = TRUE, line = 2.2)
# mtext("Probability of viability", side = 2, outer = TRUE, line = 2.2)
# mtext("Population-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)
# 
# 
# par(mfrow = c(4,5),
#     oma = c(5,4,0,0) + 0.1,
#     mar = c(0,0,1,1) + 0.1)
# for(i in 1:20){
#   plot(NA,NA,type='n',xlim=c(0,3.5),ylim=c(0,1),
#        ylab='',xlab='',xaxt='n',yaxt='n')
#   
#   points(x = 0, y = 1, pch = 21, col = 'lightgray',bg='white')
#   
#   tmp = nu0.sum[,c(i,i+20)]
#   segments(y0=tmp[1,],y1=tmp[3,],x0=c(1,2))
#   points(x=c(1,2),y=tmp[2,],pch=19)
#   
#   tmp2 = as.matrix(nu.sum[,c(i+100)])
#   segments(y0=tmp2[1,],y1=tmp2[3,],x0=c(3))
#   points(x=c(3),y=tmp2[2,],pch=19)
#   
#   segments(y0=nu0_ratio_sum1[1,i],y1=nu0_ratio_sum1[3,i],x0=c(1/3),col='black',lty='dotted')
#   points(x=1/3,(nu0_ratio_sum1[2,i]),pch=21,bg='white')
#   
#   for(j in 1:11){
#     nu0_ratio.len=nu0[,i]^(seq(0,1,by=.1)[j])
#     nu0_ratio_sum.len = quantile(nu0_ratio.len,c(.025,.5,.975))
#     segments(y0=nu0_ratio_sum.len[1],y1=nu0_ratio_sum.len[3],x0=c(seq(0,1,by=.1)[j]),col='black',lty='dotted')
#     # points(x=seq(0,1,by=.1)[j],(nu0_ratio_sum.len[2]),pch=21,bg='white')
#   }
#   
#   segments(y0=nu0_ratio_sum2[1,i],y1=nu0_ratio_sum2[3,i],x0=c(1+1/3),col='black',lty='dotted')
#   points(x=1+1/3,(nu0_ratio_sum2[2,i]),pch=21,bg='white')
#   
#   # segments(y0=nu0_ratio_sum2.1[1,i],y1=nu0_ratio_sum2.1[3,i],x0=c(1+1/3+.1),col='red',lty='dotted')
#   # points(x=1+1/3+.1,(nu0_ratio_sum2.1[2,i]),pch=21,bg='red')
#   
#   for(j in 1:11){
#     nu0_ratio.len=nu0[,i]*(nu0[,i+20]/nu0[,i])^(seq(0,1,by=.1)[j])
#     nu0_ratio_sum.len = quantile(nu0_ratio.len,c(.025,.5,.975))
#     segments(y0=nu0_ratio_sum.len[1],y1=nu0_ratio_sum.len[3],x0=c(1+seq(0,1,by=.1)[j]),col='black',lty='dotted')
#     # points(x=1+seq(0,1,by=.1)[j],(nu0_ratio_sum.len[2]),pch=21,bg='white')
#   }
#   
#   segments(y0=nu0_ratio_sum3[1,i],y1=nu0_ratio_sum3[3,i],x0=c(2+1/3),col='black',lty='dotted')
#   points(x=2+1/3,(nu0_ratio_sum3[2,i]),pch=21,bg='white')
#   
#   for(j in 1:11){
#     nu0_ratio.len=nu0[,i+20]*(nu[,i+100]/nu0[,i+20])^(seq(0,1,by=.1)[j])
#     nu0_ratio_sum.len = quantile(nu0_ratio.len,c(.025,.5,.975))
#     segments(y0=nu0_ratio_sum.len[1],y1=nu0_ratio_sum.len[3],x0=c(2+seq(0,1,by=.1)[j]),col='black',lty='dotted')
#     # points(x=2+seq(0,1,by=.1)[j],(nu0_ratio_sum.len[2]),pch=21,bg='white')
#   }
#   
#   # segments(y0=nu0_ratio_sum3.1[1,i],y1=nu0_ratio_sum3.1[3,i],x0=c(2+1/3+.1),col='red',lty='dotted')
#   # points(x=2+1/3+.1,(nu0_ratio_sum3.1[2,i]),pch=21,bg='red')
#   
#   text(x=.5,y=.05,siteNames[i])
#   ifelse(i%in%c(16:20),axis(1L, at = c(0,1,2,3)),NA)
#   ifelse(i%in%c(1,6,11,16),axis(2L, at = c(0, .2, .4, .6, .8, 1)),NA)
# }
# mtext("Age (years)", side = 1, outer = TRUE, line = 2.2)
# mtext("Probability of viability", side = 2, outer = TRUE, line = 2.2)
# mtext("Population-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)
# 
# 
# 
# viability.df = cbind(rep(1,20),
#                      nu0_ratio_sum1[2,],
#                      nu0.sum[2,1:20],
#                      nu0_ratio_sum2[2,],
#                      nu0.sum[2,21:40],
#                      nu0_ratio_sum3[2,],
#                      nu.sum[2,101:120])
# 
# 
# par(op)
# plot(x = c(0,4,12,16,24,28,36),
#      y = viability.df[1,], type = 'n',
#      ylim=c(0,1),
#      ylab='',xlab='',xaxt='n',yaxt='n')
# for(i in 1:20){
#   lines(x = c(0,4,12,16,24,28,36),
#         y = viability.df[i,],
#         type='b',
#         col = c(ifelse(i %in% c(1,18,19),'orange','lightgray'),1,1,1,1,1,1),
#         pch = c(21,21,19,21,19,21,19),
#         cex = .75)
# }
# axis(1L, at = c(0,4,8,12,16,20,24,28,32,36))
# axis(2L, at = c(0, .2, .4, .6, .8, 1))
# mtext("Time (months)", side = 1, line = 2.2)
# mtext("Probability of viability", side = 2, line = 2.2)
# 
# rect(c(0,4,12,16,24,28),0,c(4,12,16,24,28,36),-.025,
#      col=c('gray50','gray90'),lwd=0)
# text(c(0,4,12,16,24,28,36),0.025,
#      c("O","J"),cex=.75)
# dev.off()

# -------------------------------------------------------------------
# Germination
# -------------------------------------------------------------------

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

# -------------------------------------------------------------------
# Discrete survival component
# -------------------------------------------------------------------

theta_c2 = 1- gamma1
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

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  plot(NA,NA,xlim=c(0,36),ylim=c(0,1),type='l')
  surv.fun=apply(cbind(th_0[,i],th_1[,i],th_2[,i],th_3[,i],th_4[,i],th_5[,i],th_6[,i],th_7[,i],th_8[,i],th_9[,i]),2,median)
  lines(t[x]*36,surv.fun,lty='solid',col='lightgray')
  #points(t[x]*36,surv.fun,pch=16,col='black')
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

# -------------------------------------------------------------------
# Viability
# -------------------------------------------------------------------

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
nu0_ratio2=nu0[,c(1:20)]*(nu0[,c(21:40)]/nu0[,c(1:20)])^(1/3)
nu0_ratio3=nu0[,c(21:40)]*(nu[,c(101:120)]/nu0[,c(21:40)])^(1/3)
nu0_ratio2=nu0_ratio2
nu0_ratio3=nu0_ratio3

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

# -------------------------------------------------------------------
# Combine survival function + viability
# -------------------------------------------------------------------

## conditional germination 

gamma1.v = gamma1/(1-(1-nu_01.inter)*(1-gamma1))
gamma2.v = gamma2/(1-(1-nu_12.inter)*(1-gamma2))
gamma3.v = gamma3/(1-(1-nu_23.inter)*(1-gamma3))

f.wb = function(t,inv.b0,alpha){
  exp(-(t/inv.b0)^alpha)
}


## survival function formulation
# Oct_0.v = f.wb(t=t[x[1]],inv.b0=inv.b0.wb,alpha=a.wb)
# 
# Janpre_1.v = f.wb(t=t[x[2]],inv.b0=inv.b0.wb,alpha=a.wb)*(gamma1+(1-gamma1)*nu_01.inter)
# Jangerm_1.v = gamma1.v
# Janpost_1.v = f.wb(t=t[x[3]],inv.b0=inv.b0.wb,alpha=a.wb)*(1-gamma1)*nu_01.inter
# Oct_1.v = (f.wb(t=t[x[4]],inv.b0=inv.b0.wb,alpha=a.wb)*nu_1*(1-gamma1))
# 
# Janpre_2.v =( f.wb(t=t[x[5]],inv.b0=inv.b0.wb,alpha=a.wb)*(1-gamma1)*(gamma2+(1-gamma2)*nu_12.inter))
# Jangerm_2.v = gamma2.v
# Janpost_2.v = f.wb(t=t[x[6]],inv.b0=inv.b0.wb,alpha=a.wb)*(1-gamma1)*(1-gamma2)*nu_12.inter
# Oct_2.v = (f.wb(t=t[x[7]],inv.b0=inv.b0.wb,alpha=a.wb)*(1-gamma1)*(1-gamma2)*nu_2)
# 
# Janpre_3.v = ( f.wb(t=t[x[8]],inv.b0=inv.b0.wb,alpha=a.wb)*(1-gamma1)*(1-gamma2)*(gamma3+(1-gamma3)*nu_23.inter))
# Jangerm_3.v = gamma3.v
# Janpost_3.v = f.wb(t=t[x[9]],inv.b0=inv.b0.wb,alpha=a.wb)*(1-gamma1)*(1-gamma2)*(1-gamma3)*nu_23.inter
# Oct_3.v =(f.wb(t=t[x[10]],inv.b0=inv.b0.wb,alpha=a.wb)*(1-gamma1)*(1-gamma2)*(1-gamma3)*nu_3)

# phi_0 = Oct_0.v
# phi_1 = Janpre_1.v
# phi_2 = Janpost_1.v
# phi_3 = Oct_1.v
# phi_4 = Janpre_2.v
# phi_5 = Janpost_2.v
# phi_6 = Oct_2.v
# phi_7 = Janpre_3.v
# phi_8 = Janpost_3.v
# phi_9 = Oct_3.v

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


# if calculating survival by incorporating viability, end up with survival > 1
# by including viability, the trajectory is no longer a proper survival function
dev.off()
plot(rep(1,20),apply(phi_0-phi_1,2,min),col=ifelse(apply(phi_0-phi_1,2,min) <0,'orange','black'),ylim=c(-.25,1.25),xlim=c(.5,9.5))
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

# check distribution of difference of curves
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
# for(i in 1:20){hist(phi_0[,i]-phi_1[,i],main='');abline(v=0,col='red')}
# for(i in 1:20){hist(phi_1[,i]-phi_2[,i],main='');abline(v=0,col='red')}
# for(i in 1:20){hist(phi_2[,i]-phi_3[,i],main='');abline(v=c(0,1),col='red')}
for(i in 1:20){hist(phi_3[,i]-phi_4[,i],main='');abline(v=c(0,1),col='red')}
# for(i in 1:20){hist(phi_4[,i]-phi_5[,i],main='');abline(v=0,col='red')}
for(i in 1:20){hist(phi_5[,i]-phi_6[,i],main='');abline(v=0,col='red')}
for(i in 1:20){hist(phi_6[,i]-phi_7[,i],main='');abline(v=0,col='red')}
# for(i in 1:20){hist(phi_7[,i]-phi_8[,i],main='');abline(v=0,col='red')}
for(i in 1:20){hist(phi_8[,i]-phi_9[,i],main='');abline(v=0,col='red')}

# check distribution of components
# for(i in 1:20){hist(Oct_0.v[,i],main='');abline(v=c(0,1),col='red')}
# for(i in 1:20){hist(Janpre_1.v[,i],main='');abline(v=c(0,1),col='red')}
# for(i in 1:20){hist(Janpost_1.v[,i],main='');abline(v=c(0,1),col='red')}
#for(i in 1:20){hist(Oct_1.v[,i],main='');abline(v=c(0,1),col='red')}
for(i in 1:20){hist(Janpre_2.v[,i],main='');abline(v=c(0,1),col='red')}
#for(i in 1:20){hist(Janpost_2.v[,i],main='');abline(v=c(0,1),col='red')}
for(i in 1:20){hist(Oct_2.v[,i],main='');abline(v=c(0,1),col='red')}
for(i in 1:20){hist(Janpre_3.v[,i],main='');abline(v=c(0,1),col='red')}
#for(i in 1:20){hist(Janpost_3.v[,i],main='');abline(v=c(0,1),col='red')}
for(i in 1:20){hist(Oct_3.v[,i],main='');abline(v=c(0,1),col='red')}

# examine Janpre_2.v
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

#hist(x,breaks=100)
#hist(resample(x),col='red',breaks=100,add=TRUE)
Janpre_2.v=apply(Janpre_2.v,2,resample)
Oct_2.v=apply(Oct_2.v,2,resample)
Janpre_3.v=apply(Janpre_3.v,2,resample)
Oct_3.v=apply(Oct_3.v,2,resample)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){hist(Janpre_2.v[,i],main='',breaks=100)}#;abline(v=c(0,1),col='red')}
for(i in 1:20){hist(Oct_2.v[,i],main='',breaks=100)}#;abline(v=c(0,1),col='red')}
for(i in 1:20){hist(Janpre_3.v[,i],main='',breaks=100)}#;abline(v=c(0,1),col='red')}
for(i in 1:20){hist(Oct_3.v[,i],main='',breaks=100)}#;abline(v=c(0,1),col='red')}

# ALTERNATIVE?
# instead of resampling from the same distribution below one, instead sample
# from the transition distribution without viability?
# resample2=function(x,y){
#   tmp = ifelse(x<=1,x,NA)
#   len=sum(is.na(tmp))
#   #tmp2=x[x<=1]
#   out = ifelse(is.na(tmp),sample(y,len,replace=TRUE),tmp)
#   return(out)
# }
# 
# Janpre_2.v2=apply(Janpre_2.v,2,resample2,y=(f.wb(t=t[x[5]],inv.b0=inv.b0.wb,alpha=a.wb)*(gamma2+(1-gamma2)))/(f.wb(t=t[x[4]],inv.b0=inv.b0.wb,alpha=a.wb)))
# Oct_2.v2=apply(Oct_2.v,2,resample2,y=(f.wb(t=t[x[7]],inv.b0=inv.b0.wb,alpha=a.wb))/(f.wb(t=t[x[6]],inv.b0=inv.b0.wb,alpha=a.wb)))
# Janpre_3.v2=apply(Janpre_3.v,2,resample2,y=(f.wb(t=t[x[8]],inv.b0=inv.b0.wb,alpha=a.wb)*(gamma3+(1-gamma3)))/(f.wb(t=t[x[7]],inv.b0=inv.b0.wb,alpha=a.wb)))
# Oct_3.v2=apply(Oct_3.v,2,resample2,y=(f.wb(t=t[x[10]],inv.b0=inv.b0.wb,alpha=a.wb))/(f.wb(t=t[x[9]],inv.b0=inv.b0.wb,alpha=a.wb)))



# Re-Calculate survivorship schedule
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


# if calculating survival by incorporating viability, end up with survival > 1
# by including viability, the trajectory is no longer a proper survival function
dev.off()
plot(rep(1,20),apply(phi_0-phi_1,2,min),col=ifelse(apply(phi_0-phi_1,2,min) <0,'orange','black'),ylim=c(-.25,1.25),xlim=c(.5,9.5))
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
# Plot survival function: persistence+viability
# -------------------------------------------------------------------

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

# -------------------------------------------------------------------
# Plot discrete KM-style survival function plot
# -------------------------------------------------------------------

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
 
  text(x=3,y=.05,siteNames[i])
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

# check if there are any PMF values <0
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

# -------------------------------------------------------------------
# Calculate structured model parameters
# -------------------------------------------------------------------
dev.off()
# compare uncorrected and corrected estimates
plot(th_1[,1]*gamma1[,1],phi_1[,1]*gamma1.v[,1])
# unconditional calculation from both persistence and viability
hist(th_1*gamma1-phi_1*gamma1.v)
hist(th_4*gamma2-phi_4*gamma2.v)
hist(th_7*gamma3-phi_7*gamma3.v)

# age 2 and 3 unconditional germination estimates are lower than expected for viability 
# when looking at full distribution
plot(apply(th_1*gamma1,2,median),apply(phi_1*gamma1.v,2,median));abline(a=0,b=1)
plot(apply(th_4*gamma2,2,median),apply(phi_4*gamma2.v,2,median));abline(a=0,b=1)
plot(apply(th_7*gamma3,2,median),apply(phi_7*gamma3.v,2,median));abline(a=0,b=1)

g1=gamma1
g2=gamma2
g3=gamma3

g1.v=gamma1.v
g2.v=gamma2.v
g3.v=gamma3.v

# age 2 and 3 unconditional germination estimates are lower than expected for viability 
# when looking at full distribution
plot(apply(g1,2,median),apply(gamma1.v,2,median));abline(a=0,b=1)
plot(apply(g2,2,median),apply(gamma2.v,2,median));abline(a=0,b=1)
plot(apply(g3,2,median),apply(gamma3.v,2,median));abline(a=0,b=1)

plot(apply(g1.v,2,median),apply(gamma1.v,2,median));abline(a=0,b=1)
plot(apply(g2.v,2,median),apply(gamma2.v,2,median));abline(a=0,b=1)
plot(apply(g3.v,2,median),apply(gamma3.v,2,median));abline(a=0,b=1)

g = rbind(apply(g1,2,median),apply(g2,2,median),apply(g3,2,median))
par(mfrow=c(1,2))
plot(NA,NA,xlim=c(.5,3.5),ylim=c(0,1))
for(i in 1:20){
  lines(c(1,2,3),g[,i],type='b')
}

g.v = rbind(apply(g1.v,2,median),apply(g2.v,2,median),apply(g3.v,2,median))
plot(NA,NA,xlim=c(.5,3.5),ylim=c(0,1))
for(i in 1:20){
  lines(c(1,2,3),g.v[,i],type='b')
}

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  plot(NA,NA,xlim=c(.5,3.5),ylim=c(0,1),ylab='',xlab='',xaxt='n',yaxt='n')
  lines(c(1,2,3),g[,i],type='b',col='black')
  lines(c(1,2,3),g.v[,i],type='b',col='orange')
  
  text(x=1,y=.95,siteNames[i])
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Age (years)", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of germination", side = 2, outer = TRUE, line = 2.2)
mtext("Population-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)

dev.off()

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

par(mfrow=c(2,3))
plot(apply(s1.og,2,median),apply(s1,2,median),xlim=c(0,1),ylim=c(0,1),xlab="persistence",ylab="persistence+viability");abline(a=0,b=1);text(x=.1,y=1,"s1")
plot(apply(s2.og,2,median),apply(s2,2,median),xlim=c(0,1),ylim=c(0,1),xlab="persistence",ylab="persistence+viability");abline(a=0,b=1);text(x=.1,y=1,"s2")
plot(apply(s3.og,2,median),apply(s3,2,median),xlim=c(0,1),ylim=c(0,1),xlab="persistence",ylab="persistence+viability");abline(a=0,b=1);text(x=.1,y=1,"s3")
plot(apply(s4.og,2,median),apply(s4,2,median),xlim=c(0,1),ylim=c(0,1),xlab="persistence",ylab="persistence+viability");abline(a=0,b=1);text(x=.1,y=1,"s4")
plot(apply(s5.og,2,median),apply(s5,2,median),xlim=c(0,1),ylim=c(0,1),xlab="persistence",ylab="persistence+viability");abline(a=0,b=1);text(x=.1,y=1,"s5")
plot(apply(s6.og,2,median),apply(s6,2,median),xlim=c(0,1),ylim=c(0,1),xlab="persistence",ylab="persistence+viability");abline(a=0,b=1);text(x=.1,y=1,"s6")


# site 1 is a good example where viability increases in the last
# time step leading to increased survival relative to observed from
# persistence alone
par(mfrow=c(2,3))

hist(s1.og[,1],col='gray',breaks=100)
hist(s1[,1],col='orange',breaks=100,add=TRUE)

hist(s2.og[,1],col='gray',breaks=100)
hist(s2[,1],col='orange',breaks=100,add=TRUE)

hist(s3.og[,1],col='gray',breaks=100)
hist(s3[,1],col='orange',breaks=100,add=TRUE)

hist(s4.og[,1],col='gray',breaks=100)
hist(s4[,1],col='orange',breaks=100,add=TRUE)

hist(s5.og[,1],col='gray',breaks=100)
hist(s5[,1],col='orange',breaks=100,add=TRUE)

hist(s6.og[,1],col='gray',breaks=100)
hist(s6[,1],col='orange',breaks=100,add=TRUE)


par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  hist(s3.og[,i],col='gray',breaks=50,ylab='',xlab='',xaxt='n',yaxt='n',main='',freq=FALSE)
  hist(s3[,i],col='orange',breaks=50,add=TRUE,freq=FALSE)
  text(x=.65,y=6,siteNames[i])
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Age (years)", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of germination", side = 2, outer = TRUE, line = 2.2)
mtext("Population-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)



par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  hist(s5.og[,i],col='gray',breaks=50,ylab='',xlab='',xaxt='n',yaxt='n',main='',freq=FALSE)
  hist(s5[,i],col='orange',breaks=50,add=TRUE,freq=FALSE)
  text(x=.65,y=6,siteNames[i])
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Age (years)", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of germination", side = 2, outer = TRUE, line = 2.2)
mtext("Population-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)



# PLOTS TO MAKE
# population-level estimates of all parameters, by pop and easting (Figure 8)
# population level probability of germination conditional on viability and persistence (Figure 7)

gamma1.sum=apply(gamma1,2,quantile,probs=c(0.025,.25,.5,.75,.975))
gamma1.v.sum=apply(gamma1.v,2,quantile,probs=c(0.025,.25,.5,.75,.975))

plot(NA,NA,type='n',xlim=c(0,1),ylim=c(0,20),
     axes=FALSE,frame=FALSE,
     xlab="",ylab="")
y.pt = 20:1
for(i in 20:1){
  tmp<-gamma1.sum[,i]
  segments(x0=tmp[1],x1=tmp[5],y0=y.pt[i]-.2)
  segments(x0=tmp[2],x1=tmp[4],y0=y.pt[i]-.2,lwd=3)
  points(x=tmp[3],y=y.pt[i]-.2,pch=21,bg='white')
  
  tmp2<-gamma1.v.sum[,i]
  segments(x0=tmp2[1],x1=tmp2[5],y0=y.pt[i]+.2,col='red')
  segments(x0=tmp2[2],x1=tmp2[4],y0=y.pt[i]+.2,lwd=3,col='red')
  points(x=tmp2[3],y=y.pt[i]+.2,pch=21,bg='white',col='red')
}
axis(1,  seq(0,1,by=.2), col.ticks = 1)
axis(2, (1:20),
     labels = rev(siteNames), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
mtext("Germination probability",
      side=1,line=2.5,adj=.5,col='black',cex=1)
mtext("Age 1",
      side=3,line=0,adj=0,col='black',cex=1)

plot(NA,NA,type='n',xlim=c(0,1),ylim=c(0,1),
     axes=FALSE,frame=FALSE,
     xlab="",ylab="")
abline(a=0,b=1,col='gray',lty='dotted')
for(i in 1:20){
  tmp<-gamma1.sum[,i]
  tmp2<-gamma1.v.sum[,i]
  
  segments(x0=tmp[1],x1=tmp[5],y0=tmp2[3])
  segments(x0=tmp[2],x1=tmp[4],y0=tmp2[3],lwd=3)
  points(x=tmp[3],y=tmp2[3],pch=21,bg='white')
  
  segments(y0=tmp2[1],y1=tmp2[5],x0=tmp[3],col='red')
  segments(y0=tmp2[2],y1=tmp2[4],x0=tmp[3],lwd=3,col='red')
  points(y=tmp2[3],x=tmp[3],pch=21,bg='white',col='red')
}
axis(1,  seq(0,1,by=.2), col.ticks = 1)
axis(2,  seq(0,1,by=.2), col.ticks = 1)
mtext("Germination probability",
      side=1,line=2.5,adj=.5,col='black',cex=1)
mtext("Age 1",
      side=3,line=0,adj=0,col='black',cex=1)

gamma1.del.sum=apply(gamma1.v-gamma1,2,quantile,probs=c(0.025,.25,.5,.75,.975))

plot(NA,NA,type='n',xlim=c(0,.1),ylim=c(0,20),
     axes=FALSE,frame=FALSE,
     xlab="",ylab="")
y.pt = 20:1
for(i in 20:1){
  tmp<-gamma1.del.sum[,i]
  segments(x0=tmp[1],x1=tmp[5],y0=y.pt[i])
  segments(x0=tmp[2],x1=tmp[4],y0=y.pt[i],lwd=3)
  points(x=tmp[3],y=y.pt[i],pch=21,bg='white')
  
}
axis(1,  seq(0,.1,by=.02), col.ticks = 1)
axis(2, (1:20),
     labels = rev(siteNames), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
mtext("Delta",
      side=1,line=2.5,adj=.5,col='black',cex=1)
mtext("Age 1",
      side=3,line=0,adj=0,col='black',cex=1)


gamma1.del.sum<-data.frame(t(gamma1.del.sum),position)
names(gamma1.del.sum)[1:5] = c("ci.lolo","ci.lo","ci.med","ci.hi","ci.hihi")

plot(NA,NA,type='n',ylim=c(0,.1),xlim=c(340,375),
     axes=FALSE,frame=FALSE,
     xlab="",ylab="")

segments(x0=gamma1.del.sum$easting,y0=gamma1.del.sum$ci.lolo,y1=gamma1.del.sum$ci.hihi,lwd=1)
segments(x0=gamma1.del.sum$easting,y0=gamma1.del.sum$ci.lo,y1=gamma1.del.sum$ci.hi,lwd=3)
points(x=gamma1.del.sum$easting,y=gamma1.del.sum$ci.med,pch=21,bg='white')

axis(2, seq(0,.1,by=.02), col.ticks = 1)
axis(1, seq(340,375,by=5),
     labels = seq(340,375,by=5), las = 1, 
     col.ticks = 1, cex.axis = 1)
mtext("Delta",
      side=2,line=2.5,adj=.5,col='black',cex=1)
mtext("Easting (km)",
      side=1,line=2.5,adj=.5,col='black',cex=1)
mtext("Age 1",
      side=3,line=0,adj=0,col='black',cex=1)

# probability that seed persists and is viable after 1, 2, 3 years
# comparison of persistence vs. viability 
# estimated proportion of seeds that are viable at each time
