# -------------------------------------------------------------------
# Analysis of correlation between germination and RS
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE,max.print=100000)
# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(rjags) # jags interface
library(MCMCvis)
library(tidyverse)
library(reshape2)
library(HDInterval)
library(bayesplot)
library(rethinking)

# -------------------------------------------------------------------
# Functions for use when analyzing data
# -------------------------------------------------------------------
temporal_variance <- function(x,fun=var){
  apply(x,1,fun)
}

cols_fun <- function(x,fun=var){
  apply(x,2,fun)
}

# geometric var
gsd <- function(x){
  y <- exp(sd(log(x)))
  return(y)
}

# sample based calculation of gsd
gsd.am <- function(x){
#  x=x+.5
  n = length(x[!is.na(x)])
  mu = exp(mean(log(x),na.rm=TRUE))
  y <- exp(sqrt(sum((log(x/mu))^2,na.rm=TRUE)/(n-1)))
  return(y)
}

f<-function(x="param"){
  chain<-MCMCchains(zc,params = x)
  BCI <- t(apply(p,2,FUN = function(x) quantile(x, c(.025, .5, .975))))
  return(BCI)
}

posterior.mode = function(x){
  x.max=max(x)
  x.min=min(x)
  dres <- density( x ,from = x.min, to = x.max)
  modeParam <- dres$x[which.max(dres$y)]
  return(modeParam)
}

# -------------------------------------------------------------------
# Read in samples from posterior distributions
# -------------------------------------------------------------------

g1 <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/g1-pop.RDS")
sigma <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/sigma-analysis-noPool.RDS")
fec <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/tfe-analysis-noPool.RDS")
phi <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/phi-analysis-noPool.RDS")

rsPosterior = sigma*fec*phi

# get mode from full posterior and calculate var in RS based on modes
sigma.mode = apply(sigma,2,posterior.mode);
fec.mode = apply(fec,2,posterior.mode);
phi.mode = apply(phi,2,posterior.mode);

rs.mode=apply(sigma*fec*phi,2,posterior.mode)

n.iter = dim(rsPosterior)[1]
# -------------------------------------------------------------------
# Compare calculation of posterior by mode of components vs. mode of RS
# -------------------------------------------------------------------

par(mfrow = c(4,5),
        oma = c(5,4,0,0) + 0.1,
        mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
    obj = sigma.mode
    index=grep(paste0("\\[",i,","),names(obj))
  plot(sigma.mode[index]*fec.mode[index]*phi.mode[index],rs.mode[index],pch=19);
  abline(a=0,b=1,col='red');text(sigma.mode[index]*fec.mode[index]*phi.mode[index],rs.mode[index],index)
}

# -------------------------------------------------------------------
# Compare calculation of geometric SD in reproductive success
# -------------------------------------------------------------------
par(mfrow = c(1,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
plot(NA,NA,xlim=c(0,10),ylim=c(0,10),pch=19);

for(i in 1:20){
  obj = sigma.mode
  index=grep(paste0("\\[",i,","),names(obj))
  text(gsd.am(sigma.mode[index]*fec.mode[index]*phi.mode[index]),gsd.am(rs.mode[index]),i);
    abline(a=0,b=1,col='red')
}

# -------------------------------------------------------------------
# Calculate geometric SD in reproductive success
# -------------------------------------------------------------------

df.list = list()
for(i in 1:20){
  obj = rsPosterior
  index=grep(paste0("\\[",i,","),colnames(obj))
  tmp.df=apply(obj[,index],1,gsd.am)
  df.list[[i]] = tmp.df
}

gsdSummary=do.call(cbind,df.list)

# -------------------------------------------------------------------
# Analyze correlation of germination and GSD per-capita RS
# -------------------------------------------------------------------
# create empty vector for the correlation
posterior.correlation<-c()

# calculate correlation for each draw from the posterior
for(i in 1:n.iter){
  posterior.correlation[i]<-cor(g1[i,],gsdSummary[i,])
}

# calculate the 95% credible interval and HPDI for g1
CI.g1 <- apply(g1,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
HPDI.g1 <- apply(g1,2,FUN = function(x) hdi(x, .95))

# calculate the 95% credible interval and HPDI for probability of RS
CI.rs <- apply(gsdSummary,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
HPDI.rs <- apply(gsdSummary,2,FUN = function(x) hdi(x, .95))

# put medians and credible intervals into data frame
g1PosteriorSummary <- data.frame(cbind(t(HPDI.g1),CI.g1[2,]))
names(g1PosteriorSummary) <- c("lo.g1","hi.g1","med.g1")

rsPosteriorSummary<-data.frame(t(HPDI.rs),CI.rs[2,])
names(rsPosteriorSummary) <- c("lo.rs","hi.rs","med.rs")

# calculate the 95% credible interval and HPDI for the correlation
CI.correlation <- quantile(posterior.correlation, c(.025, .5, .975))
HPDI.correlation <- hdi(posterior.correlation, c(.95))
HPDI.correlation.lo <- hdi(posterior.correlation, c(.5))


pdf(
  "~/Dropbox/clarkiaSeedBanks/products/figures/analysis/correlation-germ-rs.pdf",
  height = 8, width = 6)
par(mar=c(4,4,2,1))
par(fig=c(0,10,4,10)/10)
# plot median of g1 vs. median of RS with CIs
plot(x = NA,
     y = NA,
     xlim=c(0,12),ylim=c(0,.6),
     pch=16, cex = 0.5,
     xlab = "",
     ylab = "",
     xaxt= "n", yaxt="n",
     cex.lab = 1, cex.axis = 1)

segments(x0=rsPosteriorSummary$lo.rs,x1=rsPosteriorSummary$hi.rs,
         y0=g1PosteriorSummary$med.g1, lwd=1)
segments(x0=rsPosteriorSummary$med.rs,
         y0=g1PosteriorSummary$lo.g1, y1=g1PosteriorSummary$hi.g1,
         lwd=1)
points(rsPosteriorSummary$med.rs,g1PosteriorSummary$med.g1,
       pch=21,col='black',bg='white',cex=1.25)

axis(1, seq(0,12,by=2),
     labels = seq(0,12,by=2), las = 1, line = 0,
     col = NA, col.ticks = 1, cex.axis = 1)
axis(2, seq(0,1,by=.2),
     labels = seq(0,1,by=.2), las = 1, line = 0,
     col = NA, col.ticks = 1, cex.axis = 1)
mtext("Germination probability",
      side=2,line=2.5,adj=.5,col='black',cex=1)
mtext("Geometric SD RS",
      side=1,line=2,adj=.5,col='black',cex=1)

text(x=2.5,y=.58,
     paste0("Pearson's r=",round(CI.correlation[2],2)),
     cex=1)
#abline(a=0,b=1)

par(fig=c(0,10,0,4.5)/10)
par(new=T)
# plot posterior of correlation coefficient
hist(posterior.correlation,breaks = 50, 
     main = "", xlab = "", ylab='', xaxt='n',yaxt='n', 
     xlim = c(-1, 1), ylim=c(0,2.5), border= "white",
     freq = FALSE, col = "gray75", 
     cex.lab = 1.25,cex.axis=1.5)

# as in Duskey dissertation
segments(x0=HPDI.correlation[1],x1=HPDI.correlation[2],y0=2.4,lwd=1.25)
segments(x0=HPDI.correlation.lo[1],x1=HPDI.correlation.lo[2],y0=2.4,lwd=2.5)
points(x=median(posterior.correlation),y=2.4,pch=21,col='black',bg='white',cex=1.25)

axis(1, seq(-2,2,by=.2),
     labels = seq(-2,2,by=.2), las = 1, line = 0,
     col = NA, col.ticks = 1, cex.axis = 1)
axis(2, seq(0,3,by=.4),
     labels = seq(0,3,by=.4), las = 1, line = -.5,
     col = NA, col.ticks = 1, cex.axis = 1)
segments(x0=-1,x1=1,y0=-.1,lwd=2)
segments(x0=-1.035,y0=0,y1=2.4,lwd=1.5)
mtext("Density",
      side=2,line=2.5,adj=.5,col='black',cex=1)
mtext("Correlation of germination and  GSD RS",
      side=1,line=2,adj=.5,col='black',cex=1)
dev.off()

# -------------------------------------------------------------------
# Uncomment below to label points
# -------------------------------------------------------------------

# df <- g1PosteriorSummary %>%
#   dplyr::bind_cols(site=siteNames) %>%
#   dplyr::left_join(rsPosteriorSummary %>%
#                      dplyr::bind_cols(site=siteNames),by="site")
# 
# library(ggrepel)
# 
# g1 <- ggplot(df,aes(x=med.rs,y=med.g1,label=site)) +
#   geom_point() +
#   geom_text_repel(size=3,color="black") +
#   # annotate("text", label =  paste0("Pearson's r=",round(CI.correlation[1],2)), x = 2.5, y = .29, size = 4) +
#   theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
#   # scale_x_continuous(limits = c(.99,8), expand = c(0, 0), breaks = c(1, 3, 5, 7)) +
#   # scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
#   xlab("Geometric SD of reproductive success") +
#   ylab("Mean germination probability [P(G)]") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# ggsave(filename=  "~/Dropbox/clarkiaSeedBanks/products/figures/analysis/correlation-germ-rs-labeled.pdf",
#        plot=g1,width=4,height=4)
# 

# -------------------------------------------------------------------
# Repeat with low fitness years set to 0
# -------------------------------------------------------------------

names<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>%
  dplyr::select(site)
siteNames = unique(names$site)

siteIndex <- data.frame(site=siteNames,siteIndex=1:20)
yearIndex <- data.frame(year=2006:2020,yearIndex=1:15)

lowFitnessYears <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/output/lowFitnessYears.RDS")

lowFitnessYears=lowFitnessYears %>% 
  dplyr::left_join(yearIndex) %>% 
  dplyr::left_join(siteIndex)

# set years without any plants at the population-level to zero fitness
df.list = list()
for(i in 1:20){
  obj = rsPosterior
  zeroYears=lowFitnessYears[lowFitnessYears$siteIndex==i,]$yearIndex
  index=grep(paste0("\\[",i,","),colnames(obj))
  obj=obj[,index]
  obj[,zeroYears] = 0
  tmp.df=apply(obj,1,gsd.am)
  df.list[[i]] = tmp.df
}

gsdSummary=do.call(cbind,df.list)

n.iter=dim(gsdSummary)[1]

# create empty vector for the correlation
posterior.correlation<-c()

# calculate correlation for each draw from the posterior
for(i in 1:n.iter){
  posterior.correlation[i]<-cor(g1[i,],gsdSummary[i,])
}

# calculate the 95% credible interval and HPDI for g1
CI.g1 <- apply(g1,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
HPDI.g1 <- apply(g1,2,FUN = function(x) hdi(x, .95))

# calculate the 95% credible interval and HPDI for probability of RS
CI.rs <- apply(gsdSummary,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
HPDI.rs <- apply(gsdSummary,2,FUN = function(x) hdi(x, .95))

# put medians and credible intervals into data frame
g1PosteriorSummary <- data.frame(cbind(t(HPDI.g1),CI.g1[2,]))
names(g1PosteriorSummary) <- c("lo.g1","hi.g1","med.g1")

rsPosteriorSummary<-data.frame(t(HPDI.rs),CI.rs[2,])
names(rsPosteriorSummary) <- c("lo.rs","hi.rs","med.rs")

# calculate the 95% credible interval and HPDI for the correlation
CI.correlation <- quantile(posterior.correlation, c(.025, .5, .975))
HPDI.correlation <- hdi(posterior.correlation, c(.95))
HPDI.correlation.lo <- hdi(posterior.correlation, c(.5))


pdf(
  "~/Dropbox/clarkiaSeedBanks/products/figures/analysis/correlation-germ-rs-lowfitness.pdf",
  height = 8, width = 6)
par(mar=c(4,4,2,1))
par(fig=c(0,10,4,10)/10)
# plot median of g1 vs. median of RS with CIs
plot(x = NA,
     y = NA,
     xlim=c(0,12),ylim=c(0,.6),
     pch=16, cex = 0.5,
     xlab = "",
     ylab = "",
     xaxt= "n", yaxt="n",
     cex.lab = 1, cex.axis = 1)

segments(x0=rsPosteriorSummary$lo.rs,x1=rsPosteriorSummary$hi.rs,
         y0=g1PosteriorSummary$med.g1, lwd=1)
segments(x0=rsPosteriorSummary$med.rs,
         y0=g1PosteriorSummary$lo.g1, y1=g1PosteriorSummary$hi.g1,
         lwd=1)
points(rsPosteriorSummary$med.rs,g1PosteriorSummary$med.g1,
       pch=21,col='black',bg='white',cex=1.25)

axis(1, seq(0,12,by=2),
     labels = seq(0,12,by=2), las = 1, line = 0,
     col = NA, col.ticks = 1, cex.axis = 1)
axis(2, seq(0,1,by=.2),
     labels = seq(0,1,by=.2), las = 1, line = 0,
     col = NA, col.ticks = 1, cex.axis = 1)
mtext("Germination probability",
      side=2,line=2.5,adj=.5,col='black',cex=1)
mtext("Geometric SD RS",
      side=1,line=2,adj=.5,col='black',cex=1)

text(x=2.5,y=.58,
     paste0("Pearson's r=",round(CI.correlation[2],2)),
     cex=1)
#abline(a=0,b=1)

par(fig=c(0,10,0,4.5)/10)
par(new=T)
# plot posterior of correlation coefficient
hist(posterior.correlation,breaks = 50, 
     main = "", xlab = "", ylab='', xaxt='n',yaxt='n', 
     xlim = c(-1, 1), ylim=c(0,2.5), border='white',
     freq = FALSE, col = "gray75", 
     cex.lab = 1.25,cex.axis=1.5)

# as in Duskey dissertation
segments(x0=HPDI.correlation[1],x1=HPDI.correlation[2],y0=2.4,lwd=1.25)
#segments(x0=CI.correlation[2],y0=2.3,y1=2.5,lwd=1.5)
segments(x0=HPDI.correlation.lo[1],x1=HPDI.correlation.lo[2],y0=2.4,lwd=2.5)
points(x=median(posterior.correlation),y=2.4,pch=21,col='black',bg='white',cex=1.25)
#segments(x0=CI.correlation[2],y0=2.4,y1=0,lwd=2,lty='dotted')

axis(1, seq(-2,2,by=.2),
     labels = seq(-2,2,by=.2), las = 1, line = 0,
     col = NA, col.ticks = 1, cex.axis = 1)
axis(2, seq(0,3,by=.4),
     labels = seq(0,3,by=.4), las = 1, line = -.5,
     col = NA, col.ticks = 1, cex.axis = 1)
segments(x0=-1,x1=1,y0=-.1,lwd=2)
segments(x0=-1.035,y0=0,y1=2.4,lwd=1.5)
mtext("Density",
      side=2,line=2.5,adj=.5,col='black',cex=1)
mtext("Correlation of germination and  GSD RS",
      side=1,line=2,adj=.5,col='black',cex=1)
dev.off()

# -------------------------------------------------------------------
# Uncomment below to label points
# -------------------------------------------------------------------
# 
# df <- g1PosteriorSummary %>%
#   dplyr::bind_cols(site=siteNames) %>%
#   dplyr::left_join(rsPosteriorSummary %>%
#                      dplyr::bind_cols(site=siteNames),by="site")
# 
# library(ggrepel)
# 
# g1 <- ggplot(df,aes(x=med.rs,y=med.g1,label=site)) +
#   geom_point() +
#   geom_text_repel(size=3,color="black") +
#  # annotate("text", label =  paste0("Pearson's r=",round(CI.correlation[1],2)), x = 2.5, y = .29, size = 4) +
#  theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
#  # scale_x_continuous(limits = c(.99,8), expand = c(0, 0), breaks = c(1, 3, 5, 7)) +
#  # scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
#   xlab("Geometric SD of reproductive success") +
#   ylab("Mean germination probability [P(G)]") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# ggsave(filename=  "~/Dropbox/clarkiaSeedBanks/products/figures/analysis/correlation-germ-rs-lowfitness-labeled.pdf",
#        plot=g1,width=4,height=4)

# compare RS calculations
# 
# rs <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/rsPosterior.RDS")
# lowFitnessYears <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/output/lowFitnessYearsPlots.RDS")
# lowFitnessYears=lowFitnessYears %>% dplyr::left_join(yearIndex) %>% dplyr::left_join(siteIndex)
# 
# df.list = list()
# for(i in 1:20){
#   obj = rs
#   zeroYears=lowFitnessYears[lowFitnessYears$siteIndex==i,]$yearIndex
#   index=grep(paste0("\\[",i,","),colnames(obj))
#   obj=obj[,index]
#   obj[,zeroYears] = 0
#   tmp.df=apply(obj,1,gsd.am)
#   df.list[[i]] = tmp.df
# }
# 
# gsdSummaryLowFitness=do.call(cbind,df.list)
# 
# # calculate the 95% credible interval and HPDI for probability of RS
# CI.rsLowFitness <- apply(gsdSummaryLowFitness,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
# 
# 
# df.list = list()
# for(i in 1:20){
#   obj = rs
#   index=grep(paste0("\\[",i,","),colnames(obj))
#   tmp.df=apply(obj[,index],1,gsd.am)
#   df.list[[i]] = tmp.df
# }
# 
# gsdSummary=do.call(cbind,df.list)
# 
# # calculate the 95% credible interval and HPDI for probability of RS
# CI.rs <- apply(gsdSummary,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
# 
# 
# plot(CI.rsLowFitness[2,],CI.rs[2,],
#      xlim=c(2,12),ylim=c(2,12),
#      pch=19)
# abline(a=0,b=1)
# 
# delta = CI.rsLowFitness[2,]-CI.rs[2,]
# index=order(delta)
# delta=delta[index]
# plot(NA,NA,type='n',ylim=c(0,20),xlim=c(0,8),
#      axes=FALSE,frame=FALSE,
#      xlab="",ylab="")
# 
# for(i in 1:20){
#   tmp<-delta[i]
#  # segments(y0=tmp[1],y1=tmp[5],x0=y.pt[i])
#   #segments(y0=tmp[2],y1=tmp[4],x0=y.pt[i],lwd=3)
#   points(x=tmp,y=i,pch=21,bg='white')
# }
# axis(1,  seq(0,8,by=1), col.ticks = 1)
# axis(2, 1:20,
#      labels = (siteNames[index]), las = 2, 
#      col = NA, col.ticks = 1, cex.axis = 1)
# 
# position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
#   dplyr::select(site,easting,dominant.surface.rock.type) %>%
#   dplyr::mutate(easting=easting/1000)
# 
# delta = CI.rsLowFitness[2,]-CI.rs[2,]
# df=data.frame(delta,siteIndex) %>%
#   dplyr::left_join(position,by='site')
# 
# plot(df$easting,df$delta,type='n')
# text(x=df$easting,y=df$delta,labels=df$site)
# 
# rs.var.lowFitness = CI.rsLowFitness[2,]
# df=data.frame(rs.var.lowFitness,siteIndex) %>%
#   dplyr::left_join(position,by='site')
# 
# plot(df$easting,df$rs.var.lowFitness,type='n',ylim=c(0,20))
# text(x=df$easting,y=df$rs.var.lowFitness,labels=df$site)
# 
# rs.var = CI.rs[2,]
# df=data.frame(rs.var,siteIndex) %>%
#   dplyr::left_join(position,by='site')
# 
# plot(df$easting,df$rs.var,type='n',ylim=c(0,10))
# text(x=df$easting,y=df$rs.var,labels=df$site)
# # 
# # 
# # g1.blank <- ggplot(df,aes(x=var.rs,y=med,label=site)) +
# #  # geom_point() +
# #   # geom_text_repel(size=3,color="black") +
# #  # annotate("text", label =  paste0("Pearson's r=",round(CI.correlation[1],2)), x = 1.5, y = .29, size = 4) +
# #   theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
# #   scale_x_continuous(limits = c(.99,8), expand = c(0, 0), breaks = c(1, 3, 5, 7)) +
# #   scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
# #   xlab("Geometric SD of reproductive success") +
# #   ylab("Mean germination probability [P(G)]") +
# #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# # 
# # ggsave(filename=  "~/Dropbox/clarkiaSeedBanks/products/figures/germ_rs_correlation-blank.pdf",
# #        plot=g1.blank,width=4,height=4)
# # 
# # 
# # # g1.point <- ggplot(df %>% dplyr::filter(site=="BG"),aes(x=var.rs,y=med,label=site)) +
# # #    geom_point() +
# # #   # geom_text_repel(size=3,color="black") +
# # #   theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
# # #   scale_x_continuous(limits = c(.99,8), expand = c(0, 0), breaks = c(1, 3, 5, 7)) +
# # #   scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
# # #   xlab("Geometric SD of reproductive success") +
# # #   ylab("Mean germination probability [P(G)]") +
# # #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# # # 
# # # ggsave(filename=  "~/Dropbox/clarkiaSeedBanks/products/figures/germ_rs_correlation-point.pdf",
# # #        plot=g1.point,width=4,height=4)
# # 
# # matplot(c0[,2],c0[,3:12],
# #         xlab="Mean site germination probability",
# #         ylab="Annual mean fitness",
# #         main=expression(paste("cor(",sigma,",P(G))=-.75")),
# #         pch=16,col='gray')
# # plot(c0[1:nsites,1],c0[1:nsites,2],xlab="Variance in Fitness",ylab="Probability of Germination")
# # 
# # set.seed(11)
# # 
# # # COPULA
# # 
# # library(MASS)
# # 
# # nsites = 20
# # nyears = 10
# # 
# # mu = mu = rnorm(n=nsites,mean=25,sd=0)
# # 
# # sim <- function(rho=0,nsites=20){
# #   # Defines the  sequence for stochastic trials.
# #   reps=1000
# #   
# #   # a and b are hyperparameters of the gamma distribution 
# #   # that define both the expected value and variance.   
# #   a = 6
# #   b = 1.75
# #   
# #   # alpha and beta are hyperparameters of the beta distribution that define both the expected value and
# #   # variance.  
# #   alpha =  1
# #   beta =  1
# #   
# #   # Defines the temporal correlation between the two parameters.
# #   rho = rho
# #   
# #   # Generates standard multivariate normal data with correlation structure defined by rho.
# #   Z <- mvrnorm(n=reps,mu=c(0,0), 
# #                matrix(data=c(1,rho,rho,1),
# #                       nrow=2,ncol=2))
# #   
# #   # Apply the Normal CDF function to Z to obtain data that is uniform on the interval [0,1], but still correlated.
# #   U = pnorm(Z)
# #   
# #   # x is gamma distributed
# #   X = qgamma(U[,1],shape=a,rate=b) 
# #  # X = qunif(U[,1],0,.75) 
# #   
# #   # y is beta distributed
# #   #X <- cbind(X,qbeta(U[,2],shape1=alpha,shape2=beta) )
# #   X <- cbind(X,qunif(U[,2],0,.3) )
# #   
# #   # gamma marginal of multivariate X.
# #   # hist(X[,1])
# #   # beta marginal of multivariate X
# #   # hist(X[,2])
# #   
# #   # plot(X[,1],X[,2])
# #   
# #   nsites = nsites
# #   nyears = 10
# #   X[1:nsites,1:2]
# #   
# #   # generate samples for each site with a different variance
# #   d<-lapply(X[1:nsites,1],rnorm,n=nyears,mean=rnorm(n=nsites,mean=mu,sd=0))
# #   # check that the variances are appropriately sampled
# #   cbind(unlist(lapply(d, sd)),X[1:nsites,1])
# #   
# #   d<-data.frame(d)
# #   names(d) <- 1:20
# #   
# #   # check variances again
# #   cbind(apply(d,2, sd),X[1:nsites,1])
# #   
# #   dim(d)
# #   
# #   dt<-cbind(X[1:nsites,1:2],t(d))
# #   return(dt)
# # }
# # 
# # 
# # 
# # c0<-sim(rho=-.75,nsites=20)
# # 
# # df.sim <- data.frame(df$site,c0)
# # g1.hypothesis<-ggplot(df.sim,aes(x=X,y=V2)) +
# #   geom_point(color='gray') +
# #   # geom_text_repel(size=3,color="black") +
# #    annotate("text", label =  paste0("Pearson's r=",-.75), x = 6, y = .29, size = 4,color='gray') +
# #   theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
# #   scale_x_continuous(limits = c(.99,8), expand = c(0, 0), breaks = c(1, 3, 5, 7)) +
# #   scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
# #   xlab("Geometric SD of reproductive success") +
# #   ylab("Mean germination probability [P(G)]") +
# #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# # 
# # ggsave(filename=  "~/Dropbox/clarkiaSeedBanks/products/figures/germ_rs_correlation-hypothesis.pdf",
# #        plot=g1.hypothesis,width=4,height=4)
# # 
# # 
# # g1.point<-ggplot(df.sim[1,],aes(x=X,y=V2)) +
# #   geom_segment(aes(x=X,y=0,xend=X,yend=V2),color='#E69F00') +
# #   geom_segment(aes(x=1,y=V2,xend=X,yend=V2),color='#CC79A7') +
# #   
# #   geom_point(color='gray',size=5) +
# #   # geom_text_repel(size=3,color="black") +
# #  # annotate("text", label =  paste0("Pearson's r=",-.75), x = 6, y = .29, size = 4,color='gray') +
# #   theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
# #   scale_x_continuous(limits = c(.99,8), expand = c(0, 0), breaks = c(1, 3, 5, 7)) +
# #   scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
# #   xlab("Geometric SD of reproductive success") +
# #   ylab("Mean germination probability [P(G)]") +
# #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# # 
# # ggsave(filename=  "~/Dropbox/clarkiaSeedBanks/products/figures/germ_rs_correlation-point.pdf",
# #        plot=g1.point,width=4,height=4)
# # 
# # 
# # # # Germination
# # # # extract parameters for analysis
# # #  posterior.g1<-MCMCchains(belowground,params = "g1")
# # # 
# # # # calculate the 95% credible interval and HPDI for g1
# # # CI.g1 <- apply(posterior.g1,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
# # # HPDI.g1 <- apply(posterior.g1,2,FUN = function(x) rethinking::HPDI(x, .95))
# # # 
# # # ## Reproductive success
# # # lapply(aboveground,)
# # # 
# # # d<-lapply(aboveground,temporal_variance,fun=gsd)
# # # reproductiveSuccess<-matrix(unlist(d), ncol = 20, byrow = FALSE)
# # # # # exclude problem sites
# # # # reproductiveSuccess<-reproductiveSuccess[,-probs]
# # # 
# # # # don't currently have the full distribution?
# # # CI.reproductiveSuccess <- apply(reproductiveSuccess,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
# # # HPDI.reproductiveSuccess <- apply(reproductiveSuccess,2,FUN = function(x) hdi(x, .95))
# # # 
# # # g1PosteriorSummary<-data.frame(t(CI.g1))
# # # names(g1PosteriorSummary) <- c("lo.g1","med.g1","hi.g1")
# # # 
# # # ## need to finish this section
# # # reproductiveSuccessPosteriorSummary<-data.frame(t(CI.reproductiveSuccess))
# # # names(reproductiveSuccessPosteriorSummary) <- c("lo.rs","med.rs","hi.rs")
# # # 
# # # # calculate correlation for each draw from the posterior
# # # n.iter=3000
# # # posterior.correlation<-c()
# # # 
# # # for(i in 1:n.iter){
# # #   posterior.correlation[i]<-cor(posterior.g1[i,],reproductiveSuccess[i,])
# # # }
# # # 
# # # # calculate the 95% credible interval and HPDI for the correlation
# # # CI.correlation <- quantile(posterior.correlation, c(.025, .5, .975))
# # # HPDI.correlation <- hdi(posterior.correlation, .95)
# # # 
# # # 
# # #   
# # # 
# # # ggsave(filename="~/Dropbox/clarkiaSeedBanks/products/figures/germ_rs_correlation-labeled.pdf",
# # #        plot=g1,width=6,height=6)
# # # 
# # # 
# # # 
# # # hist(posterior.correlation,breaks = 50, main = "", xlab = "", xlim = c(-1, 1),
# # #      freq = FALSE, col = "azure1", cex.lab = 1.5,cex.axis=1.5)
# # # 
# # # title(xlab="Correlation of germination and geometric SD of fitness \n (Pearson's r)", line=4, cex.lab=1.5)
# # # 
# # # abline(v=CI.correlation[c(1,3)],lty='dashed',lwd='2')
# # # abline(v=CI.correlation[2],lty='solid',lwd='2',col='red')
# 
# 
# # -------------------------------------------------------------------
# # Code below was comparing calculation of variance in reproductive success
# # -------------------------------------------------------------------
# 
# f.rs = function(s,f,p){
#   tmp = s*f*p
#   df.list = list()
#   for(i in 1:20){
#     obj = tmp
#     index=grep(paste0("\\[",i,","),names(obj))
#     tmp.df=gsd.am(obj[index])
#     df.list[[i]] = tmp.df
#   }
#   unlist(df.list)
# }
# 
# # calculate rs. mode based on modes of components
# rs.mode = sigma.mode*fec.mode*phi.mode
# 
# df.list = list()
# for(i in 1:20){
#   obj = rs.mode
#   index=grep(paste0("\\[",i,","),names(obj))
#   tmp.df=gsd.am(obj[index])
#   df.list[[i]] = tmp.df
# }
# 
# gsdSummary=unlist(df.list)
# 
# 
# for(i in 1:20){
#   obj = sigma.mode*fec.mode*phi.mode
#   index=grep(paste0("\\[",i,","),colnames(obj))
#   tmp.df=apply(obj[,index],1,gsd.am)
#   df.list[[i]] = tmp.df
# }
# 
# n.iter = dim(sigma)[1]
# gsdSummaryComponent=matrix(NA,ncol=20,nrow=n.iter)
# 
# for(i in 1:20){gsdSummaryComponent[,i]=df.list[[i]]}
# 
# 
# par(mfrow = c(4,5),
#     oma = c(5,4,0,0) + 0.1,
#     mar = c(0,0,1,1) + 0.1)
# for(i in 1:20){
#   hist(df.list[[i]],breaks=100);abline(v=gsdSummary[i],col='red',lwd=2)
# }
# 
# # get mode from resampling posterior and calculate var in RS based on modes
# n.iter=dim(sigma)[1]
# n.sub = 1000
# n.rep=100
# 
# # sig.samples=sample(n.iter,n.rep)
# # fec.samples=sample(n.iter,n.rep)
# # phi.samples=sample(n.iter,n.rep)
# 
# # repeat calculation from above using subsamples from the posterior
# rs.mat=matrix(NA,nrow=n.rep,ncol=20)
# for(i in 1:n.rep){
#   sig.samples=sample(1:n.iter,n.sub)
#   fec.samples=sample(1:n.iter,n.sub)
#   phi.samples=sample(1:n.iter,n.sub)
#   
#   sigma.mode.tmp = apply(sigma[sig.samples,],2,posterior.mode);
#   fec.mode.tmp = apply(fec[fec.samples,],2,posterior.mode);
#   phi.mode.tmp = apply(phi[phi.samples,],2,posterior.mode);
#   
#   rs.mat[i,]=f.rs(sigma.mode.tmp,fec.mode.tmp,phi.mode.tmp)
# }
# 
# # variance in reproductive success based on modes alone is higher than
# # variance in reproductive success based on full distribution
# # this is likely because when using the full distribution there is some chance
# # of capturing high probability years
# # issue in site 10, 12, 20, 17, 18
# plot(rs.mat)
# 
# # par(mfrow=c(4,5),oma=c(0,0,0,0),mar=c(0,0,0,0))
# # for(i in 1:20){
# # plot(rs.mat[,i],ylim=c(0,10));abline(h=gsdSummary[i],col="red",lwd=2)
# # }
# # 
# # for(i in 1:20){hist(rs.mat[,i],breaks=25,main="");abline(v=gsdSummary[i],col="red",lwd=2)}
# # n.iter=dim(gsdSummary)[1]