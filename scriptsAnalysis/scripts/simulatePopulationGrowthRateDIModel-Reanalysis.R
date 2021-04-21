# -------------------------------------------------------------------
# Density-independent model of germination
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE,max.print=100000)
# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(rjags) # jags interface
library(MCMCvis)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(HDInterval)
library(bayesplot)

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


f<-function(x="parm",chain){
  chain<-MCMCchains(chain=belowground,params = x)
  p<-boot::inv.logit(chain)
  BCI <- t(apply(p,2,FUN = function(x) quantile(x, c(.025, .5, .975))))
  return(BCI)
}

# geometric mean from 
# https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
gm_mean = function(x, na.rm=TRUE){
exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x), na.rm=na.rm) / length(x))
}

# geometric mean from +.5 to 0
# gm_mean = function(x, na.rm=TRUE){
#   exp(sum(log(x+.5)) / length(x))
# }


# -------------------------------------------------------------------

# read in samples from posterior distributions
s0 <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/s0-pop.RDS")
s1 <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/s1-pop.RDS")
g1 <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/g1-pop.RDS")
s2 <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/s2-pop.RDS")
s3 <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/s3-pop.RDS")
rs <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/rsPosterior.RDS")

names<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>%
  dplyr::select(site)
siteNames = unique(names$site)

siteIndex <- data.frame(site=siteNames,siteIndex=1:20)
yearIndex <- data.frame(year=2006:2020,yearIndex=1:15)
index=expand.grid(1:20,1:15)

# -------------------------------------------------------------------

fitness <- function(g=g1,s0=s0,s1=s1,s2=s2,s3=s3,rs=rs){
  p1 = g*rs*s1*s0
  p2 = (1-g)*(s2*s3)
  return(as.numeric(p1+p2))
}

posterior.median <- function(x){return(apply(x,2,median))}

posterior.mode = function(x){
  x.max=max(x)
  x.min=min(x)
  dres <- density( x ,from = x.min, to = x.max)
  modeParam <- dres$x[which.max(dres$y)]
  return(modeParam)
}

modeEst <- function(x){return(apply(x,2,posterior.mode))}

# use mode as Bayesian estimator
s0.hat = modeEst(s0)
g1.hat = modeEst(g1)
s1.hat = modeEst(s1)
s2.hat = modeEst(s2)
s3.hat = modeEst(s3)
rs.hat = modeEst(rs)

# -------------------------------------------------------------------
# calculate autocorrelation of per-capita reproductive success (mode)
# no significant autocorrelation, include in appendix
# -------------------------------------------------------------------

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for( k in 1:20){
  pop.index=index[,1]==k
  rs.tmp = rs.hat[pop.index]
  yt = rs.tmp
  acf(yt)
}

# -------------------------------------------------------------------
# calculate harmonic mean for each population
# -------------------------------------------------------------------

f = function(x){1/mean(1/x)}

pop.harmonicMean =c ()
for( k in 1:20){
  pop.index=index[,1]==k
  rs.tmp = c(rs.hat[pop.index])
  hm = f(rs.tmp)
  pop.harmonicMean[k] = hm
}

par(mfrow=c(1,1))
plot(pop.harmonicMean,s2.hat*s3.hat,type='n');
abline(a=0,b=1);
text(pop.harmonicMean,s2.hat*s3.hat,siteNames)

# -------------------------------------------------------------------
# Simulation to calculate optimal value of germination
# -------------------------------------------------------------------

g = seq(0,1,by=.01)
g.seq = g
fit = c()
g.mat = matrix(NA,nrow=1000,ncol=length(g))
g.sites = list()
yt = c()
n.iter=dim(rs)[1]

# -------------------------------------------------------------------
# for each population k (use mode)
# get the median of per-capita reproductive success in each year
# sample a sequence of 1000 years to use in the analysis for each value of g
# for each g, calculate the growth rate for 1000 years
# -------------------------------------------------------------------

for( k in 1:20){
  rs.tmp = rs[,grep(paste0("\\[",k,","),colnames(rs))]
  rs.tmp = rs.tmp[sample(1:4500,1),]
#  rs.tmp = rs.hat[grep(paste0("\\[",k,","),names(rs.hat))]
 # pop.index=index[,1]==k
 # rs.tmp = rs.hat[pop.index]
  yt = sample(rs.tmp,1000,replace=TRUE)
    for( i in 1:length(g)){
        fit<-fitness(g=g[i],s0.hat[k],s1.hat[k],s2.hat[k],s3.hat[k],rs=yt)
        #logfit<-log(fit)
        g.mat[,i] <- fit
        #g.mat[,i] <- logfit
    }
  g.sites[[k]] <- g.mat
}

par(mfrow=c(1,1))
#plot(g.seq,apply(g.mat,2,gm_mean))

apply(g.mat,2,prod)
# -------------------------------------------------------------------
# Plots of geometric mean fitness
# note that without reproductive failure fitness
# remains high at high values of g
# -------------------------------------------------------------------

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(k in 1:20){
  plot(g.seq,apply(g.sites[[k]],2,gm_mean),type='l')
}

# -------------------------------------------------------------------
# Plots of arithmetic mean fitness (log)
# -------------------------------------------------------------------

for(k in 1:20){
  plot(g.seq,log(apply(g.sites[[k]],2,mean)),type='l')
}

# -------------------------------------------------------------------
# plots of standard deviation of pop growth rate
# -------------------------------------------------------------------

for(k in 1:20){
  plot(g.seq,apply(g.sites[[k]],2,sd),type='l')
}

# -------------------------------------------------------------------
# Calculate optimal germination fraction
# -------------------------------------------------------------------
gmean <- function(x){apply(x,2,gm_mean)}
site.optima<-lapply(g.sites,gmean)

maxfun <- function(x){g[which(x %in% max(x))]}
optima<-unlist(lapply(site.optima,maxfun))

dev.off()

# -------------------------------------------------------------------
# Make plots
# -------------------------------------------------------------------
pdf(
  "~/Dropbox/clarkiaSeedBanks/products/figures/analysis/obs-pred-germ.pdf",width=6,height=4)

par(mar=c(4,5,4,1),mfrow=c(1,1))

plot(NA,NA,xlim=c(0,1),ylim=c(0,.4),
      cex.lab = 1.25, cex.axis = 1,
     xlab="",
     ylab="Observed germination fraction")

# Now, define a custom axis
points(x=optima,y=g1.hat,pch=21,col='black',bg='white',cex=1)
abline(a=0,b=1,lty='dashed')

mtext("Predicted germination fraction",
      side=1,line=2.1,adj=.5,col='black',cex=1.25)
# 
# g1.sum=apply(g1,2,quantile,probs=c(0.025,.25,.5,.75,.975))
# 
# par(fig=c(0,10,0,4.5)/10)
# par(new=T)
# plot(NA,NA,type='n',xlim=c(0,20),ylim=c(0,.6),
#      axes=FALSE,frame=FALSE,
#      xlab="",ylab="")
# y.pt = 1:20
# for(i in 1:20){
#   tmp<-g1.sum[,i]
#   segments(y0=tmp[1],y1=tmp[5],x0=y.pt[i])
#   segments(y0=tmp[2],y1=tmp[4],x0=y.pt[i],lwd=3)
#   points(y=tmp[3],x=y.pt[i],pch=21,bg='white')
# }
# axis(2,  seq(0,.6,by=.1), col.ticks = 1)
# axis(1, (1:20),
#      labels = (siteNames), las = 2, 
#      col = NA, col.ticks = 1, cex.axis = 1)
# mtext("Germination probability",
#       side=2,line=2.5,adj=.5,col='black',cex=1)
dev.off()


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

# calculate posterior mode
df.list = list()
for(i in 1:20){
  obj = rs
  zeroYears=lowFitnessYears[lowFitnessYears$siteIndex==i,]$yearIndex
  index=grep(paste0("\\[",i,","),colnames(obj))
  obj=obj[,index]
  obj[,zeroYears] = 0
  tmp.df=apply(obj,2,posterior.mode)
  df.list[[i]] = tmp.df
}

rs.hat=unlist(df.list)
index=expand.grid(1:20,1:15)

g = seq(0,1,length.out=101)
g.seq = g
fit = c()
g.mat = matrix(NA,nrow=1000,ncol=length(g))
g.sites = list()
yt = c()
n.iter=dim(rs)[1]

for( k in 1:20){
  rs.tmp = rs[,grep(paste0("\\[",k,","),colnames(rs))]
  rs.tmp = rs.tmp[sample(1:4500,1),]
  #  rs.tmp = rs.hat[grep(paste0("\\[",k,","),names(rs.hat))]
  
  yt = sample(rs.tmp,1000,replace=TRUE)
  for( i in 1:length(g)){
    fit<-fitness(g=g[i],s0.hat[k],s1.hat[k],s2.hat[k],s3.hat[k],rs=yt)
    logfit<-log(fit)
    g.mat[,i] <- fit
    #g.mat[,i] <- logfit
  }
  g.sites[[k]] <- g.mat
}

par(mfrow=c(1,1))
plot(g.seq,apply(g.mat,2,gm_mean))

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(k in 1:20){
  plot(g.seq,apply(g.sites[[k]],2,gm_mean),type='l')
  lines(g.seq,apply(g.sites[[k]],2,mean,na.rm=TRUE),type='l',lty='dotted')
}

par(mfrow=c(1,1))
plot(g.seq,apply(g.sites[[1]],2,gm_mean),type='n',ylim=c(0,20))
for(k in 1){
  lines(g.seq,apply(g.sites[[k]],2,gm_mean),type='l')
  lines(g.seq,apply(g.sites[[k]],2,mean),type='l',lty='dotted')
}



plot(g.seq,apply(g.mat,2,gm_mean))

gmean <- function(x){apply(x,2,gm_mean)}
site.optima<-lapply(g.sites,gmean)

maxfun <- function(x){g[which(x %in% max(x))]}
optima<-unlist(lapply(site.optima,maxfun))

dev.off()
pdf(
  "~/Dropbox/clarkiaSeedBanks/products/figures/analysis/obs-pred-germ-lowFitness.pdf",width=6,height=4)

par(mar=c(4,5,4,1),mfrow=c(1,1))

plot(NA,NA,xlim=c(0,1),ylim=c(0,.4),
     cex.lab = 1.25, cex.axis = 1,
     xlab="",
     ylab="Observed germination fraction")
# Now, define a custom axis
points(x=optima,y=g1.hat,pch=21,col='black',bg='white',cex=1)
abline(a=0,b=1,lty='dashed')

mtext("Predicted germination fraction",
      side=1,line=2,adj=.5,col='black',cex=1.25)

dev.off()

