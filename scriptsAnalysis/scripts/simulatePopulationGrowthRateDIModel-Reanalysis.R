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

set.seed(10)
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

# rs = rsMedians %>% tidyr::pivot_wider(values_from=c(med_sigma,med_fec,med_phi))
# rs<-split(rs,rs$site)

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/"
dataFiles <- paste0(directory,list.files(directory))

countSeedPerUndamagedFruit <- readRDS(dataFiles[[grep("countSeedPerUndamagedFruit",dataFiles)]])
# recover site indices and year indices
siteIndex <- data.frame(siteIndex=unique(countSeedPerUndamagedFruit$site3),site=1:20)
yearIndex <- data.frame(yearIndex=unique(countSeedPerUndamagedFruit$year3),
                        year=1:13) 
siteNames = unique(countSeedPerUndamagedFruit$site3)
index=expand.grid(1:20,1:13)

# -------------------------------------------------------------------

fitness <- function(g=g1,s0=s0,s1=s1,s2=s2,s3=s3,rs=rs){
  p1 = g*rs*s1*s0
  p2 = (1-g)*(s2*s3)
  return(as.numeric(p1+p2))
}

medianEst <- function(x){return(apply(x,2,median))}
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
#rs.hat = medianEst(rs)
rs.hat = modeEst(rs)

#rs.mean = meanEst(rs)
#dev.off()
#plot(rs.hat,rs.mean);abline(a=0,b=1)
#logfit<-c()

# -------------------------------------------------------------------
# calculate autocorrelation of per-capita reproductive success (mode)
# no significant autocorrelation, include in appendix

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

#(pop.harmonicMean,siteNames) %>% View


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

for( k in 1:20){
  pop.index=index[,1]==k
  rs.tmp = rs.hat[pop.index]
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
plot(g.seq,apply(g.mat,2,gm_mean))

# make plots of geometric mean fitness
# note that without reproductive failure fitness
# remains high at high values of g
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(k in 1:20){
  plot(g.seq,apply(g.sites[[k]],2,gm_mean),type='l')
#  lines(g.seq,apply(g.sites[[k]],2,mean,na.rm=TRUE),type='l',lty='dotted')
}

# plots of arithmetic mean fitness (log)
for(k in 1:20){
  plot(g.seq,log(apply(g.sites[[k]],2,mean)),type='l')
  #lines(g.seq,apply(g.sites[[k]],2,mean,na.rm=TRUE),type='l',lty='dotted')
}

# plots of standard deviation of pop growth rate
for(k in 1:20){
  plot(g.seq,apply(g.sites[[k]],2,sd),type='l')
}


# calculate optimal germination fraction
gmean <- function(x){apply(x,2,gm_mean)}
site.optima<-lapply(g.sites,gmean)

maxfun <- function(x){g[which(x %in% max(x))]}
optima<-unlist(lapply(site.optima,maxfun))

dev.off()
pdf(
  "~/Dropbox/clarkiaSeedBanks/products/figures/analysis/obs-pred-germ.pdf",width=6,height=8)

par(mar=c(4,5,2,1))
par(fig=c(0,10,4,10)/10)

plot(NA,NA,xlim=c(0,1),ylim=c(0,1),
      cex.lab = 1.5, cex.axis = 1,
     xlab="",
     ylab="Observed germination fraction")
# Now, define a custom axis
points(x=optima,y=g1.hat,pch=21,col='black',bg='white',cex=.75)
abline(a=0,b=1,lty='dashed')

mtext("Predicted germination fraction",
      side=1,line=2.1,adj=.5,col='black',cex=1)

g1.sum=apply(g1,2,quantile,probs=c(0.025,.25,.5,.75,.975))

par(fig=c(0,10,0,4.5)/10)
par(new=T)
plot(NA,NA,type='n',xlim=c(0,20),ylim=c(0,.6),
     axes=FALSE,frame=FALSE,
     xlab="",ylab="")
y.pt = 1:20
for(i in 1:20){
  tmp<-g1.sum[,i]
  segments(y0=tmp[1],y1=tmp[5],x0=y.pt[i])
  segments(y0=tmp[2],y1=tmp[4],x0=y.pt[i],lwd=3)
  points(y=tmp[3],x=y.pt[i],pch=21,bg='white')
}
axis(2,  seq(0,.6,by=.1), col.ticks = 1)
axis(1, (1:20),
     labels = (siteNames), las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)
mtext("Germination probability",
      side=2,line=2.5,adj=.5,col='black',cex=1)
dev.off()



### Repeat analysis with lower fitness data
directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/"
dataFiles <- paste0(directory,list.files(directory))
countSeedPerUndamagedFruit <- readRDS(dataFiles[[grep("countSeedPerUndamagedFruit",dataFiles)]])
# recover site indices and year indices
siteIndex <- data.frame(site=unique(countSeedPerUndamagedFruit$site3),siteIndex=1:20)
yearIndex <- data.frame(year=as.numeric(unique(countSeedPerUndamagedFruit$year3)),
                        yearIndex=1:13) 
siteNames = unique(countSeedPerUndamagedFruit$site3)

rs <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/rsPosterior.RDS")
lowFitnessYears <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/output/lowFitnessYearsPlots.RDS")
lowFitnessYears=lowFitnessYears %>% dplyr::left_join(yearIndex) %>% dplyr::left_join(siteIndex)

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
length(rs.hat)

index=expand.grid(1:20,1:13)

g = seq(0,.99,length.out=101)
g.seq = g
fit = c()
g.mat = matrix(NA,nrow=1000,ncol=length(g))
g.sites = list()
yt = c()
n.iter=dim(rs)[1]

for( k in 1:20){
  pop.index=index[,1]==k
  rs.tmp = rs.hat[pop.index]
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
  "~/Dropbox/clarkiaSeedBanks/products/figures/analysis/obs-pred-germ-lowFitness.pdf",width=6,height=8)

par(mar=c(4,4,2,1))
par(fig=c(0,10,4,10)/10)

plot(NA,NA,xlim=c(0,1),ylim=c(0,1),
     cex.lab = 1, cex.axis = 1,
     xlab="",
     ylab="Observed germination fraction")
# Now, define a custom axis
points(x=optima,y=g1.hat,pch=21,col='black',bg='white',cex=.75)
abline(a=0,b=1,lty='dashed')

mtext("Predicted germination fraction",
      side=1,line=2,adj=.5,col='black',cex=1)

g1.sum=apply(g1,2,quantile,probs=c(0.025,.25,.5,.75,.975))

par(fig=c(0,10,0,4.5)/10)
par(new=T)
plot(NA,NA,type='n',xlim=c(0,20),ylim=c(0,.6),
     axes=FALSE,frame=FALSE,
     xlab="",ylab="")
y.pt = 1:20
for(i in 1:20){
  tmp<-g1.sum[,i]
  segments(y0=tmp[1],y1=tmp[5],x0=y.pt[i])
  segments(y0=tmp[2],y1=tmp[4],x0=y.pt[i],lwd=3)
  points(y=tmp[3],x=y.pt[i],pch=21,bg='white')
}
axis(2,  seq(0,.6,by=.1), col.ticks = 1)
axis(1, (1:20),
     labels = (siteNames), las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)
mtext("Germination probability",
      side=2,line=2.5,adj=.5,col='black',cex=1)
dev.off()


# 
# g = seq(0,1,by=.1)
# g.seq = g
# fit = c()
# g.mat = matrix(NA,nrow=1000,ncol=length(g))
# g.sites = list()
# yt = c()
# n.iter=dim(rs)[1]
# 
# for( k in 1:20){
#   # get index for population
#   pop.index=index[,1]==k
#   # get vector of years of fitness 
#   yt = sample(1:13,1000,replace=TRUE)
#   # get sample from posterior
#   sample.ind=sample(1:n.iter,1000,replace=TRUE)
#   # replace year with draw from that year
#   for(i in 1:length(yt)){
#     yt[i]=(rs[,yt[i]])[sample.ind[i]]
#   }
#   s0.ind=sample(1:n.iter,1000,replace=TRUE)
#   s1.ind=sample(1:n.iter,1000,replace=TRUE)
#   s2.ind=sample(1:n.iter,1000,replace=TRUE)
#   s3.ind=sample(1:n.iter,1000,replace=TRUE)
#   for( i in 1:length(g)){
#     fit<-fitness(g=g[i],s0[s0.ind,k],s1[s1.ind,k],s2[s2.ind,k],s3[s3.ind,k],rs=yt)
#     logfit<-log(fit)
#     g.mat[,i] <- fit
#     #g.mat[,i] <- logfit
#   }
#   g.sites[[k]] <- g.mat
# }
# 
# # par(mfrow=c(1,1))
# # plot(g.seq,apply(g.mat,2,gm_mean))
# 
# par(mfrow = c(4,5),
#     oma = c(5,4,0,0) + 0.1,
#     mar = c(0,0,1,1) + 0.1)
# for(k in 1:20){
#   plot(g.seq,apply(g.sites[[k]],2,gm_mean),type='l')
#   lines(g.seq,apply(g.sites[[k]],2,mean,na.rm=TRUE),type='l',lty='dotted')
# }
# 
# par(mfrow=c(1,1))
# plot(g.seq,apply(g.sites[[1]],2,gm_mean),type='n',ylim=c(0,20))
# for(k in 1){
#   lines(g.seq,apply(g.sites[[k]],2,gm_mean),type='l')
#   lines(g.seq,apply(g.sites[[k]],2,mean),type='l',lty='dotted')
# }
# 
# 
# 
# plot(g.seq,apply(g.mat,2,gm_mean))
# 
# gmean <- function(x){apply(x,2,gm_mean)}
# site.optima<-lapply(g.sites,gmean)
# 
# maxfun <- function(x){g[which(x %in% max(x))]}
# optima<-unlist(lapply(site.optima,maxfun))
# 
# #dev.off()
# # pdf(
# #   "~/Dropbox/clarkiaSeedBanks/products/figures/analysis/obs-pred-germ-sample.pdf",width=4,height=4)
# 
# plot(optima,g1.hat,xlim=c(0,1),ylim=c(0,1),pch=16,
#      cex.lab = 1.5, cex.axis = 1.5,
#      xlab="Predicted germination fraction",
#      ylab="Observed germination fraction")
# # Now, define a custom axis
# abline(a=0,b=1,lty='dashed')
# dev.off()
# 


## Cohen 1966

#dev.off()
# 
# 
# names<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>%
#   dplyr::select(site)
# df <- data.frame(names,optima,g1)
# 
# library(ggrepel)
# 
# g1 <- ggplot(df,aes(x=optima,y=g1,label=site)) +
#   geom_abline(intercept=0,slope=1,alpha=.5) +
#   geom_point() +
#   geom_text_repel(size=3,color="black") +
#   theme_bw() + xlim(c(0,1)) + ylim(c(0,1)) +
#   xlab("Predicted germination fraction") +
#   ylab("Observed germination fraction")
# 
# g1
# # 
# # ggsave(filename="~/Dropbox/clarkiaSeedBanks/products/figures/obs_pred_germ-labeled.pdf",
# #        plot=g1,width=6,height=6)
# # 
# # 
# # 
# # optimaCalc<-function(s1,s2,s3,rs){
# #   
# #   g = seq(0,1,by=.01)
# #   fit = c()
# #   g.mat = matrix(NA,nrow=100,ncol=length(g))
# #   g.sites = list()
# #   yt=c()
# #   
# #   for( k in 1:20){
# #     yt = rs[[k]][sample(7,100,replace=TRUE), ]
# # #    y = rs[[k]][sample(30000,1000), ]
# #   #  year = sample(1:(dim(y)[2]),1000,replace=TRUE)
# #     # for(i in 1:1000){
# #     #   yt[i] <- y[i,year[i]]
# #     # }
# #     for( i in 1:length(g)){
# #       for( j in 1:100){
# #         fit[j]<-fitness(g=g[i],s1[k],s2[k],s3[k],rs=sample(yt[i,],1))
# #         #logfit[i]<-log(fit)
# #       }
# #       g.mat[,i] <- fit
# #     }
# #     g.sites[[k]] <- g.mat
# #   }
# #   
# #   gmean <- function(x){apply(x,2,gm_mean)}
# #   site.optima<-lapply(g.sites,gmean)
# #   
# #   maxfun <- function(x){g[which(x %in% max(x))]}
# #   optima<-unlist(lapply(site.optima,maxfun))
# # }
# # 
# # # uncomment to re-estimate mean.optima
# # # need to save mean.optima to create figure
# # trials <- replicate(
# #   50,
# #   optimaCalc(s1 = s1, s2 = s2, s3 = s3, rs= rs),
# #   simplify = FALSE
# # )
# # # 
# # # set.seed(10)
# # # mean.optima = apply(matrix(unlist(trials),ncol=17,byrow=TRUE),2,mean)
# # # 
# # # optimalG<-matrix(unlist(trials),ncol=17,byrow=TRUE)
# # #save(optimalG,file="~/Dropbox/modelsF2019/output/optimalG")
# # 
# # # load("~/Dropbox/modelsF2019/output/optimalG")
# # 
# # o.mu <- apply(optimalG,2,mean)
# # o.se <- apply(optimalG,2,sd)/sqrt(50)
# # 
# # pdf(
# #   "~/Dropbox/clarkiaSeedBanks/products/figures/obs_pred_germ.pdf",
# #   onefile=TRUE,
# #   paper="USr",
# #   height = 7.5, width = 10)
# # 
# # plot(mean.optima,g1,xlim=c(0,1),ylim=c(0,1),pch=16,
# #      cex.lab = 1.5, cex.axis = 1.5,
# #      xlab="Predicted germination fraction",
# #      ylab="Observed germination fraction")
# # # Now, define a custom axis
# # abline(a=0,b=1,lty='dashed')
# # ## Cohen 1966
# # segments(y0=g1,x0=o.mu-o.se,x1=o.mu+o.se)
# # 
# # dev.off()
# # 
# # names<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
# #   dplyr::select(site) 
# # df <- data.frame(names,optima,g1)
# # 
# # library(ggrepel)
# # 
# # g1 <- ggplot(df,aes(x=optima,y=g1,label=site)) +
# #   geom_abline(intercept=0,slope=1,alpha=.5) +
# #   geom_point() +
# #   geom_text_repel(size=3,color="black") +
# #   theme_bw() + xlim(c(0,1)) + ylim(c(0,1)) +
# #   xlab("Predicted germination fraction") +
# #   ylab("Observed germination fraction")
# # 
# # ggsave(filename="~/Dropbox/clarkiaSeedBanks/products/figures/obs_pred_germ-labeled.pdf",
# #        plot=g1,width=6,height=6)
