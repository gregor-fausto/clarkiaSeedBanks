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

posterior.mode = function(x){
  x.max=max(x)
  x.min=min(x)
  dres <- density( x ,from = x.min, to = x.max,adjust=2.5)
  modeParam <- dres$x[which.max(dres$y)]
  return(modeParam)
}

# -------------------------------------------------------------------
# Loading posterior
# -------------------------------------------------------------------
directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/"
sigma = readRDS(file="~/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/sigma-analysis.RDS")
fec = readRDS(file="~/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/tfe-analysis.RDS")
phi = readRDS(file="~/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/phi-analysis.RDS")

position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>%
  dplyr::select(site,easting,dominant.surface.rock.type,elevation) %>%
  dplyr::mutate(easting=easting/1000)
siteNames = position$site

# -------------------------------------------------------------------
# Randomly sample a site
# -------------------------------------------------------------------
siteNumber=sample(1:20,1)

sigma = sigma[,grep(paste0("\\[",siteNumber,","),colnames(sigma))]
fec = fec[,grep(paste0("\\[",siteNumber,","),colnames(fec))]
phi = phi[,grep(paste0("\\[",siteNumber,","),colnames(phi))]

# -------------------------------------------------------------------
# Make plots for each vital rate
# ------------------------------------------------------------------- 
f=function(x){
  tmp=density(x,from=min(x),to=max(x),adjust=2.5)
  df=cbind(tmp$x,tmp$y)
  df=rbind(c(tmp$x[1],0),df)
  df=rbind(df,c(tmp$x[length(tmp$x)],0))
  return(df)
}

f.boxplot=function(x,y,width=.2){
  d=boxplot(y,plot=FALSE)
  segments(x0=x,y0=d$stats[2],y1=d$stats[4],lty='solid',lwd=2)
  points(x=x,y=d$stats[3],pch=21,bg='white',cex=1)
}



sigma.sum=apply(sigma,2,quantile,c(0.025,.25,.5,.75,.975))
fec.sum=apply(fec,2,quantile,c(0.025,.25,.5,.75,.975))
phi.sum=apply(phi,2,quantile,c(0.025,.25,.5,.75,.975))

rs=sigma*fec*phi
rs.sum=apply(rs,2,quantile,c(0.025,.25,.5,.75,.975))

par(mfrow=c(4,1),mar=c(0,.25,.5,0),
    oma=c(4,4,1,1))
# sigma
plot(NA,NA,
     type='b',pch=19,
     xlim=c(2006,2020),ylim=c(0,1),
     frame=TRUE,xaxt='none',
     ylab="Probability of seedling survival to fruiting",
     xlab="")
segments(x0=2006:2020,y0=sigma.sum[1,],y1=sigma.sum[5,])
segments(x0=2006:2020,y0=sigma.sum[2,],y1=sigma.sum[4,],lwd=2.5)
points(2006:2020,sigma.sum[3,],pch=21,cex=1,bg='white',type='b')
# axis(2, seq(0,1,by=.2),
#      labels = seq(0,1,by=.2), las = 2, 
#      col = NA, col.ticks = 1, cex.axis = 1)
mtext(paste0("Randomly chosen population: ", siteNames[siteNumber]),adj=0)
mtext("Probability(seedling survival to fruiting)", side=2,line = 2.2,cex=.5)

#fec
plot(NA,NA,
     type='b',pch=19,
     xlim=c(2006,2020),ylim=c(0,max(fec.sum[5,])),
     frame=TRUE,xaxt='none',
     ylab="Fruits per plant",
     xlab="")
segments(x0=2006:2020,y0=fec.sum[1,],y1=fec.sum[5,])
segments(x0=2006:2020,y0=fec.sum[2,],y1=fec.sum[4,],lwd=2.5)
points(2006:2020,fec.sum[3,],pch=21,cex=1,bg='white',type='b')
mtext("Fruits per plant", side=2,line = 2.2,cex=.5)


#phi
plot(NA,NA,
     type='b',pch=19,
     xlim=c(2006,2020),ylim=c(min(phi.sum[1,]),max(phi.sum[5,])),
     frame=TRUE,xaxt='none',
     ylab="Seeds per fruit",
     xlab="")
segments(x0=2006:2020,y0=phi.sum[1,],y1=phi.sum[5,])
segments(x0=2006:2020,y0=phi.sum[2,],y1=phi.sum[4,],lwd=2.5)
points(2006:2020,phi.sum[3,],pch=21,cex=1,bg='white',type='b')
mtext("Seeds per fruit", side=2,line = 2.2,cex=.5)


#rs
plot(NA,NA,
     type='b',pch=19,
     xlim=c(2006,2020),ylim=c(min(rs.sum[1,]),max(rs.sum[5,])),
     frame=TRUE,
     ylab="Per-capita reproductive success",
     xlab="")
segments(x0=2006:2020,y0=rs.sum[1,],y1=rs.sum[5,])
segments(x0=2006:2020,y0=rs.sum[2,],y1=rs.sum[4,],lwd=2.5)
points(2006:2020,rs.sum[3,],pch=21,cex=1,bg='white',type='b')

mtext("Per-capita reproductive success", side=2,line = 2.2,cex=.5)


directory = "/Users/Gregor/Dropbox/dataLibrary/postProcessingData-2021/"
dataFiles <- paste0(directory,list.files(directory))
censusSeedlingsFruitingPlants = readRDS(dataFiles[[grep("censusSeedlingsFruitingPlants.RDS",dataFiles)]])

sigma.all = readRDS(file="~/Dropbox/dataLibrary/mcmcSamplesThinned/seedlingSurvivalSamples.RDS")
mu_pop=MCMCchains(sigma.all,params='mu0')[,siteNumber]
sigma_pop = MCMCchains(sigma.all,params="sigma0")[,siteNumber]
dist=rnorm(4500,mean=mu_pop,sd=sigma_pop)

siteNames =censusSeedlingsFruitingPlants$site %>% unique
df=censusSeedlingsFruitingPlants[censusSeedlingsFruitingPlants$site==siteNames[siteNumber],]
years=2006:2020
# seedling survival to fruiting
dev.off()
library(plotrix)

pdf("~/Dropbox/clarkiaSeedBanks/products/manuscript/figures-overview/interannual-model.pdf",width=10,height=7)
par(mfrow=c(1,1))
plot(x=NA,NA,
     type='n',
     xlim=c(2006,2022.5),ylim=c(0,1),
     #axes=FALSE,
     frame=FALSE,xaxt='n',yaxt='n',
     ylab="Probability of seedling survival to fruiting",
     xlab="Year",cex=3,cex.axis=3)
year=2006:2020

# put a break at the default axis and position


polygon(x=c(2007,2007,2008,2008)-.25,y=c(-1,2,2,-1),col='gray98',border=0)
# polygon(x=c(2009,2009,2010,2010)-.25,y=c(-1,2,2,-1),col='gray98',border=0)
# polygon(x=c(2011,2011,2012,2012)-.25,y=c(-1,2,2,-1),col='gray98',border=0)
# polygon(x=c(2013,2013,2014,2014)-.25,y=c(-1,2,2,-1),col='gray98',border=0)
# polygon(x=c(2015,2015,2016,2016)-.25,y=c(-1,2,2,-1),col='gray98',border=0)
# polygon(x=c(2017,2017,2018,2018)-.25,y=c(-1,2,2,-1),col='gray98',border=0)
# polygon(x=c(2019,2019,2020,2020)-.25,y=c(-1,2,2,-1),col='gray98',border=0)

abline(v=2005.3)
for(i in c(1,2,3)){
  df.tmp=sigma[,i]
  dat.tmp = df[df$year==years[i],] %>%
    dplyr::mutate(p = fruitplNumber/seedlingNumber) %>%
    dplyr::filter(!is.na(p))
  upper.limit=max(f(df.tmp)[,2])*1.5
  polygon(x=year[i]+f(df.tmp)[,2]/upper.limit,y=f(df.tmp)[,1],col='orange',border='orange')
  f.boxplot(year[i],df.tmp)
  n.obs = length(dat.tmp$seedlingNumber)
  size = dat.tmp$seedlingNumber/max(df$seedlingNumber,na.rm=TRUE)
   points(rep(i+2005-0.1,n.obs)+rnorm(n.obs,0,.025),dat.tmp$p, pch = 16, cex = .5)
}

for(i in 14:15){

df.tmp=sigma[,i]
dat.tmp = df[df$year==years[i],] %>%
  dplyr::mutate(p = fruitplNumber/seedlingNumber) %>%
  dplyr::filter(!is.na(p))
upper.limit=max(f(df.tmp)[,2])*1.5
polygon(x=year[i]-10+f(df.tmp)[,2]/upper.limit,y=f(df.tmp)[,1],col='orange',border='orange')
f.boxplot(year[i]-10,df.tmp)
n.obs = length(dat.tmp$seedlingNumber)
size = dat.tmp$seedlingNumber/max(df$seedlingNumber,na.rm=TRUE)
points(rep(i-10+2005-0.1,n.obs)+rnorm(n.obs,0,.025),dat.tmp$p, pch = 16, cex = .5)

}

axis(2, seq(0,1,by=.2),
     labels = seq(0,1,by=.2), las = 2,
     col = NA, col.ticks = 1, cex.axis = 1)
axis(1, c(seq(2006,2020,by=2)),
     labels = c(seq(2006,2020,by=2)), las = 1,
     col = NA, col.ticks = 1, cex.axis = 1)

rect(2021,-1,2024,2,col='gray80',border='gray80')
  df.tmp=boot::inv.logit(dist)
  upper.limit=max(f(df.tmp)[,2])*1
  polygon(x=2021.5+f(df.tmp)[,2]/upper.limit,y=f(df.tmp)[,1],col='purple',border='purple')
  f.boxplot(2021.5,df.tmp)

#abline(h=c(0,1))
box()
axis.break(1,2009,style='zigzag')

dev.off()


# RS
par(mfrow=c(1,1))
plot(x=NA,NA,
     type='n',
     xlim=c(2005,2021),ylim=c(0,500),
     axes=FALSE,frame=FALSE,
     ylab="Probability of seedling survival to fruiting",
     xlab="")
year=2006:2020
for(i in 1:15){
  df.tmp=rs[,i]
  upper.limit=max(f(df.tmp)[,2])*1.2
  polygon(x=year[i]-f(df.tmp)[,2]/upper.limit,y=f(df.tmp)[,1],col='gray80',border=0)
  f.boxplot(year[i]-.025,df.tmp)
 # points(rep(i+.25,length(df.tmp$rs))+rnorm(length(df.tmp$rs),0,.025),df.tmp$rs, pch = 16, cex = .5)
}
axis(2, seq(0,1,by=.2),
     labels = seq(0,1,by=.2), las = 2,
     col = NA, col.ticks = 1, cex.axis = 1)
axis(1, (1:15),
     labels = 2006:2020, las = 2,
     col = NA, col.ticks = 1, cex.axis = 1)



index=rev(order(apply(rs,2,posterior.mode)))
vec=c()
for(i in 1:15){
  df.tmp=rs[,index[i]]
  vec[i]=max(f(df.tmp)[,1])
}

par(mfrow=c(1,1))

plot(x=NA,NA,
     type='n',
     xlim=c(0,max(vec)),ylim=c(0,15),
     #axes=FALSE,frame=FALSE,
     yaxt='none',
     ylab="Probability of seedling survival to fruiting",
     xlab="")
abline(v=0,col='black',lty='dotted')
year=2006:2020
for(i in 1:15){
  df.tmp=rs[,index[i]]
  upper.limit=max(f(df.tmp)[,2])*1.2
  polygon(y=i+(f(df.tmp)[,2])/upper.limit-1,x=f(df.tmp)[,1],col=rgb(.9, .9, .9,0.5),border="black")
  val=f(df.tmp)[,1][max((f(df.tmp)[,2])/upper.limit)==(f(df.tmp)[,2])/upper.limit]
  
  # mode
 # segments(x0=posterior.mode(df.tmp),y0=i-1,y1=max(i+(f(df.tmp)[,2])/upper.limit-1),lwd=1,lty='dotted')
  points(y=i-1,x=posterior.mode(df.tmp),col='black',pch=16,cex=.5)
  
  # mean
 # segments(x0=mean(df.tmp),y0=i-1,y1=1,lwd=1,lty='dotted',col='orange')
  points(y=i-1,x=mean(df.tmp),col='orange',pch=16,cex=.5)
  
  # median
#  segments(x0=median(df.tmp),y0=i-1,y1=max(i+(f(df.tmp)[,2])/upper.limit-1),lwd=1,lty='dotted',col='blue')
  points(y=i-1,x=median(df.tmp),col='blue',pch=16,cex=.5)
  }

rs.mode=apply(rs,2,posterior.mode)
sample.seq=sample(rs.mode,1000,replace=TRUE)
hist(sample.seq,col='gray90',border='white',freq=FALSE)

plot(rs.mode,pch=16,type='n',ylim=c(0,500))
iter.seq=sample(1:dim(rs)[1],50,replace=TRUE)
for(i in 1:50){
  lines(rs[i,],type='l',col=rgb(.9, .9, .9,0.5))
}
lines(rs.mode,col='black',lwd=1.5)
lines(rs[iter.seq[1],],col='purple',lwd=1.5)
lines(rs[iter.seq[2],],col='orange',lwd=1.5)


# plot(x=NA,NA,
#      type='n',
#      xlim=c(0,25),ylim=c(0,15),
#      #axes=FALSE,frame=FALSE,
#      yaxt='none',
#      ylab="Probability of seedling survival to fruiting",
#      xlab="")
# abline(v=0,col='black',lty='dotted')
# year=2006:2020
# for(i in 1:15){
#   df.tmp=rs[,index[i]]
#   upper.limit=max(f(df.tmp)[,2])*1.2
#   polygon(y=i+(f(df.tmp)[,2])/upper.limit-1,x=f(df.tmp)[,1],col=rgb(.9, .9, .9,0.5),border="black")
#   val=f(df.tmp)[,1][max((f(df.tmp)[,2])/upper.limit)==(f(df.tmp)[,2])/upper.limit]
#   
#   # mode
#   # segments(x0=posterior.mode(df.tmp),y0=i-1,y1=max(i+(f(df.tmp)[,2])/upper.limit-1),lwd=1,lty='dotted')
#   points(y=i-1,x=posterior.mode(df.tmp),col='black',pch=16,cex=.5)
#   
#   # mean
#   # segments(x0=mean(df.tmp),y0=i-1,y1=1,lwd=1,lty='dotted',col='orange')
#   points(y=i-1,x=mean(df.tmp),col='orange',pch=16,cex=.5)
#   
#   # median
#   #  segments(x0=median(df.tmp),y0=i-1,y1=max(i+(f(df.tmp)[,2])/upper.limit-1),lwd=1,lty='dotted',col='blue')
#   points(y=i-1,x=median(df.tmp),col='blue',pch=16,cex=.5)
# }


