#################################################################################
################################################################################
################################################################################
# Code for figures to compare the following modeling approaches for the seedling survivorship data
# 
# Scripts by Gregor Siegmund
# fausto.siegmund@gmail.com
# last updated 03-07-2021
#################################################################################
#################################################################################
#################################################################################
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)

library(MCMCvis)
library(tidybayes)
library(tidyverse)
library(magrittr)
library(bayesplot)
library(rethinking)

directory = "/Users/Gregor/Dropbox/dataLibrary/mcmcSamplesThinned/"
mcmcSampleDirectory <- paste0(directory,list.files(directory))
mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("seedlingSurvivalSamples.rds",mcmcSampleDirectory)]])

dirPars = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/parameterSummary/"
dir.create(file.path(dirPars), showWarnings = FALSE)

################################################################################
# Data directory
#################################################################################

data <- readRDS(mcmcSampleDirectory[[grep("seedlingSurvivalData.rds",mcmcSampleDirectory)]])

# -------------------------------------------------------------------
# Get site names and position
# -------------------------------------------------------------------

directory = "/Users/Gregor/Dropbox/dataLibrary/postProcessingData-2021/"
dataFiles <- paste0(directory,list.files(directory))
censusSeedlingsFruitingPlants = readRDS(dataFiles[[grep("censusSeedlingsFruitingPlants.RDS",dataFiles)]])

siteNames = unique(censusSeedlingsFruitingPlants$site)
years = as.numeric(unique(censusSeedlingsFruitingPlants$year))
yearNames = years[order(years)]

position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,dominant.surface.rock.type) %>%
  dplyr::mutate(easting=easting/1000)

################################################################################
# Summarize seedlings survival parameters in CSV
#################################################################################

siteIndex <- data.frame(siteIndex=unique(censusSeedlingsFruitingPlants$site),site=1:20)
yearIndex <- data.frame(yearIndex=as.numeric(unique(censusSeedlingsFruitingPlants$year)),
                        year=unique(unique(data$year)))

summary.fun = function(x){
  quant=quantile(x,c(.025,.975))
  mu = mean(x)
  med = median(x)
  hpdi.interval = HPDI(x, prob = .95)
  return(c(mu,quant,med,hpdi.interval))
}

parm.mu=MCMCchains(mcmcSamples,params='mu')
parm.prob=apply(parm.mu,2,boot::inv.logit)
parm.prob.sum = apply(parm.prob,2,summary.fun)

df.list = list()
for(i in 1:20){
  obj = parm.prob.sum
  index=grep(paste0("\\[",i,","),colnames(obj))
  tmp = signif(parm.prob.sum[,index],3)
  tmp.df=data.frame(site=siteNames[i],year=yearNames,(t(tmp)))
  names(tmp.df) = c("site","year","mu","ci.lo95","ci.hi95","med","hpdi.lo95","hpdi.hi95")
  rownames(tmp.df) = NULL
  df.list[[i]] = tmp.df
}

summary.df=do.call(rbind,df.list)
# 
# sigma_py<-mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(mu[site,year]) %>%
#   dplyr::mutate(p_1 = boot::inv.logit(mu)) %>%
#   dplyr::group_by(site,year) %>%
#   dplyr::summarise(mu = mean(p_1),
#                    ci.lo95 = quantile(p_1,probs=0.025), 
#                    ci.hi95 = quantile(p_1,probs=0.975),
#                    med = median(p_1), 
#                    hpdi.lo95 = HPDI(p_1,.95)[1], 
#                    hpdi.hi95 = HPDI(p_1,.95)[2]
#   )
# 
# sigma_py<-sigma_py %>%
#   dplyr::left_join(siteIndex,by="site") %>%
#   dplyr::left_join(yearIndex,by="year") %>%
#   dplyr::ungroup() %>%
#   dplyr::select(-c(site,year)) %>%
#   dplyr::rename(site = siteIndex) %>%
#   dplyr::rename(year = yearIndex) %>%
#   dplyr::select(site,year,mu,ci.lo95,ci.hi95,med,hpdi.lo95,hpdi.hi95)

# saveRDS(sigma_py,file=paste0(dirPars,"sigmaSummary.RDS"))
sigma_py = summary.df 
write.csv(sigma_py,file=paste0(dirPars,"sigmaSummary.csv"))


################################################################################
# Calculate summary objects
#################################################################################


summary.fun = function(x){
  med = median(x)
  hpdi.interval = HPDI(x, prob = .89)
  hpdi.interval2 = HPDI(x, prob = .5)
  return(c(med,hpdi.interval,hpdi.interval2))
}

parm.mu0=MCMCchains(mcmcSamples,params='mu0')
parm.prob0=apply(parm.mu0,2,boot::inv.logit)
parm.prob0.sum = apply(parm.prob0,2,summary.fun)

df.list = list()
for(i in 1:20){
  obj = parm.prob0.sum
  index=grep(paste0("\\[",i,"\\]"),colnames(obj))
  tmp = signif(parm.prob0.sum[,index],3)
  tmp.df=data.frame(site=siteNames[i],(t(tmp)))
  names(tmp.df) = c("site","med","ci.lo","ci.hi","ci.lo2","ci.hi2")
  rownames(tmp.df) = NULL
  df.list[[i]] = tmp.df
}

summary.pop.df=do.call(rbind,df.list)


parm.prob=apply(parm.mu,2,boot::inv.logit)
parm.prob.sum = apply(parm.prob,2,summary.fun)

df.list = list()
for(i in 1:20){
  obj = parm.prob.sum
  index=grep(paste0("\\[",i,","),colnames(obj))
  tmp = signif(parm.prob.sum[,index],3)
  tmp.df=data.frame(site=siteNames[i],year=yearNames,(t(tmp)))
  names(tmp.df) = c("site","year","med","ci.lo","ci.hi","ci.lo2","ci.hi2")
  rownames(tmp.df) = NULL
  df.list[[i]] = tmp.df
}

summary.df=do.call(rbind,df.list)


# rename
sigma_p = summary.pop.df %>% dplyr::select(site,med)
sigma_py = summary.df 


################################################################################
# Make summary plots
#################################################################################


tmp = sigma_py %>%
  dplyr::rename(lambda.med = med) %>%
  dplyr::left_join(sigma_p,by = "site") 


interannualSigmaDF =  tmp %>%
  dplyr::left_join(position,by='site') %>%
  dplyr::arrange(easting) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(desc(lambda.med),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  mutate(id = row_number())

interannualSigma <- interannualSigmaDF %>%
  ggplot(aes(x = id , y = lambda.med)) + 
  geom_hline(aes(yintercept=med),linetype='dotted') +
  geom_point(aes(color=year)) +
  geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +
  #coord_flip() +
  facet_grid(. ~ site, scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  xlab('') +
  ylab("Probability of seedling survival to fruiting") +
  ylim(c(0,1)) +
  labs(color="Year") 

sigma_p = summary.pop.df

spatialSigma <- sigma_p %>%
  dplyr::left_join(position,by="site") %>%
  ggplot(aes(x = easting , y = med)) + 
  geom_point() +
  geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  ylab("Probability of seedling survival to fruiting") +
  xlab("Easting (km)") +
  ylim(c(0,1))

dirFigures = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/"

ggsave(filename=paste0(dirFigures,"interannualSigma.pdf"),
       plot=interannualSigma,width=12,height=6)

ggsave(filename=paste0(dirFigures,"spatialSigma.pdf"),
       plot=spatialSigma,width=6,height=12)


################################################################################
# Make summary plots
#################################################################################

f.param = function(x.v,parm="g1"){
  
  x.sum=apply(boot::inv.logit(x.v),2,quantile,probs=c(0.025,.25,.5,.75,.975))
  
  plot(NA,NA,type='n',xlim=c(0,1),ylim=c(0,20),
       axes=FALSE,frame=FALSE,
       xlab="",ylab="")
  
  position.index=order(position$easting)
  
  y.pt = 20:1
  for(i in 20:1){
    tmp<-x.sum[,position.index[i]]
    segments(x0=tmp[1],x1=tmp[5],y0=y.pt[i])
    segments(x0=tmp[2],x1=tmp[4],y0=y.pt[i],lwd=3)
    points(x=tmp[3],y=y.pt[i],pch=21,bg='white')
  }
  
  axis(1,  seq(0,1,by=.2), col.ticks = 1)
  axis(2, (1:20),
       labels = (siteNames)[position.index], las = 1, 
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

par(mfrow=c(1,2),mar=c(0,2,2,0),
    oma=c(4,1,1,1))
f.param(parm.mu0,parm="sigma")


### INDIVIDUAL PLOTS
dev.off()
par(mfrow=c(1,1))
x.v = parm.mu0
x.sum=apply(boot::inv.logit(x.v),2,quantile,probs=c(0.025,.25,.5,.75,.975))

plot(NA,NA,type='n',xlim=c(0,20),ylim=c(0,1),
     axes=FALSE,frame=FALSE,
     xlab="",ylab="")

position.index=order(position$easting)

y.pt = 20:1
for(i in 20:1){
  tmp<-x.sum[,position.index[i]]
  segments(y0=tmp[1],y1=tmp[5],x0=y.pt[i])
  segments(y0=tmp[2],y1=tmp[4],x0=y.pt[i],lwd=3)
  points(y=tmp[3],x=y.pt[i],pch=21,bg='white')
}

axis(2,  seq(0,1,by=.2), col.ticks = 1, las =1)
axis(1, (1:20),
     labels = (siteNames)[position.index], las = 3, 
     col = NA, col.ticks = 1, cex.axis = 1)
mtext("Probability of seedling survival to fruiting",
      side=2,line=2.5,adj=.5,col='black',cex=1)

x.sum.df<-data.frame(t(x.sum),position)
names(x.sum.df)[1:5] = c("ci.lolo","ci.lo","ci.med","ci.hi","ci.hihi")

plot(NA,NA,type='n',ylim=seq(0,1),xlim=c(340,375),
     axes=FALSE,frame=FALSE,
     xlab="",ylab="")
#abline(h=0,col='gray')

segments(x0=x.sum.df$easting,y0=x.sum.df$ci.lolo,y1=x.sum.df$ci.hihi,lwd=1)
segments(x0=x.sum.df$easting,y0=x.sum.df$ci.lo,y1=x.sum.df$ci.hi,lwd=3)
points(x=x.sum.df$easting,y=x.sum.df$ci.med,pch=21,bg='white')

axis(2, seq(0,1,by=.2), col.ticks = 1, las = 1)
axis(1, seq(340,375,by=5),
     labels = seq(340,375,by=5), las = 1,
     col.ticks = 1, cex.axis = 1)
mtext("Probability of seedling survival to fruiting",
      side=2,line=2.5,adj=.5,col='black',cex=1)
mtext("Easting (km)",
      side=1,line=2.5,adj=.5,col='black',cex=1)



position.index = order(position$easting)






sigma_py
sigma_p = summary.pop.df

dev.off()


par(mfrow=c(4,5),mar=c(0,.25,.25,0),
    oma=c(4,4,1,1))

time.sample = 1:15+2005

for(i in 1:20){

  tmp=sigma_py[sigma_py$site==siteNames[position.index[i]],]
  tmp.pop=sigma_p[sigma_p$site==siteNames[position.index[i]],]
  
  tmp.vec=c()
  for(j in 1:15){
    tmp.vec[j]=sum(data$seedlingNumber[data$site==position.index[i]&data$year==j])
  }

  plot(NA,NA,
       ylim=c(0,1),pch=16,xlim=c(2006,2020),
       ylab='',xlab='',xaxt='n',yaxt='n')
  
  index=tmp.vec==0
  if (sum(index)>0) {for(j in 1:length(time.sample[index])){
    ts=time.sample[index]
    polygon(x=c(ts[j]-.5,ts[j]-.5,
                ts[j]+.5, ts[j]+.5),
            y=c(-.1,1.1,1.1,-.1),col='gray95',border='gray95')
  }} else {NA}
  
  
  polygon(x=c(2005,2021,2021,2005),
          y=c(tmp.pop$ci.lo,tmp.pop$ci.lo,tmp.pop$ci.hi,tmp.pop$ci.hi),
          col='gray95',border='gray95')
  
  abline(h=tmp.pop$med,col='gray')
  segments(time.sample,y0=tmp$ci.lo,y1=tmp$ci.hi)
  points(time.sample,
         tmp$med,pch=21,cex=1,
         col="black",bg='white')
  
  
  text(x=2005.5,y=.95,siteNames[position.index[i]],pos=4)
  ifelse(i%in%c(16:20),axis(1L,c(2006,2010,2014,2018),las=3),NA)
  ifelse(i%in%c(1,6,11,16),axis(2, las = 1),NA)
  ifelse(i%in%c(5), legend(x = 15, y = 1,
                           col = c('gray','orange'),
                           lty = c(1,1),
                           legend = c("Persistence only","Persistence & viability"),
                           cex=.55,
                           box.lty=0), NA)
}
mtext("Year", side = 1, outer = TRUE, line = 2.2)
mtext("Probability of seedling survival to fruiting", side = 2, outer = TRUE, line = 2.2)
#mtext("Population*year-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)

sigma0=MCMCchains(mcmcSamples,params="sigma0")

x.sum=apply(sigma0,2,quantile,probs=c(0.025,.25,.5,.75,.975))
index=order(x.sum[3,],decreasing=TRUE)
x.sum=x.sum[,index]
par(mfrow=c(1,1))
plot(NA,NA,type='n',xlim=c(0,3),ylim=c(0,20),
     axes=FALSE,frame=FALSE,
     xlab="",ylab="")
y.pt = 1:20
for(i in 1:20){
  tmp<-x.sum[,i]
  segments(x0=tmp[1],x1=tmp[5],y0=y.pt[i])
  segments(x0=tmp[2],x1=tmp[4],y0=y.pt[i],lwd=3)
  points(x=tmp[3],y=y.pt[i],pch=21,bg='white')
  
}
axis(1,  seq(0,3,by=.5), col.ticks = 1)
axis(2, (1:20),
     labels = siteNames[index], las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
mtext("sigma0",
      side=1,line=2.5,adj=.5,col='black',cex=1)

df.sum=censusSeedlingsFruitingPlants %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(var(seedlingNumber,na.rm=TRUE))

sigma=MCMCchains(mcmcSamples,params="sigma")
x.sum=apply(sigma,2,quantile,probs=c(0.025,.25,.5,.75,.975))

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

time.sample = 1:15+2005
for(i in 1:20){
  
  index=grep(paste0("\\[",i,","),colnames(x.sum))
  tmp=x.sum[,index]
  
  plot(NA,NA,
       ylim=c(0,4),pch=16,xlim=c(2006,2020),
       ylab='',xlab='',xaxt='n',yaxt='n')
  
  segments(time.sample,y0=tmp[1,],y1=tmp[5,])
  points(time.sample,
         tmp[3,],pch=21,cex=1,
         col="black",bg='white')
  
  
  text(x=2005.5,y=3.75,siteNames[i],pos=4)
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
mtext("sigma", side = 2, outer = TRUE, line = 2.2)
#mtext("Population*year-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)

sigma0=MCMCchains(mcmcSamples,params="sigma0")

