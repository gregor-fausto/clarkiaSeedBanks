#################################################################################
################################################################################
################################################################################
# Code for figures for seeds per fruit
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

directory = "/Users/Gregor/Dropbox/dataLibrary/posteriors/"
simFiles <- paste0(directory,list.files(directory))

mcmcSamples <- readRDS(simFiles[[grep("fitnessSamples.rds",simFiles)]])

# dirPars = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/parameterSummary/"
# dir.create(file.path(dirPars), showWarnings = FALSE)

################################################################################
# Data directory
#################################################################################

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/fitness-weakly-informative/"
dataFiles <- paste0(directory,list.files(directory))

data <- readRDS(dataFiles[[grep("data.rds",dataFiles)]])


################################################################################
# Create composite
#################################################################################

summary.fun = function(x){
  quant=quantile(x,c(.025,.975))
  mu = mean(x)
  med = median(x)
  hpdi.interval = HPDI(x, prob = .95)
  return(c(mu,quant,med,hpdi.interval))
}

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/fitness/"
dataFiles <- paste0(directory,list.files(directory))

countSeedPerUndamagedFruit <- readRDS(dataFiles[[grep("countSeedPerUndamagedFruit",dataFiles)]])

# recover site indices and year indices
siteIndex <- data.frame(siteIndex=unique(countSeedPerUndamagedFruit$site3),site=unique(data$site3))
yearIndex <- data.frame(yearIndex=unique(countSeedPerUndamagedFruit$year3),
                        year=unique(data$year3)) 

siteNames=siteIndex$siteIndex

parm.mu_seeds=MCMCchains(mcmcSamples, params=c("mu_seeds"))
parm.mu_seeds.sum = apply(parm.mu_seeds,2,summary.fun)

df.list = list()
for(i in 1:20){
  obj = parm.mu_seeds.sum
  index=grep(paste0("\\[",i,","),colnames(obj))
  tmp = signif(obj[,index],3)
  tmp.df=data.frame(site=(siteIndex[,1])[i],year=yearIndex[,1],(t(tmp)))
  names(tmp.df) = c("site","year","mu","ci.lo95","ci.hi95","med","hpdi.lo95","hpdi.hi95")
  rownames(tmp.df) = NULL
  df.list[[i]] = tmp.df
}

summary.parm.mu_seeds.df=do.call(rbind,df.list)

################################################################################
# Summarize
#################################################################################

dirPars = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/parameterSummary/"
dir.create(file.path(dirPars), showWarnings = FALSE)

write.csv(summary.parm.mu_seeds.df,file=paste0(dirPars,"seedsSummary.csv"))

################################################################################
# Make summary plots
#################################################################################

summary.fun = function(x){
  med = median(x)
  hpdi.interval = HPDI(x, prob = .89)
  hpdi.interval2 = HPDI(x, prob = .5)
  return(c(med,hpdi.interval,hpdi.interval2))
}

parm.mu0=MCMCchains(mcmcSamples,params='nu_seeds')
parm.prob0=apply(parm.mu0,2,exp)
parm.prob0.sum = apply(parm.prob0,2,summary.fun)

df.list = list()
for(i in 1:20){
  obj = parm.prob0.sum
  index=grep(paste0("\\[",i,"\\]"),colnames(obj))
  tmp = signif(obj[,index],3)
  tmp.df=data.frame(site=siteNames[i],(t(tmp)))
  names(tmp.df) = c("site","med","ci.lo","ci.hi","ci.lo2","ci.hi2")
  rownames(tmp.df) = NULL
  df.list[[i]] = tmp.df
}

summary.pop.df=do.call(rbind,df.list)


parm.mu_seeds.sum = apply(parm.mu_seeds,2,summary.fun)

df.list = list()
for(i in 1:20){
  obj = parm.mu_seeds.sum
  index=grep(paste0("\\[",i,","),colnames(obj))
  tmp = signif(obj[,index],3)
  tmp.df=data.frame(site=(siteIndex[,1])[i],year=yearIndex[,1],(t(tmp)))
  names(tmp.df) = c("site","year","med","ci.lo","ci.hi","ci.lo2","ci.hi2")
  rownames(tmp.df) = NULL
  df.list[[i]] = tmp.df
}

summary.parm.mu_seeds.df=do.call(rbind,df.list)

position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,dominant.surface.rock.type) %>%
  dplyr::mutate(easting=easting/1000)



# rename
sigma_p = summary.pop.df %>% dplyr::select(site,med)
sigma_py = summary.parm.mu_seeds.df 

tmp = sigma_py %>%
  dplyr::rename(lambda.med = med) %>%
  dplyr::left_join(sigma_p,by = "site")


interannualSeedsDF =  tmp %>%
  dplyr::arrange(med) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(desc(lambda.med),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  mutate(id = row_number())

interannualSeeds <- interannualSeedsDF %>%
  ggplot(aes(x = id , y = lambda.med)) + 
  geom_hline(aes(yintercept=med),linetype='dotted') +
  geom_point(aes(color=year)) +
  geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +
  coord_flip() +
  facet_grid(site ~ ., scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ylab("Seeds per fruit (undamaged fruits)") +
  # ylim(c(0,1)) +
  labs(color="Year") 

sigma_p = summary.pop.df

spatialSeeds <- sigma_p %>%
  dplyr::left_join(position,by="site") %>%
  ggplot(aes(x = easting , y = med)) + 
  geom_point() +
  geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  ylab("Seeds per fruit (undamaged fruits)") +
  xlab("Easting (km)")

dirFigures = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/"

ggsave(filename=paste0(dirFigures,"interannualSeeds.pdf"),
       plot=interannualSeeds,width=6,height=12)

ggsave(filename=paste0(dirFigures,"spatialSeeds.pdf"),
       plot=spatialSeeds,width=6,height=6)


################################################################################
# Make summary plots
#################################################################################

summary_py = sigma_py
summary_p = summary.pop.df

pdf("~/Dropbox/clarkiaSeedBanks/products/figures/seedsPerFruit-population.pdf",width=8,height=6)
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1    )

time.sample = 1:13+2005
for(i in 1:20){
  
  tmp=summary_py[summary_py$site==siteNames[i],]
  tmp.pop=summary_p[summary_p$site==siteNames[i],]
  
  tmp.vec=c()
  for(j in 1:length(time.sample)){
    tmp.vec[j]=sum(data$sdno[data$site3==i&data$year3==j],na.rm=TRUE)
  }

  tmp.vec=ifelse(is.na(tmp.vec),0,tmp.vec)
  
  
  
  plot(NA,NA,
       ylim=c(0,max(tmp[,3:7])),pch=16,xlim=c(2006,2019),
       ylab='',xlab='',xaxt='n',yaxt='n')
  
  index=tmp.vec==0
  if (sum(index)>0) {for(j in 1:length(time.sample[index])){
    ts=time.sample[index]
    polygon(x=c(ts[j]-.5,ts[j]-.5,
                ts[j]+.5, ts[j]+.5),
            y=c(-.1,max(tmp[,3:7])*2,max(tmp[,3:7])*2,-.1),col='gray95',border='gray95')
  }} else {NA}
  
  
  polygon(x=c(2005,2020,2020,2005),
          y=c(tmp.pop$ci.lo,tmp.pop$ci.lo,tmp.pop$ci.hi,tmp.pop$ci.hi),
          col='gray95',border='gray95')
  
  abline(h=tmp.pop$med,col='gray')
  segments(time.sample,y0=tmp$ci.lo,y1=tmp$ci.hi)
  
  points(time.sample,
         tmp$med,pch=21,cex=1,
         col="black",bg='white')
  
  
  text(x=2005.5,y=.1*max(tmp[,3:7]),siteNames[i],pos=4)
  ifelse(i%in%c(16:20),axis(1L),NA)
  axis(2, seq(0,max(tmp[,3:7]),by=5), tick=FALSE,
       labels = seq(0,max(tmp[,3:7]),by=5), las = 1, 
       cex.axis = 1,line=0,mgp=c(3,.25,0))
  #ifelse(i%in%c(1,6,11,16),,NA)
  ifelse(i%in%c(5), legend(x = 15, y = 1,
                           col = c('gray','orange'),
                           lty = c(1,1),
                           legend = c("Persistence only","Persistence & viability"),
                           cex=.55,
                           box.lty=0), NA)
}
mtext("Year", side = 1, outer = TRUE, line = 2.2)
mtext("Seeds per fruit", side = 2, outer = TRUE, line = 2.2)
#mtext("Population*year-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)
dev.off()
