fileDirectory = "/Volumes/RUGGEDKEY/mcmcSamplesThinned/"
list.files(fileDirectory)

mcmcSamples = readRDS(paste0(fileDirectory,"seedSamples.rds"))
library(MCMCvis)
mu0_g=MCMCchains(mcmcSamples,params="mu0_g")
g=apply(mu0_g,2,boot::inv.logit)
g
str(g)
g.sum = apply(g,2,quantile,c(.025,.5,.975))
g.sum
par(mfrow=c(4,5),mar=c(2,2,2,2))
for(i in 1:20){
tmp=g.sum[,grep(paste0("\\[",i,","),colnames(g.sum))]
plot(tmp[2,],type='b',ylim=c(0,1))
}

climate = readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/climateData.RDS")
climate <- climate %>% dplyr::filter(intenseDemography==1)
siteNames = unique(climate$site)

t.mu=climate %>% dplyr::group_by(site) %>% dplyr::filter(season=="winter") %>% dplyr::summarise(t.mean = mean(t,na.rm=TRUE))

dev.off()

par(mfrow=c(1,3))
g.1=g.sum[,grep(paste0(",",1,"\\]"),colnames(g.sum))]
plot(t.mu$t.mean,g.1[2,])

g.2=g.sum[,grep(paste0(",",2,"\\]"),colnames(g.sum))]
plot(t.mu$t.mean,g.2[2,])

g.3=g.sum[,grep(paste0(",",3,"\\]"),colnames(g.sum))]
plot(t.mu$t.mean,g.3[2,])




p.mu=climate %>% dplyr::group_by(site) %>% dplyr::filter(season=="winter") %>% dplyr::summarise(p.mean = mean(p,na.rm=TRUE))

dev.off()

par(mfrow=c(1,3))
g.1=g.sum[,grep(paste0(",",1,"\\]"),colnames(g.sum))]
plot(p.mu$p.mean,g.1[2,])

g.2=g.sum[,grep(paste0(",",2,"\\]"),colnames(g.sum))]
plot(p.mu$p.mean,g.2[2,])

g.3=g.sum[,grep(paste0(",",3,"\\]"),colnames(g.sum))]
plot(p.mu$p.mean,g.3[2,])



library(MCMCvis)
mu_g=MCMCchains(mcmcSamples,params="mu_g")
g=apply(mu_g,2,boot::inv.logit)
g1.1=grep(",1\\]",colnames(g))
g1.2=grep(",4\\]",colnames(g))
g1.3=grep(",6\\]",colnames(g))

par(mfcol=c(2,3))
pop.g1.1=apply(g[,g1.1],2,quantile,c(.025,.5,.975))
year1=climate %>% dplyr::filter(year==2005) %>% dplyr::filter(season=="winter")
plot(year1$t,pop.g1.1[2,]);mod=lm(pop.g1.1[2,]~year1$t);abline(a=coef(mod)[1],b=coef(mod)[2])
plot(year1$p,pop.g1.1[2,]);mod=lm(pop.g1.1[2,]~year1$p);abline(a=coef(mod)[1],b=coef(mod)[2])

pop.g1.2=apply(g[,g1.2],2,quantile,c(.025,.5,.975))
year2=climate %>% dplyr::filter(year==2006) %>% dplyr::filter(season=="winter")
plot(year2$t,pop.g1.2[2,]);mod=lm(pop.g1.2[2,]~year2$t);abline(a=coef(mod)[1],b=coef(mod)[2])
plot(year2$p,pop.g1.2[2,]);mod=lm(pop.g1.2[2,]~year2$p);abline(a=coef(mod)[1],b=coef(mod)[2])

pop.g1.3=apply(g[,g1.3],2,quantile,c(.025,.5,.975))
year3=climate %>% dplyr::filter(year==2007) %>% dplyr::filter(season=="winter")
plot(year3$t,pop.g1.3[2,]);mod=lm(pop.g1.3[2,]~year3$t);abline(a=coef(mod)[1],b=coef(mod)[2])
plot(year3$p,pop.g1.3[2,]);mod=lm(pop.g1.3[2,]~year3$p);abline(a=coef(mod)[1],b=coef(mod)[2])



par(mfcol=c(1,2))
plot(year1$t,pop.g1.1[2,],xlim=c(5,12),ylim=c(0,.5),pch=16,col='red');mod=lm(pop.g1.1[2,]~year1$t);abline(a=coef(mod)[1],b=coef(mod)[2],col='red')
points(year2$t,pop.g1.2[2,],pch=16,col='orange');mod=lm(pop.g1.2[2,]~year2$t);abline(a=coef(mod)[1],b=coef(mod)[2],col='orange')
points(year3$t,pop.g1.3[2,],pch=16,col='yellow');mod=lm(pop.g1.3[2,]~year3$t);abline(a=coef(mod)[1],b=coef(mod)[2],col='yellow')

plot(year1$p,pop.g1.1[2,],xlim=c(50,200),ylim=c(0,.5),pch=16,col='red');mod=lm(pop.g1.1[2,]~year1$p);abline(a=coef(mod)[1],b=coef(mod)[2],col='red')
points(year2$p,pop.g1.2[2,],pch=16,col='orange');mod=lm(pop.g1.2[2,]~year2$p);abline(a=coef(mod)[1],b=coef(mod)[2],col='orange')
points(year3$p,pop.g1.3[2,],pch=16,col='yellow');mod=lm(pop.g1.3[2,]~year3$p);abline(a=coef(mod)[1],b=coef(mod)[2],col='yellow')


plot(NA,NA,xlim=c(5,12),ylim=c(0,.5))

for(i in 1:20){

  g1.1=grep(paste0("\\[",i,",1\\]"),colnames(g))
  g1.2=grep(paste0("\\[",i,",4\\]"),colnames(g))
  g1.3=grep(paste0("\\[",i,",6\\]"),colnames(g))
  
  tmp=apply(g[,c(g1.1,g1.2,g1.3)],2,quantile,c(.025,.5,.975))
  tmp.clim=climate[climate$site==siteNames[i],] %>% dplyr::filter(year%in% 2005:2007&season=="winter")

  points(tmp.clim$t,tmp[2,],pch=16,col='orange');mod=lm(tmp[2,]~tmp.clim$t);abline(a=coef(mod)[1],b=coef(mod)[2],col='orange')
  
}

plot(NA,NA,xlim=c(50,200),ylim=c(0,.5))

for(i in 1:20){
  
  g1.1=grep(paste0("\\[",i,",1\\]"),colnames(g))
  g1.2=grep(paste0("\\[",i,",4\\]"),colnames(g))
  g1.3=grep(paste0("\\[",i,",6\\]"),colnames(g))
  
  tmp=apply(g[,c(g1.1,g1.2,g1.3)],2,quantile,c(.025,.5,.975))
  tmp.clim=climate[climate$site==siteNames[i],] %>% dplyr::filter(year%in% 2005:2007&season=="winter")
  
  points(tmp.clim$p,tmp[2,],pch=16,col='orange');mod=lm(tmp[2,]~tmp.clim$p);abline(a=coef(mod)[1],b=coef(mod)[2],col='orange')
  
}


par(mfrow=c(4,5),mar=c(1,1,1,1))

for(i in 1:20){
  
  g1.1=grep(paste0("\\[",i,",1\\]"),colnames(g))
  g1.2=grep(paste0("\\[",i,",4\\]"),colnames(g))
  g1.3=grep(paste0("\\[",i,",6\\]"),colnames(g))
  
  tmp=apply(g[,c(g1.1,g1.2,g1.3)],2,quantile,c(.025,.5,.975))
  tmp.clim=climate[climate$site==siteNames[i],] %>% dplyr::filter(year%in% 2005:2007&season=="winter")
  
  plot(tmp.clim$t,tmp[2,],pch=16,col='orange',xlim=c(5,12),ylim=c(0,.5));mod=lm(tmp[2,]~tmp.clim$t);abline(a=coef(mod)[1],b=coef(mod)[2],col='orange')
  
}



for(i in 1:20){
  
  g1.1=grep(paste0("\\[",i,",1\\]"),colnames(g))
  g1.2=grep(paste0("\\[",i,",4\\]"),colnames(g))
  g1.3=grep(paste0("\\[",i,",6\\]"),colnames(g))
  
  tmp=apply(g[,c(g1.1,g1.2,g1.3)],2,quantile,c(.025,.5,.975))
  tmp.clim=climate[climate$site==siteNames[i],] %>% dplyr::filter(year%in% 2005:2007&season=="winter")
  
  plot(tmp.clim$p,tmp[2,],pch=16,col='orange',xlim=c(50,200),ylim=c(0,.5));mod=lm(tmp[2,]~tmp.clim$p);abline(a=coef(mod)[1],b=coef(mod)[2],col='orange')
  
}

par(mfcol=c(1,2))

plot(NA,NA,xlim=c(min(climate$easting),max(climate$easting)),ylim=c(-.1,.1))

for(i in 1:20){
  
  g1.1=grep(paste0("\\[",i,",1\\]"),colnames(g))
  g1.2=grep(paste0("\\[",i,",4\\]"),colnames(g))
  g1.3=grep(paste0("\\[",i,",6\\]"),colnames(g))
  
  tmp=apply(g[,c(g1.1,g1.2,g1.3)],2,quantile,c(.025,.5,.975))
  tmp.easting=climate[climate$site==siteNames[i],]$easting %>% unique
  mod=lm(tmp[2,]~tmp.clim$t);
  points(tmp.easting,coef(mod)[2],pch=16,col='orange');
  text(tmp.easting,coef(mod)[2],siteNames[i],cex=.75)
  
}

plot(NA,NA,xlim=c(min(climate$easting),max(climate$easting)),ylim=c(-.025,.025))

for(i in 1:20){
  
  g1.1=grep(paste0("\\[",i,",1\\]"),colnames(g))
  g1.2=grep(paste0("\\[",i,",4\\]"),colnames(g))
  g1.3=grep(paste0("\\[",i,",6\\]"),colnames(g))
  
  tmp=apply(g[,c(g1.1,g1.2,g1.3)],2,quantile,c(.025,.5,.975))
  tmp.easting=climate[climate$site==siteNames[i],]$easting %>% unique
  mod=lm(tmp[2,]~tmp.clim$p);
  points(tmp.easting,coef(mod)[2],pch=16,col='orange');  
  text(tmp.easting,coef(mod)[2],siteNames[i],cex=.75)
}

