# script to follow data visualization workflow from Gabry
# with seed bag data
# rainclouds from Allen et al. 2019

library(tidyverse)
library(reshape2)

climate = read.csv(file="/Users/Gregor/Dropbox/clarkia-LTREB/weather data/datafile_demography_site_environmental_data_2005-2020.csv")

climate <- climate %>%
  dplyr::select(c(SITE,intense.demography.,Twinter05.06:Pspring20)) %>%
  dplyr::rename(site=SITE,demography=intense.demography.)

envMolten<-melt(climate, id.vars = c("site", "demography"))

envMolten$season=NA
envMolten$season=ifelse(grepl("Pspring",envMolten$variable),"P_spring",envMolten$season)
envMolten$season=ifelse(grepl("Pwinter",envMolten$variable),"P_winter",envMolten$season)
envMolten$season=ifelse(grepl("Psummer",envMolten$variable),"P_summer",envMolten$season)
envMolten$season=ifelse(grepl("Tspring",envMolten$variable),"T_spring",envMolten$season)
envMolten$season=ifelse(grepl("Twinter",envMolten$variable),"T_winter",envMolten$season)
envMolten$season=ifelse(grepl("Tsummer",envMolten$variable),"T_summer",envMolten$season)

envMolten$year=as.numeric(paste0(20,stringr::str_sub(envMolten$variable, start= -2)))

envMolten$variable2=tolower(stringr::str_sub(envMolten$variable,start=1,end=1))
envMolten$season=stringr::str_sub(envMolten$season,start=3)

climate=envMolten %>%
  dplyr::select(-variable) %>%
  dplyr::select(site,year,demography,season,variable2,value) %>%
  dplyr::rename(intenseDemography=demography) %>%
  dplyr::rename(variable=variable2)

saveRDS(climate,file="/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/climateData-2021.RDS")

position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>%
  dplyr::select(site,easting,dominant.surface.rock.type,elevation) %>%
  dplyr::mutate(easting=easting/1000)

climate <- climate %>%
  dplyr::left_join(position,by=c("site"))

#climate=readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/climateData.RDS")
climate=climate %>% dplyr::ungroup()
climate <- climate %>% dplyr::filter(intenseDemography==1)

climate = climate %>% tidyr::spread(variable,value)

library(tidyverse)

climate.w = climate %>% dplyr::filter(season=="winter")
climate.w=climate.w[!is.na(climate.w$t ),]
climate.sp = climate %>% dplyr::filter(season=="spring")
climate.sp=climate.sp[!is.na(climate.sp$t ),]


#climate.w$site=droplevels(climate.w$site)
climate.w=climate.w[order(climate.w$easting),]

#climate.sp$site=droplevels(climate.sp$site)
climate.sp=climate.sp[order(climate.sp$easting),]

climate.w$site <- factor(climate.w$site , levels=unique(climate.w$site[order(climate.w$easting)]))
climate.sp$site <- factor(climate.sp$site , levels=unique(climate.sp$site[order(climate.sp$easting)]))

plot(climate.w$easting,climate.w$t)
plot(climate.w$easting,climate.w$p)

f=function(x){
  tmp=density(x,from=min(x),to=max(x))
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



df.sum=climate.w %>% dplyr::group_by(site) %>% 
  dplyr::summarise(mu.t = mean(t))
values = df.sum$mu.t
## Use n equally spaced breaks to assign each value to n-1 equal sized bins 
ii <- cut(values, breaks = seq(min(values), max(values), len = 20), 
          include.lowest = TRUE)
## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
colors <- colorRampPalette(c("purple", "orange"))(19)[ii]

## This call then also produces the plot below
image(seq_along(values), 1, as.matrix(seq_along(values)), col = colors,
      axes = F)

polygon(f(climate.w$t)[,1],f(climate.w$t)[,2],col='gray95',border=0)
siteNames=unique(climate.w$site)
dev.off()

pdf(
  "~/Dropbox/clarkiaSeedBanks/products/figures/analysis/climate-plots.pdf",
  height = 6, width = 8)
par(mfrow=c(1,1))
# Figure showing winter temperature gradient
plot(x=f(climate.w$t)[,1],y=f(climate.w$t)[,2],
     type='n',
     xlim=c(0,21),ylim=c(5,13),
     axes=FALSE,frame=FALSE,
     ylab="Winter temperature (\u00B0C)",
     xlab="")
for(i in 1:20){
  climate.tmp=climate.w[climate.w$site==siteNames[i],]
  polygon(i-f(climate.tmp$t)[,2]*1.5,f(climate.tmp$t)[,1],col=colors[i],border=0)
  f.boxplot(i-.025,climate.tmp$t)
  points(rep(i+.25,length(climate.tmp$t))+rnorm(length(climate.tmp$t),0,.025),climate.tmp$t, pch = 16, cex = .5)
}
axis(2, seq(4,14,by=1),
     labels = seq(4,14,by=1), las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(1, (1:20),
     labels = siteNames, las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)


# Spring temperature 
df.sum=climate.sp %>% dplyr::group_by(site) %>% 
  dplyr::summarise(mu.t = mean(t))
values = df.sum$mu.t
## Use n equally spaced breaks to assign each value to n-1 equal sized bins 
ii <- cut(values, breaks = seq(min(values), max(values), len = 20), 
          include.lowest = TRUE)
## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
colors <- colorRampPalette(c("purple", "orange"))(19)[ii]

## This call then also produces the plot below
# image(seq_along(values), 1, as.matrix(seq_along(values)), col = colors,
#       axes = F)

siteNames=unique(climate.sp$site)

# Figure showing spring temperature gradient
plot(x=f(climate.sp$t)[,1],y=f(climate.sp$t)[,2],
     type='n',
     xlim=c(0,21),ylim=c(12,20),
     axes=FALSE,frame=FALSE,
     ylab="Spring temperature (\u00B0C)",
     xlab="")
for(i in 1:20){
  climate.tmp=climate.sp[climate.sp$site==siteNames[i],]
  polygon(i-f(climate.tmp$t)[,2]*1.5,f(climate.tmp$t)[,1],col=colors[i],border=0)
  f.boxplot(i-.025,climate.tmp$t)
  points(rep(i+.25,length(climate.tmp$t))+rnorm(length(climate.tmp$t),0,.025),climate.tmp$t, pch = 16, cex = .5)
}
axis(2, seq(12,20,by=1),
     labels = seq(12,20,by=1), las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(1, (1:20),
     labels = siteNames, las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)
# 


# 
# 
# # Winter precipitation gradient
df.sum=climate.w %>% dplyr::group_by(site) %>%
  dplyr::summarise(mu.p = mean(p))
values = df.sum$mu.p
## Use n equally spaced breaks to assign each value to n-1 equal sized bins
ii <- cut(values, breaks = seq(min(values), max(values), len = 20),
          include.lowest = TRUE)
## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
colors <- colorRampPalette(c("gold", "#6495ed"))(19)[ii]

## This call then also produces the plot below
# image(seq_along(values), 1, as.matrix(seq_along(values)), col = colors,
#       axes = F)

polygon(f(climate.w$t)[,1],f(climate.w$t)[,2],col='gray95',border=0)
siteNames=unique(climate.w$site)

par(mfrow=c(1,1))
plot(f(climate.w$p)[,1],f(climate.w$p)[,2],
     type='n',
     xlim=c(0,21),ylim=c(0,500),
     axes=FALSE,frame=FALSE,
     ylab="Winter precipitation (mm)",
     xlab="")
for(i in 1:20){
  climate.tmp=climate.w[climate.w$site==siteNames[i],]
  polygon(i-f(climate.tmp$p)[,2]*100,f(climate.tmp$p)[,1],col=colors[i],border=0)
  f.boxplot(i-.025,climate.tmp$p)
  points(rep(i+.25,length(climate.tmp$p)),climate.tmp$p, pch = 16, cex = .5)
}
axis(2, seq(0,500,by=50),
     labels = seq(0,500,by=50), las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(1, (1:20),
     labels = siteNames, las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)


df.sum=climate.sp %>% dplyr::group_by(site) %>%
  dplyr::summarise(mu.p = mean(p))
values = df.sum$mu.p
## Use n equally spaced breaks to assign each value to n-1 equal sized bins
ii <- cut(values, breaks = seq(min(values), max(values), len = 20),
          include.lowest = TRUE)
## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
colors <- colorRampPalette(c("gold", "#6495ed"))(19)[ii]

siteNames=unique(climate.sp$site)

plot(f(climate.sp$p)[,1],f(climate.sp$p)[,2],
     type='n',
     xlim=c(0,21),ylim=c(0,200),
     axes=FALSE,frame=FALSE,
     ylab="Spring precipitation (mm)",
     xlab="")
for(i in 1:20){
  climate.tmp=climate.sp[climate.sp$site==siteNames[i],]
  polygon(i-f(climate.tmp$p)[,2]*35,f(climate.tmp$p)[,1],col=colors[i],border=0)
  f.boxplot(i-.025,climate.tmp$p)
  points(rep(i+.25,length(climate.tmp$p)),climate.tmp$p, pch = 16, cex = .5)
}
axis(2, seq(0,200,by=50),
     labels = seq(0,200,by=50), las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(1, (1:20),
     labels = siteNames, las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)


clim.bound=climate.w %>% dplyr::left_join(climate.sp,by=c('site','year'))
par(mfrow=c(1,1))
plot(clim.bound$t.x,clim.bound$p.x, pch = 16, xlab="Winter temperature (\u00B0C)", ylab="Winter precipitation (mm)")
text(8.75,580,paste0("Pearson's r=",signif(cor(clim.bound$t.x,clim.bound$p.x,use="complete.obs"),2)),pos=4)

plot(clim.bound$t.y,clim.bound$p.y, pch = 16, xlab="Spring temperature (\u00B0C)", ylab="Spring precipitation (mm)")
text(15.5,340,paste0("Pearson's r=",signif(cor(clim.bound$t.y,clim.bound$p.y,use="complete.obs"),2)),pos=4)

plot(clim.bound$t.x,clim.bound$t.y, pch = 16, xlab="Winter temperature (\u00B0C)", ylab="Spring temperature (\u00B0C)")
text(4.75,19,paste0("Pearson's r=",signif(cor(clim.bound$t.x,clim.bound$t.y,use="complete.obs"),2)),pos=4)

plot(clim.bound$p.x,clim.bound$p.y, pch = 16, xlab="Winter precipitation (mm)", ylab="Spring precipitation (mm)")
text(300,340,paste0("Pearson's r=",signif(cor(clim.bound$p.x,clim.bound$p.y,use="complete.obs"),2)),pos=4)

plot(clim.bound$t.x,clim.bound$p.y, pch = 16, xlab="Winter temperature (\u00B0C)", ylab="Spring precipitation (mm)")
text(4.75,320,paste0("Pearson's r=",signif(cor(clim.bound$t.x,clim.bound$p.y,use="complete.obs"),2)),pos=4)


dev.off()


ggplot(climate.sp) +
  geom_point(aes(x=easting,y=p)) +
  facet_wrap(~year,scales='free_y') +
  theme_bw()

ggplot(climate.sp %>% 
         dplyr::filter(year<2011)) +
  geom_point(aes(x=easting,y=p)) +
  facet_wrap(~year,scales='free_y') +
  theme_bw()

df=climate.sp %>% 
  dplyr::group_by(site) %>% 
  dplyr::summarise(mu = mean(p),
                   sd= sd(p))



ggplot(climate.sp %>% 
         dplyr::left_join(df,by='site')) +
  geom_point(aes(x=easting,y=mu,group=year)) +
  geom_label(aes(x=easting,y=mu,label=site)) +
  theme_bw()

f.cv = function(x){ sd(x)/mean(x) }
# sliding window calculation
library("slider")

roll.cv=climate.sp %>%
  group_by(site) %>%
  mutate(cv_roll = slide_dbl(p, f.cv, .before = Inf, .complete = TRUE)) %>%
  mutate(mu_roll = slide_dbl(p, mean, .before = Inf, .complete = TRUE)) %>%
  mutate(sd_roll = slide_dbl(p, sd, .before = Inf, .complete = TRUE))

ggplot(climate.w%>% dplyr::filter(site%in%c("LO","GCN")))+
         geom_line(aes(x=year,y=p,group=site,color=site)) +
  theme_bw() +xlim(c(2005,2020))

ggplot(data=roll.cv %>% dplyr::filter(year>2009) %>%
         dplyr::filter(site%in%c("LO","CF","CP3"))) +
  geom_line(aes(x=year,y=mu_roll,group=site,color=easting)) +
  geom_label(aes(x=year,y=mu_roll,label=site)) +
  theme_bw()


ggplot(data=roll.cv) +
  geom_line(aes(x=year,y=cv_roll,group=site)) +
  geom_label(data=roll.cv %>% 
               dplyr::filter(year%in%c(2020)),
             aes(x=year,y=cv_roll,label=site)) +
  theme_bw()
  
roll.cv %>% dplyr::filter(year %in% c(2010,2020))

ggplot(data=roll.cv %>% dplyr::filter(year>2009) ) +
  geom_line(aes(x=year,y=mu_roll,group=site,color=easting)) +
  theme_bw()

ggplot(data=roll.cv %>% dplyr::filter(year>2009) ) +
  geom_line(aes(x=year,y=sd_roll,group=site,color=easting)) +
  theme_bw()

ggplot(data=roll.cv %>% dplyr::filter(year>2009) ) +
  geom_line(aes(x=year,y=cv_roll,group=site,color=easting)) +
  theme_bw() +
  geom_label(data=roll.cv %>% 
               dplyr::filter(year%in%c(2020)),
             aes(x=year,y=cv_roll,label=site)) 

ggplot(data=roll.cv %>% dplyr::filter(year>2009)) +
  geom_point(aes(x=easting,y=cv_roll)) +
  facet_wrap(~year) + theme_bw()


ggplot(data=roll.cv %>% dplyr::filter(year==2020)) +
  geom_point(aes(x=easting,y=cv_roll)) +
  facet_wrap(~year,scales="free_y") + theme_bw() +
  geom_label(aes(x=easting,y=cv_roll,label=site)) 

# winter

roll.cv=climate.w %>%
  group_by(site) %>%
  mutate(cv_roll = slide_dbl(p, f.cv, .before = Inf, .complete = TRUE)) %>%
  mutate(mu_roll = slide_dbl(p, mean, .before = Inf, .complete = TRUE)) %>%
  mutate(sd_roll = slide_dbl(p, sd, .before = Inf, .complete = TRUE))


ggplot(data=roll.cv %>% dplyr::filter(site=="BR")) +
  geom_line(aes(x=year,y=cv_roll,group=site)) +
  geom_label(data=roll.cv %>% 
               dplyr::filter(year%in%c(2020)),
             aes(x=year,y=cv_roll,label=site)) +
  theme_bw()


ggplot(data=roll.cv %>% dplyr::filter(year>2009) %>%
         dplyr::filter(site%in%c("LO","CF","S22"))) +
  geom_line(aes(x=year,y=mu_roll,group=site,color=easting)) +
  geom_label(aes(x=year,y=mu_roll,label=site)) +
  theme_bw()

ggplot(data=roll.cv %>% dplyr::filter(year>2009) ) +
  geom_line(aes(x=year,y=sd_roll,group=site,color=easting)) +
  theme_bw()

ggplot(data=roll.cv %>% dplyr::filter(year>2009) ) +
  geom_line(aes(x=year,y=cv_roll,group=site,color=easting)) +
  theme_bw() +
  geom_label(data=roll.cv %>% 
               dplyr::filter(year%in%c(2020)),
             aes(x=year,y=cv_roll,label=site)) 

ggplot(data=roll.cv %>% dplyr::filter(year>2009)) +
  geom_point(aes(x=easting,y=cv_roll)) +
  facet_wrap(~year,scales="free_y") + theme_bw() +
  geom_label(data=roll.cv%>% dplyr::filter(year>2009) ,
             aes(x=easting,y=cv_roll,label=site)) 

ggplot(data=roll.cv %>% dplyr::filter(year==2020)) +
  geom_point(aes(x=easting,y=cv_roll)) +
  facet_wrap(~year,scales="free_y") + theme_bw() +
  geom_label(aes(x=easting,y=cv_roll,label=site)) 



pdf(
  "~/Dropbox/clarkiaSeedBanks/products/figures/analysis/precip.pdf",
  height = 6, width = 8)
par(mfrow=c(2,1))
par(mar=c(4,4,2,1))
par(fig=c(0,10,5.5,10)/10)

cv.df=climate.sp %>%
  dplyr::group_by(site,easting) %>%
  dplyr::summarise(sd.p=sd(p,na.rm=TRUE),mu.p=mean(p,na.rm=TRUE)) %>%
  dplyr::mutate(cv=sd.p/mu.p)

plot(cv.df$easting/1000,cv.df$cv,
     ylim=c(.2,.7),
     #axes=FALSE,frame=FALSE,
     ylab="CV spring precipitation",
     xlab="Easting",pch = 21, bg = 'white', lwd=1.5)

# # Figure showing spring precipitation gradient

df.sum=climate.sp %>% dplyr::group_by(site) %>%
  dplyr::summarise(mu.p = mean(p))
values = df.sum$mu.p
## Use n equally spaced breaks to assign each value to n-1 equal sized bins
ii <- cut(values, breaks = seq(min(values), max(values), len = 20),
          include.lowest = TRUE)
## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
colors <- colorRampPalette(c("gold", "#6495ed"))(19)[ii]

siteNames=unique(climate.sp$site)

par(fig=c(0,10,0,6)/10)

par(new=T)
plot(f(climate.sp$p)[,1],f(climate.sp$p)[,2],
     type='n',
     xlim=c(0,21),ylim=c(0,200),
     axes=FALSE,frame=FALSE,
     ylab="Spring precipitation (mm)",
     xlab="")
for(i in 1:20){
  climate.tmp=climate.sp[climate.sp$site==siteNames[i],]
  polygon(i-f(climate.tmp$p)[,2]*35,f(climate.tmp$p)[,1],col=colors[i],border=0)
  f.boxplot(i-.025,climate.tmp$p)
  points(rep(i+.25,length(climate.tmp$p)),climate.tmp$p, pch = 16, cex = .5)
}
axis(2, seq(0,200,by=50),
     labels = seq(0,200,by=50), las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(1, (1:20),
     labels = siteNames, las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)
dev.off()


f.climate=function(x,y,width=.2){
  d = quantile(y,c(.025,.25,.5,.75,.975))
  segments(x0=x,y0=d[1],y1=d[5],lty='solid')
  segments(x0=x,y0=d[2],y1=d[4],lty='solid',lwd=2)
  # rect(xleft=x-width/2,xright=x+width/2,
  #      ybottom=d$stats[2],ytop=d$stats[4],col='white')
  points(x=x,y=d[3],pch=21,bg='white',cex=1)
}



climate=readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/climateData-2021.RDS")
climate=climate %>% dplyr::ungroup()
climate <- climate %>% dplyr::filter(intenseDemography==1)
siteNames = unique(climate$site)

climate = climate %>% tidyr::spread(variable,value)


library(tidyverse)

climate.w = climate %>% dplyr::filter(season=="winter")
climate.w=climate.w[!is.na(climate.w$t ),]
climate.sp = climate %>% dplyr::filter(season=="spring")
climate.sp=climate.sp[!is.na(climate.sp$t ),]

# climate.w$site <- factor(climate.w$site , levels=unique(siteNames))
# climate.sp$site <- factor(climate.sp$site , levels=unique(siteNames))

pdf(
  "~/Dropbox/clarkiaSeedBanks/products/figures/analysis/climate-gradients.pdf",
  height = 6, width = 8)
par(mfrow=c(2,2),mar=c(2,4,0,2),oma=c(2,2,2,2))
plot(x=NA,y=NA,
     type='n',
     xlim=c(340,375),ylim=c(5,13),
     #  axes=FALSE,frame=FALSE,
     ylab="Winter temperature (\u00B0C)",
     xlab="Easting (km)",xaxt='n')
for(i in 1:20){
  climate.tmp=climate.w[climate.w$site==siteNames[i],]
  f.climate(position$easting[i],climate.tmp$t)
}

plot(x=NA,y=NA,
     type='n',
     xlim=c(340,375),ylim=c(12,20),
     #  axes=FALSE,frame=FALSE,
     ylab="Spring temperature (\u00B0C)",
     xlab="Easting (km)",xaxt='n')
for(i in 1:20){
  climate.tmp=climate.sp[climate.sp$site==siteNames[i],]
  f.climate(position$easting[i],climate.tmp$t)
}

plot(x=NA,y=NA,
     type='n',
     xlim=c(340,375),ylim=c(0,600),
     #  axes=FALSE,frame=FALSE,
     ylab="Winter precipitation (mm)",
     xlab="Easting (km)")
for(i in 1:20){
  climate.tmp=climate.w[climate.w$site==siteNames[i],]
  f.climate(position$easting[i],climate.tmp$p)
}
mtext("Easting (km)",line=2,side=1,cex=.75)

plot(x=NA,y=NA,
     type='n',
     xlim=c(340,375),ylim=c(0,275),
     #  axes=FALSE,frame=FALSE,
     ylab="Spring precipitation (mm)",
     xlab="Easting (km)")
for(i in 1:20){
  climate.tmp=climate.sp[climate.sp$site==siteNames[i],]
  f.climate(position$easting[i],climate.tmp$p)
}
mtext("Easting (km)",line=2,side=1,cex=.75)

par(mfrow=c(2,2),mar=c(2,4,0,2),oma=c(2,2,2,2))
clim.sum = climate %>% dplyr::group_by(site,season) %>%
  dplyr::summarise(mu.t = mean(t,na.rm=TRUE),
                   mu.p = mean(p,na.rm=TRUE)) %>%
  dplyr::left_join(position,by='site')

plot(clim.sum[clim.sum$season=="winter",]$easting,clim.sum[clim.sum$season=="winter",]$mu.t,xaxt='n',
     ylab="Winter temperature (\u00B0C)")
plot(clim.sum[clim.sum$season=="spring",]$easting,clim.sum[clim.sum$season=="spring",]$mu.t,xaxt='n',
     ylab="Spring temperature (\u00B0C)")

plot(clim.sum[clim.sum$season=="winter",]$easting,clim.sum[clim.sum$season=="winter",]$mu.p,
     ylab="Winter precipitation (mm)",xlab="Easting (km)")
mtext("Easting (km)",line=2,side=1,cex=.75)
#text(clim.sum[clim.sum$season=="winter",]$easting,clim.sum[clim.sum$season=="winter",]$mu.p,clim.sum[clim.sum$season=="winter",]$site)


plot(clim.sum[clim.sum$season=="spring",]$easting,clim.sum[clim.sum$season=="spring",]$mu.p,
     ylab="Spring precipitation (mm)",xlab="Easting (km)")
mtext("Easting (km)",line=2,side=1,cex=.75)
#text(clim.sum[clim.sum$season=="spring",]$easting,clim.sum[clim.sum$season=="spring",]$mu.p,clim.sum[clim.sum$season=="spring",]$site)




par(mfrow=c(2,2),mar=c(2,4,0,2),oma=c(2,2,2,2))
plot(x=NA,y=NA,
     type='n',
     xlim=c(400,1200),ylim=c(5,13),
     #  axes=FALSE,frame=FALSE,
     ylab="Winter temperature (\u00B0C)",
     xlab="Elevation (m)",xaxt='n')
for(i in 1:20){
  climate.tmp=climate.w[climate.w$site==siteNames[i],]
  f.climate(position$elevation[i],climate.tmp$t)
}

plot(x=NA,y=NA,
     type='n',
     xlim=c(400,1200),ylim=c(12,20),
     #  axes=FALSE,frame=FALSE,
     ylab="Spring temperature (\u00B0C)",
     xlab="Elevation (m)",xaxt='n')
for(i in 1:20){
  climate.tmp=climate.sp[climate.sp$site==siteNames[i],]
  f.climate(position$elevation[i],climate.tmp$t)
}

plot(x=NA,y=NA,
     type='n',
     xlim=c(400,1200),ylim=c(0,600),
     #  axes=FALSE,frame=FALSE,
     ylab="Winter precipitation (mm)",
     xlab="Elevation (m)")
for(i in 1:20){
  climate.tmp=climate.w[climate.w$site==siteNames[i],]
  f.climate(position$elevation[i],climate.tmp$p)
}
mtext("Elevation (m)",line=2,side=1,cex=.75)

plot(x=NA,y=NA,
     type='n',
     xlim=c(400,1200),ylim=c(0,275),
     #  axes=FALSE,frame=FALSE,
     ylab="Spring precipitation (mm)",
     xlab="Elevation (m)")
for(i in 1:20){
  climate.tmp=climate.sp[climate.sp$site==siteNames[i],]
  f.climate(position$elevation[i],climate.tmp$p)
}
mtext("Elevation (m)",line=2,side=1,cex=.75)


par(mfrow=c(2,2),mar=c(2,4,0,2),oma=c(2,2,2,2))
clim.sum = climate %>% dplyr::group_by(site,season) %>%
  dplyr::summarise(mu.t = mean(t,na.rm=TRUE),
                   mu.p = mean(p,na.rm=TRUE)) %>%
  dplyr::left_join(position,by='site')

plot(clim.sum[clim.sum$season=="winter",]$elevation,clim.sum[clim.sum$season=="winter",]$mu.t,xaxt='n',
     ylab="Winter temperature (\u00B0C)")
plot(clim.sum[clim.sum$season=="spring",]$elevation,clim.sum[clim.sum$season=="spring",]$mu.t,xaxt='n',
     ylab="Spring temperature (\u00B0C)")

plot(clim.sum[clim.sum$season=="winter",]$elevation,clim.sum[clim.sum$season=="winter",]$mu.p,
     ylab="Winter precipitation (mm)")
mtext("Elevation (m)",line=2,side=1,cex=.75)
#text(clim.sum[clim.sum$season=="winter",]$elevation,clim.sum[clim.sum$season=="winter",]$mu.p,clim.sum[clim.sum$season=="winter",]$site)

plot(clim.sum[clim.sum$season=="spring",]$elevation,clim.sum[clim.sum$season=="spring",]$mu.p,
     ylab="Spring precipitation (mm)")
mtext("Elevation (m)",line=2,side=1,cex=.75)
#text(clim.sum[clim.sum$season=="spring",]$elevation,clim.sum[clim.sum$season=="spring",]$mu.p,clim.sum[clim.sum$season=="spring",]$site)


par(mfrow=c(1,1),mar=c(4,4,3,3))
plot(clim.sum$easting,clim.sum$elevation,pch=21,lwd=1,
     xlab="Easting (km)",ylab="Elevation (m)")
text(342,1100,paste0("Pearson's r=",signif(cor(clim.sum$easting,clim.sum$elevation),2)),pos=4)
dev.off()

plot(clim.sum$easting,clim.sum$elevation,pch=21,lwd=1,
     xlab="Easting (km)",ylab="Elevation (m)",type='n')
text(342,1100,paste0("Pearson's r=",signif(cor(clim.sum$easting,clim.sum$elevation),2)),pos=4)
text(clim.sum$easting,clim.sum$elevation,clim.sum$site)

# climate=climate %>% dplyr::filter(year %in% 2006:2018)
# climate.w = climate %>% dplyr::filter(season=="winter")
# plot(climate.w$t,climate.w$p)
# climate.sp = climate %>% dplyr::filter(season=="spring")
# plot(climate.sp$t,climate.sp$p)
# climate.su = climate %>% dplyr::filter(season=="summer")
# plot(climate.su$t,climate.su$p)
# 
# par(mfrow=c(3,3),mar=c(2,2,2,2))
# climate
# plot(climate.w$t,climate.sp$t)
# plot(climate.w$t,climate.w$p,type='n')
# plot(climate.w$t,climate.w$p,type='n')
# plot(climate.sp$t,climate.w$p)
# plot(climate.sp$t,climate.sp$p)
# plot(climate.sp$t,climate.w$p,type='n')
# plot(climate.su$t,climate.su$p)
# plot(climate.su$t,climate.sp$p)
# plot(climate.su$t,climate.w$p)
