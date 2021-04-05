# script to follow data visualization workflow from Gabry
# with seed bag data
# rainclouds from Allen et al. 2019

#source("/Users/Gregor/Dropbox/clarkiaSeedBanks/library/geomFlatViolin.R")
climate=readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/climateData.RDS")
climate=climate %>% dplyr::ungroup()
climate <- climate %>% dplyr::filter(intenseDemography==1)
#climate <- climate %>% dplyr::filter(year %in% 2005:2007&season=="winter") %>% dplyr::mutate(year = year+1)

library(tidyverse)

climate.w = climate %>% dplyr::filter(season=="winter")
climate.w=climate.w[!is.na(climate.w$t ),]
climate.sp = climate %>% dplyr::filter(season=="spring")
climate.sp=climate.sp[!is.na(climate.sp$t ),]


climate.w$site=droplevels(climate.w$site)
climate.w=climate.w[order(climate.w$easting),]

climate.sp$site=droplevels(climate.sp$site)
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
  segments(x0=x,y0=d$stats[2],y1=d$stats[4],lty='solid')
  # rect(xleft=x-width/2,xright=x+width/2,
  #      ybottom=d$stats[2],ytop=d$stats[4],col='white')
  points(x=x,y=d$stats[3],pch=21,bg='white',cex=.5)

  
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

par(mfrow=c(2,1))
# Figure showing winter temperature gradient
plot(x=f(climate.w$t)[,1],y=f(climate.w$t)[,2],
     type='n',
     xlim=c(0,21),ylim=c(5,13),
     axes=FALSE,frame=FALSE,
     ylab="Winter temperature",
     xlab="")
for(i in 1:20){
  climate.tmp=climate.w[climate.w$site==siteNames[i],]
  polygon(i-f(climate.tmp$t)[,2]*1.5,f(climate.tmp$t)[,1],col=colors[i],border=0)
  points(rep(i+.25,length(climate.tmp$t))+rnorm(length(climate.tmp$t),0,.025),climate.tmp$t, pch = 16, cex = .5)
}
axis(2, seq(4,14,by=1),
     labels = seq(4,14,by=1), las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(1, (1:20),
     labels = siteNames, las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)


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
     ylab="Spring temperature",
     xlab="")
for(i in 1:20){
  climate.tmp=climate.sp[climate.sp$site==siteNames[i],]
  polygon(i-f(climate.tmp$t)[,2]*1.5,f(climate.tmp$t)[,1],col=colors[i],border=0)
  points(rep(i+.25,length(climate.tmp$t))+rnorm(length(climate.tmp$t),0,.025),climate.tmp$t, pch = 16, cex = .5)
}
axis(2, seq(12,20,by=1),
     labels = seq(12,20,by=1), las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(1, (1:20),
     labels = siteNames, las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)

clim.bound=climate.w %>% dplyr::left_join(climate.sp,by=c('site','year'))
plot(clim.bound$t.x,clim.bound$t.y)
cor(clim.bound$t.x,clim.bound$t.y,use="complete.obs")

# 
# 
# # Figure showing winter precipitation gradient

df.sum=climate.w %>% dplyr::group_by(site) %>%
  dplyr::summarise(mu.p = mean(p))
values = df.sum$mu.p
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

par(mfrow=c(2,1))
plot(f(climate.w$p)[,1],f(climate.w$p)[,2],
     type='n',
     xlim=c(0,21),ylim=c(0,500),
     axes=FALSE,frame=FALSE,
     ylab="Winter precipitation",
     xlab="")
for(i in 1:20){
  climate.tmp=climate.w[climate.w$site==siteNames[i],]
  polygon(i-f(climate.tmp$p)[,2]*100,f(climate.tmp$p)[,1],col=colors[i],border=0)
  points(rep(i+.25,length(climate.tmp$p)),climate.tmp$p, pch = 16, cex = .5)
}
axis(2, seq(0,500,by=50),
     labels = seq(0,500,by=50), las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(1, (1:20),
     labels = siteNames, las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)


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
     ylim=c(.2,.6),
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
colors <- colorRampPalette(c("orange", "purple"))(19)[ii]

siteNames=unique(climate.sp$site)

par(fig=c(0,10,0,6)/10)

par(new=T)
plot(f(climate.sp$p)[,1],f(climate.sp$p)[,2],
     type='n',
     xlim=c(0,21),ylim=c(0,200),
     axes=FALSE,frame=FALSE,
     ylab="Spring precipitation",
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
