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
  dres <- density( x ,from = x.min, to = x.max)
  modeParam <- dres$x[which.max(dres$y)]
  return(modeParam)
}
# -------------------------------------------------------------------

# read in samples from posterior distributions
rs <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/rsPosterior.RDS")
climate <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/climateData-2021.RDS")
siteNames=climate %>% dplyr::filter(intenseDemography==1) %>% 
  dplyr::select(site) %>% unique
siteNames = siteNames$site
rs.sum=apply(rs,2,posterior.mode)

df.list=list()
for(i in 1:20){
  index=grep(paste0("\\[",i,","),names(rs.sum))
  tmp=rs.sum[index]
  tmp.clim=climate[climate$site==siteNames[i],] %>%
    dplyr::filter(season=="spring"&variable=='p')
 df=data.frame(tmp.clim,rs=tmp)
 df.list[[i]] = df
}

df=do.call(rbind,df.list)

# par(mfrow=c(4,5),mar=c(1,1,1,1))
# for(i in 1:20){
#   tmp=df.list[[i]]
#   plot(tmp$value,tmp$rs)
# }

b.vec=c()
for(i in 1:20){
    tmp=df.list[[i]]
    tmp$value.log=log(tmp$value)
    tmp$rs.log=log(tmp$rs)
    m1=lm((tmp$rs.log)~(tmp$value.log))
    summary(m1)
    b = coef(m1)[2];
    b.vec[i]=b
}

pdf(
  "~/Dropbox/clarkiaSeedBanks/products/figures/analysis/rs-climate-sensitivity.pdf",
  height = 8, width = 8)
index=order(b.vec)
par(mfrow=c(4,5),mar=c(0,.25,.25,0),
    oma=c(4,4,1,1))
p.val=c()
for(i in 1:20){
  tmp=df.list[[index[i]]]
  tmp$value.log=log(tmp$value)
  tmp$rs.log=log(tmp$rs)
  m1=lm((tmp$rs.log)~(tmp$value.log))
  b = coef(m1)[2];  a = coef(m1)[1]
  r.2 = signif(summary(m1)$r.squared,2)
  plot((tmp$value.log),(tmp$rs.log),xlim=c(3,6),ylim=c(-6,7),
       xaxt='n',yaxt='n',pch=16)
  segments(x0=min(tmp$value.log),x1=max(tmp$value.log),
           y0=a+b*min(tmp$value.log),y1=a+b*max(tmp$value.log),
           lwd=1,lty='dotted')
  ifelse(i%in%c(1,6,11,16),axis(2L,las=1),NA)
  ifelse(i%in%c(16:20),axis(1L),NA)
  
  text(3,-3.5, cex=.9,siteNames[index[i]] ,adj=c(0,0))
  text(3,-4.75, cex=.9,paste0("b=",signif(b,2)) ,adj=c(0,0))
  
  p.val[1]=ifelse((summary(m1))$coefficients[,4][2]<.05/20,1,0)
  p.val=as.numeric(p.val)
  p.val=rep("*",sum(p.val))
  p.val=ifelse(length(p.val)==0,"",p.val)
  text(3,-6, cex=.9, bquote(R^2==.(r.2)~.(p.val)),adj=c(0,0))
  
}
mtext("Log(per-capita reproductive success)", side = 2, outer = TRUE, line = 2)
mtext("Log(spring precipitation)", side = 1, outer = TRUE, line = 2.2)

dev.off()
