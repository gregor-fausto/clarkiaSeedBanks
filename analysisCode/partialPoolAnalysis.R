# -------------------------------------------------------------------
# Analysis of spatial pattern in seed vital rates
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

f<-function(x="alphaS1"){
  chain<-MCMCchains(zc,params = x)
  p<-boot::inv.logit(chain)
  BCI <- t(apply(p,2,FUN = function(x) quantile(x, c(.025, .5, .975))))
  return(BCI)
}

# -------------------------------------------------------------------
# -------------------------------------------------------------------

load("~/Dropbox/modelsF2019/output/partialPoolPlotsFit")

# Analysis of spatial pattern in vital rates
setwd("~/Dropbox/projects/clarkiaScripts/data/reshapeData")
position<-read.csv(file="siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,northing,elevation) %>%
  dplyr::mutate(easting=easting/1000,northing=northing/1000)

# change par 
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 

# S1
ci <- position %>% 
  dplyr::bind_cols(data.frame(f(x="mu.sig")))
names(ci)[5:7] <- c("lo","med","hi")

plot(ci$easting,ci$med,ylim=c(0,1),
     ylab="Seed survival probability [P(S1)]",
     xlab="Easting (km)", 
     cex.lab = 1.5, cex.axis = 1.5,
     pch=16)
segments(x0=position$easting, y0=ci$lo, y1=ci$hi)

head(dat)
glm.dat=cbind(dat$totalJan,(dat$seedStart-dat$totalJan))
library(lme4)
glm1<-glmer(glm.dat ~ 0 + site + (1|yearStart), data = dat, family = binomial)

ci2<-position %>% 
  dplyr::bind_cols(data.frame(est=boot::inv.logit(fixef(glm1))))
#names(ci)[5:7] <- c("lo","med","hi")

# glmer and bayesian model provide similar estimates for site level means
# difference from here and 2011 paper may be
# s1 does not incorporate viability
# s1 is total rather than intact
plot(ci2$easting,ci2$est,ylim=c(0,1),
     ylab="Seed survival probability in first winter [P(S1)]",
     xlab="Easting (km)", 
     cex.lab = 1.5, cex.axis = 1.5,
     pch=16)
points(ci$easting,ci$med,pch=16,col='red')


head(dat)
glm.dat=cbind(dat$seedling,(dat$totalJan-dat$seedling))
library(lme4)
glm1<-glmer(glm.dat ~ 0 + site + (1|yearStart), data = dat, family = binomial)

ci2<-position %>% 
  dplyr::bind_cols(data.frame(est=boot::inv.logit(fixef(glm1))))

plot(ci2$easting,ci2$est,ylim=c(0,1),
     ylab="Seedling germination probability [P(G)]",
     xlab="Easting (km)", 
     cex.lab = 1.5, cex.axis = 1.5,
     pch=16)
ci2
lm(ci2$est~ci2$easting) %>% summary
