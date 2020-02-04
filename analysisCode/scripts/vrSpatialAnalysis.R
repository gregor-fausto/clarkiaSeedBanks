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

load("~/Dropbox/modelsF2019/output/seedbagfit")

# Analysis of spatial pattern in vital rates
setwd("~/Dropbox/projects/clarkiaScripts/data/reshapeData")
position<-read.csv(file="siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,northing,elevation) %>%
  dplyr::mutate(easting=easting/1000,northing=northing/1000)


pdf(
  "~/Dropbox/modelsF2019/figures/vr_spatial.pdf",
  onefile=TRUE,
  paper="USr",
  height = 7.5, width = 10)

# change par 
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 

# G1
ci <- position %>% 
  dplyr::bind_cols(data.frame(f(x="alphaG1")))
names(ci)[5:7] <- c("lo","med","hi")

plot(ci$easting,ci$med,ylim=c(0.05,.3),
     ylab="Germination probability [P(G)]",
     xlab="Easting (km)", 
     cex.lab = 1.5, cex.axis = 1.5,
     pch=16)
segments(x0=position$easting, y0=ci$lo, y1=ci$hi)

# s1

ci <- position %>% 
  dplyr::bind_cols(data.frame(f(x="alphaS1")))
names(ci)[5:7] <- c("lo","med","hi")

plot(ci$easting,ci$med,ylim=c(0.5,1),
     ylab="Seed survival probability in first winter [P(S1)]",
     xlab="Easting (km)", 
     cex.lab = 1.5, cex.axis = 1.5,
     pch=16)
segments(x0=position$easting, y0=ci$lo, y1=ci$hi)

# s2

ci <- position %>% 
  dplyr::bind_cols(data.frame(f(x="alphaS2")))
names(ci)[5:7] <- c("lo","med","hi")

plot(ci$easting,ci$med,ylim=c(0.5,1),
     ylab="Seed survival probability in summer [P(S2)]",
     xlab="Easting (km)", 
     cex.lab = 1.5, cex.axis = 1.5,
     pch=16)
segments(x0=position$easting, y0=ci$lo, y1=ci$hi)

# s3

ci <- position %>% 
  dplyr::bind_cols(data.frame(f(x="alphaS3")))
names(ci)[5:7] <- c("lo","med","hi")

plot(ci$easting,ci$med,ylim=c(.5,1),
     ylab="Seed survival probability in second winter [P(S3)]",
     xlab="Easting (km)", 
     cex.lab = 1.5, cex.axis = 1.5,
     pch=16)
segments(x0=position$easting, y0=ci$lo, y1=ci$hi)

dev.off()

# 
# lm1<-lm(ci$med~ci$easting)
# abline(a=coef(lm1)[1],b=coef(lm1)[2])
# text(350,.275,paste0("y=",round(coef(lm1)[1],2),"+",round(coef(lm1)[2],3),"x"))
# 
# summary(lm(ci$med~ci$easting))
# 
# ci <- position[-probs,] %>% 
#   dplyr::bind_cols(data.frame(t(rs.BCI)))
# names(ci)[5:7] <- c("lo","med","hi")
# plot(ci$easting,ci$med,
#      ylab="Geometric SD of reproductive success",
#      xlab="Easting (km)",ylim=c(0,16))
# 
# 