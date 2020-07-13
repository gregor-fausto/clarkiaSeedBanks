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

# -------------------------------------------------------------------
# Functions for use when analyzing data
# -------------------------------------------------------------------

f<-function(x="alphaS1",ci=.95){
  chain<-MCMCchains(zc,params = x)
  qs = c(.5-ci/2,.5,.5+ci/2)
  BCI <- t(apply(chain,2,FUN = function(x) quantile(x, qs)))
  return(BCI)
}

tidyCI <- function(par="g1",ci=.95){
  ci <- position %>%
    dplyr::bind_cols(data.frame(f(x=par,ci=ci)))
  names(ci)[5:7] <- c("lo","med","hi")
  return(ci)
}

# -------------------------------------------------------------------
# -------------------------------------------------------------------

#load("~/Dropbox/modelsF2019/output/seedbagfit")
zc <- readRDS("~/Dropbox/dataLibrary/posteriors/belowgroundSamplesAllYears.RDS")

# Analysis of spatial pattern in vital rates
#setwd("")
position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>%
  dplyr::select(site,easting,northing,elevation) %>%
  dplyr::mutate(easting=easting/1000,northing=northing/1000)


pdf("~/Dropbox/clarkiaSeedBanks/products/figures/vr_spatial.pdf",
 onefile=TRUE,
 paper="USr",
 height = 7.5, width = 10)

# change par
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0))

# G1
ci <- tidyCI("g1",.95)

plot(ci$easting,ci$med,
     ylim=c(0,1),
     ylab="Germination probability [P(G)]",
     xlab="Easting (km)",
     cex.lab = 1.5, cex.axis = 1.5,
     pch=16)
segments(x0=position$easting, y0=ci$lo, y1=ci$hi)

ci <- tidyCI("g1",.5)
segments(x0=position$easting, y0=ci$lo, y1=ci$hi,lwd=2)

# s1
ci <- tidyCI("s1",.95)

plot(ci$easting,ci$med,
     ylim=c(0,1),
     ylab="Seed survival probability in first winter [P(S1)]",
     xlab="Easting (km)",
     cex.lab = 1.5, cex.axis = 1.5,
     pch=16)
segments(x0=position$easting, y0=ci$lo, y1=ci$hi)

ci <- tidyCI("s1",.5)

segments(x0=position$easting, y0=ci$lo, y1=ci$hi,lwd=2)

# s2

ci <- tidyCI("s2",.95)

plot(ci$easting,ci$med,
     ylim=c(0,1),
     ylab="Seed survival probability in summer [P(S2)]",
     xlab="Easting (km)",
     cex.lab = 1.5, cex.axis = 1.5,
     pch=16)
segments(x0=position$easting, y0=ci$lo, y1=ci$hi)

ci <- tidyCI("s2",.5)

segments(x0=position$easting, y0=ci$lo, y1=ci$hi,lwd=2)

# s3
ci <- tidyCI("s3",.95)

plot(ci$easting,ci$med,
     ylim=c(0,1),
     ylab="Seed survival probability in second winter [P(S3)]",
     xlab="Easting (km)",
     cex.lab = 1.5, cex.axis = 1.5,
     pch=16)
segments(x0=position$easting, y0=ci$lo, y1=ci$hi)

ci <- tidyCI("s3",.5)

segments(x0=position$easting, y0=ci$lo, y1=ci$hi,lwd=2)

dev.off()
