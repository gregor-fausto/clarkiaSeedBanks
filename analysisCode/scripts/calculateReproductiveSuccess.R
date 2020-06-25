# -------------------------------------------------------------------
# Analysis of fitness models
# Outputs reproductive success estimates
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
# -------------------------------------------------------------------
# Import and organize seed bag data
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Seedling survival to fruiting
# -------------------------------------------------------------------

load("~/Dropbox/modelsF2019/output/sigmaFit")

# summary table
# caterpillar plots
par(mfrow=c(1,1))

MCMCplot(zmP, params = c("alphaS"), horiz = FALSE,
         sz_labels = .4, sz_med = .75, sz_thin=.7, sz_thick = 2,ylim=c(-30,30))
MCMCplot(zmP, params = c("betaS"), horiz = FALSE,
         sz_labels = .4, sz_med = .75, sz_thin=.7, sz_thick = 2,ylim=c(-15,15))
MCMCplot(zmP, params = c("gammaS"), horiz = FALSE,
         sz_labels = .4, sz_med = .75, sz_thin=.7, sz_thick = 2,ylim=c(-30,15))

a <- MCMCchains(zmP,params="alphaS")
b <- MCMCchains(zmP,params="betaS")
g <- MCMCchains(zmP,params="gammaS")

gamma_names <- dimnames(g)[[2]]

library(stringr)
vals<-str_extract_all(gamma_names, "[0-9]+")
siteGamma<-as.numeric(unlist(lapply(vals, `[[`, 1)))
yearGamma<-as.numeric(unlist(lapply(vals, `[[`, 2)))

# site-by-year estimates for survival
nsite <- dim(a)[2]
nyear <- dim(b)[2]

sigmaEstimates<-list()
sigma<-matrix(NA,nrow=dim(a)[1],ncol=dim(b)[2])

for(j in 1:nsite){
  for(k in 1:nyear){
    alphaIndex <- j
    betaIndex <- k
    gammaIndex <- intersect(which(siteGamma %in% j), which(yearGamma %in% k))
    
    sigma[,k] = boot::inv.logit(a[,alphaIndex] + b[,betaIndex] + g[,gammaIndex])
  }
  sigmaEstimates[[j]] <- sigma
}

# -------------------------------------------------------------------
# Fruits per plant
# -------------------------------------------------------------------

load("~/Dropbox/modelsF2019/output/fecFit")

# summary table
# caterpillar plots
par(mfrow=c(1,1))

MCMCplot(zmP, params = c("alphaF"), horiz = FALSE,
         sz_labels = .4, sz_med = .75, sz_thin=.7, sz_thick = 2,ylim=c(-30,30))
MCMCplot(zmP, params = c("betaF"), horiz = FALSE,
         sz_labels = .4, sz_med = .75, sz_thin=.7, sz_thick = 2,ylim=c(-15,15))
MCMCplot(zmP, params = c("gammaF"), horiz = FALSE,
         sz_labels = .4, sz_med = .75, sz_thin=.7, sz_thick = 2,ylim=c(-30,15))

a <- MCMCchains(zmP,params="alphaF")
b <- MCMCchains(zmP,params="betaF")
g <- MCMCchains(zmP,params="gammaF")

gamma_names <- dimnames(g)[[2]]

library(stringr)
vals<-str_extract_all(gamma_names, "[0-9]+")
siteGamma<-as.numeric(unlist(lapply(vals, `[[`, 1)))
yearGamma<-as.numeric(unlist(lapply(vals, `[[`, 2)))

# site-by-year estimates for survival
nsite <- dim(a)[2]
nyear <- dim(b)[2]

fecEstimates<-list()
fec<-matrix(NA,nrow=dim(a)[1],ncol=dim(b)[2])

for(j in 1:nsite){
  for(k in 1:nyear){
    alphaIndex <- j
    betaIndex <- k
    gammaIndex <- intersect(which(siteGamma %in% j), which(yearGamma %in% k))
    
    fec[,k] = exp(a[,alphaIndex] + b[,betaIndex] + g[,gammaIndex])
  }
  fecEstimates[[j]] <- fec
}

# -------------------------------------------------------------------
# Seeds per fruit
# -------------------------------------------------------------------

load("~/Dropbox/modelsF2019/output/phiFit")

# summary table
# caterpillar plots
par(mfrow=c(1,1))

MCMCplot(zmP, params = c("alphaP"), horiz = FALSE,
         sz_labels = .4, sz_med = .75, sz_thin=.7, sz_thick = 2,ylim=c(-30,30))
MCMCplot(zmP, params = c("betaP"), horiz = FALSE,
         sz_labels = .4, sz_med = .75, sz_thin=.7, sz_thick = 2,ylim=c(-15,15))
MCMCplot(zmP, params = c("gammaP"), horiz = FALSE,
         sz_labels = .4, sz_med = .75, sz_thin=.7, sz_thick = 2,ylim=c(-30,15))

a <- MCMCchains(zmP,params="alphaP")
b <- MCMCchains(zmP,params="betaP")
g <- MCMCchains(zmP,params="gammaP")

gamma_names <- dimnames(g)[[2]]

library(stringr)
vals<-str_extract_all(gamma_names, "[0-9]+")
siteGamma<-as.numeric(unlist(lapply(vals, `[[`, 1)))
yearGamma<-as.numeric(unlist(lapply(vals, `[[`, 2)))

# site-by-year estimates for survival
nsite <- dim(a)[2]
nyear <- dim(b)[2]

phiEstimates<-list()
phi<-matrix(NA,nrow=dim(a)[1],ncol=dim(b)[2])

for(j in 1:nsite){
  for(k in 1:nyear){
    alphaIndex <- j
    betaIndex <- k
    gammaIndex <- intersect(which(siteGamma %in% j), which(yearGamma %in% k))
    
    phi[,k] = exp(a[,alphaIndex] + b[,betaIndex] + g[,gammaIndex])
  }
  phiEstimates[[j]] <- phi
}

# -------------------------------------------------------------------
# Reproductive success
# -------------------------------------------------------------------
nsite <- 20
nyear <- 6

rsEstimates<-list()
rs<-matrix(NA,nrow=30000,ncol=nyear)

for(j in 1:nsite){
  s <- sigmaEstimates[[j]]
  f <- fecEstimates[[j]]
  p <- phiEstimates[[j]]
  
  for(k in 1:nyear){
    
    rs[,k]<-s[,k]*f[,k]*p[,k]
  }
  rsEstimates[[j]] <- rs
}

# load file with index site*year oof missing data
load("~/Dropbox/modelsF2019/output/missingness")

for(i in 1:(dim(missing_dat)[1])){
  v<-as.numeric(missing_dat[i,])
  rsEstimates[[v[1]]]<-rsEstimates[[(v[1])]][,-(v[2])]
}

# save estimates of RS with sites/years where there is missing data excluded
save(rsEstimates,file="~/Dropbox/modelsF2019/output/rsVarPosterior")

# 
# 
# temporal_variance <- function(x,fun=var){
#   apply(x,1,fun)
# }
# 
# cols_fun <- function(x,fun=var){
#   apply(x,2,fun)
# }
# 
# #lapply(rsEstimates,cols_fun,fun=median)
# 
# # # temporal_variance 
# # d<-lapply(rsEstimates,temporal_variance,fun=var)
# # d2<-matrix(unlist(d), ncol = 20, byrow = FALSE)
# # 
# # BCI <- apply(d2,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
# # HPDI <- apply(d2,2,FUN = function(x) hdi(x, .95))
# 
# # geometric var
# gsd <- function(x){
#   y <- exp(sd(log(x)))
#   return(y)
# }
# 
# d<-lapply(rsEstimates,temporal_variance,fun=gsd)
# d2<-matrix(unlist(d), ncol = 20, byrow = FALSE)
# 
# BCI <- apply(d2,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
# HPDI <- apply(d2,2,FUN = function(x) hdi(x, .95))
# 
# plot(t(BCI)[,2],ylim=c(0,25))
# 
# df<-data.frame(t(BCI))
# names(df) <- c('lo','med','hi')
# df<-cbind(df,siteNo=seq(1:20))
# df %>% dplyr::filter(med>20)
# 
# #probs <- c(4,9,11,12,13,16,18)
# 
# probs <- c(9,11,18)
# 
# lapply(sigmaEstimates[probs],cols_fun)
# lapply(fecEstimates[probs],cols_fun)
# lapply(phiEstimates[probs],cols_fun)


