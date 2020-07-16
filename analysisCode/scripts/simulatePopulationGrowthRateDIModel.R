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

set.seed(10)
# -------------------------------------------------------------------
# Functions for use when analyzing data
# -------------------------------------------------------------------
temporal_variance <- function(x,fun=var){
  apply(x,1,fun)
}

cols_fun <- function(x,fun=var){
  apply(x,2,fun)
}

# geometric var
gsd <- function(x){
  y <- exp(sd(log(x)))
  return(y)
}


f<-function(x="alphaS1",chain){
  chain<-MCMCchains(chain=belowground,params = x)
  p<-boot::inv.logit(chain)
  BCI <- t(apply(p,2,FUN = function(x) quantile(x, c(.025, .5, .975))))
  return(BCI)
}

# geometric mean from 
# https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
gm_mean = function(x, na.rm=TRUE){
exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# -------------------------------------------------------------------

# read in samples from posterior distributions
belowground <- readRDS("~/Dropbox/dataLibrary/posteriors/belowgroundSamplesAllYears.RDS")
aboveground <- readRDS("~/Dropbox/clarkiaSeedBanks/products/dataFiles2/rsVarFullPosterior.RDS")



# Exclude sites with very high temporal variance (orders of magnitude higher)
# 
# d<-lapply(rsEstimates,temporal_variance,fun=gsd)
# d2<-matrix(unlist(d), ncol = 20, byrow = FALSE)
# 
# BCI <- apply(d2,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
# HPDI <- apply(d2,2,FUN = function(x) hdi(x, .95))
# 
# df<-data.frame(t(BCI))
# names(df) <- c('lo','med','hi')
# df<-cbind(df,siteNo=seq(1:20))
# df.exclude<-df %>% dplyr::filter(med>20)
# 
# # problem sites index
# probs <- df.exclude$siteNo

# G1
g1<-MCMCchains(belowground,params = "g1")
# exclude problem sites
# g1<-g1[,-probs]

g1.med <- apply(g1,2,FUN = function(x) quantile(x, c(.5)))

# S1
s1<-MCMCchains(belowground,params = "s1")
# exclude problem sites
# s1<-s1[,-probs]

s1.med <- apply(s1,2,FUN = function(x) quantile(x, c(.5)))

# S2
s2<-MCMCchains(belowground,params = "s2")
# exclude problem sites
# s2<-s2[,-probs]

s2.med <- apply(s2,2,FUN = function(x) quantile(x, c(.5)))


# S2
s3<-MCMCchains(belowground,params = "s3")
# exclude problem sites
# s3<-s3[,-probs]

s3.med <- apply(s3,2,FUN = function(x) quantile(x, c(.5)))


# RS
# rs<-aboveground %>%
#   dplyr::filter(!is.na(rs)) %>%
#   dplyr::select(site,year,rs)

rs<-lapply(aboveground,cols_fun,median)

#rs<-split(data.frame(rs),rs$site)

fitness <- function(g=g1,s1=s1,s2=s2,s3=s3,rs=rs){
  p1 = g*rs*s1*(s2^(3/8))
  p2 = (1-g)*(s2*s3)
  return(as.numeric(p1+p2))
}

medianEst <- function(x){return(apply(x,2,median))}
g1 = medianEst(g1)
s1 = medianEst(s1)
s2 = medianEst(s2)
s3 = medianEst(s3)

plot(s1,s3,xlim=c(0,1),ylim=c(0,1))
plot(s1,s2,xlim=c(0,1),ylim=c(0,1))
plot(s2,s3,xlim=c(0,1),ylim=c(0,1))
plot(g1,s1,xlim=c(0,1),ylim=c(0,1))
plot(g1,s2,xlim=c(0,1),ylim=c(0,1))
plot(g1,s3,xlim=c(0,1),ylim=c(0,1))
plot(g1,s2*s3,xlim=c(0,1),ylim=c(0,1))
plot(1-g1,s2*s3,xlim=c(0,1),ylim=c(0,1))

logfit<-c()

g = seq(0,1,by=.01)
fit = c()
g.mat = matrix(NA,nrow=1000,ncol=length(g))
g.sites = list()
yt = c()

for( k in 1:20){
  yt = sample(rs[[k]],1000,replace=TRUE)
  #yt = sample(yt,1000,replace=TRUE)
    for( i in 1:length(g)){
        fit<-fitness(g=g[i],s1[k],s2[k],s3[k],rs=yt)
      #logfit[i]<-log(fit)
      g.mat[,i] <- fit
    }
  g.sites[[k]] <- g.mat
}

plot(seq(0,1,by=0.01),apply(g.mat,2,gm_mean))

gmean <- function(x){apply(x,2,gm_mean)}
site.optima<-lapply(g.sites,gmean)

maxfun <- function(x){g[which(x %in% max(x))]}
optima<-unlist(lapply(site.optima,maxfun))


plot(optima,g1,xlim=c(0,1),ylim=c(0,1),pch=16,
      cex.lab = 1.5, cex.axis = 1.5,
     xlab="Predicted germination fraction",
     ylab="Observed germination fraction")
# Now, define a custom axis
abline(a=0,b=1,lty='dashed')
## Cohen 1966

names<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site) 
df <- data.frame(names,optima,g1)

library(ggrepel)

g1 <- ggplot(df,aes(x=optima,y=g1,label=site)) +
  geom_abline(intercept=0,slope=1,alpha=.5) +
  geom_point() +
  geom_text_repel(size=3,color="black") +
  theme_bw() + xlim(c(0,1)) + ylim(c(0,1)) +
  xlab("Predicted germination fraction") +
  ylab("Observed germination fraction")

ggsave(filename="~/Dropbox/clarkiaSeedBanks/products/figures/obs_pred_germ-labeled.pdf",
       plot=g1,width=6,height=6)



optimaCalc<-function(s1,s2,s3,rs){
  
  g = seq(0,1,by=.01)
  fit = c()
  g.mat = matrix(NA,nrow=100,ncol=length(g))
  g.sites = list()
  yt=c()
  
  for( k in 1:20){
    yt = rs[[k]][sample(7,100,replace=TRUE), ]
#    y = rs[[k]][sample(30000,1000), ]
  #  year = sample(1:(dim(y)[2]),1000,replace=TRUE)
    # for(i in 1:1000){
    #   yt[i] <- y[i,year[i]]
    # }
    for( i in 1:length(g)){
      for( j in 1:100){
        fit[j]<-fitness(g=g[i],s1[k],s2[k],s3[k],rs=sample(yt[i,],1))
        #logfit[i]<-log(fit)
      }
      g.mat[,i] <- fit
    }
    g.sites[[k]] <- g.mat
  }
  
  gmean <- function(x){apply(x,2,gm_mean)}
  site.optima<-lapply(g.sites,gmean)
  
  maxfun <- function(x){g[which(x %in% max(x))]}
  optima<-unlist(lapply(site.optima,maxfun))
}

# uncomment to re-estimate mean.optima
# need to save mean.optima to create figure
trials <- replicate(
  50,
  optimaCalc(s1 = s1, s2 = s2, s3 = s3, rs= rs),
  simplify = FALSE
)
# 
# set.seed(10)
# mean.optima = apply(matrix(unlist(trials),ncol=17,byrow=TRUE),2,mean)
# 
# optimalG<-matrix(unlist(trials),ncol=17,byrow=TRUE)
#save(optimalG,file="~/Dropbox/modelsF2019/output/optimalG")

# load("~/Dropbox/modelsF2019/output/optimalG")

o.mu <- apply(optimalG,2,mean)
o.se <- apply(optimalG,2,sd)/sqrt(50)

pdf(
  "~/Dropbox/clarkiaSeedBanks/products/figures/obs_pred_germ.pdf",
  onefile=TRUE,
  paper="USr",
  height = 7.5, width = 10)

plot(mean.optima,g1,xlim=c(0,1),ylim=c(0,1),pch=16,
     cex.lab = 1.5, cex.axis = 1.5,
     xlab="Predicted germination fraction",
     ylab="Observed germination fraction")
# Now, define a custom axis
abline(a=0,b=1,lty='dashed')
## Cohen 1966
segments(y0=g1,x0=o.mu-o.se,x1=o.mu+o.se)

dev.off()

names<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site) 
df <- data.frame(names,optima,g1)

library(ggrepel)

g1 <- ggplot(df,aes(x=optima,y=g1,label=site)) +
  geom_abline(intercept=0,slope=1,alpha=.5) +
  geom_point() +
  geom_text_repel(size=3,color="black") +
  theme_bw() + xlim(c(0,1)) + ylim(c(0,1)) +
  xlab("Predicted germination fraction") +
  ylab("Observed germination fraction")

ggsave(filename="~/Dropbox/clarkiaSeedBanks/products/figures/obs_pred_germ-labeled.pdf",
       plot=g1,width=6,height=6)
