# -------------------------------------------------------------------
# Analysis of correlation between germination and RS
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


f<-function(x="alphaS1"){
  chain<-MCMCchains(zc,params = x)
  p<-boot::inv.logit(chain)
  BCI <- t(apply(p,2,FUN = function(x) quantile(x, c(.025, .5, .975))))
  return(BCI)
}

# -------------------------------------------------------------------
# -------------------------------------------------------------------

load("~/Dropbox/modelsF2019/output/seedbagfit")
load("~/Dropbox/modelsF2019/output/rsVarPosterior")

# Exclude sites with very high temporal variance (orders of magnitude higher)

d<-lapply(rsEstimates,temporal_variance,fun=gsd)
d2<-matrix(unlist(d), ncol = 20, byrow = FALSE)

BCI <- apply(d2,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
HPDI <- apply(d2,2,FUN = function(x) hdi(x, .95))

df<-data.frame(t(BCI))
names(df) <- c('lo','med','hi')
df<-cbind(df,siteNo=seq(1:20))
df.exclude<-df %>% dplyr::filter(med>20)

# problem sites index
probs <- df.exclude$siteNo

# G1
g1<-boot::inv.logit(MCMCchains(zc,params = "alphaG1"))
# exclude problem sites
g1<-g1[,-probs]

g1.BCI <- apply(g1,2,FUN = function(x) quantile(x, c(.05, .5, .95)))
g1.HPDI <- apply(g1,2,FUN = function(x) hdi(x, .95))

# RS
d<-lapply(rsEstimates,temporal_variance,fun=gsd)
rs<-matrix(unlist(d), ncol = 20, byrow = FALSE)
# exclude problem sites
rs<-rs[,-probs]

rs.BCI <- apply(rs,2,FUN = function(x) quantile(x, c(.05, .5, .95)))
rs.HPDI <- apply(rs,2,FUN = function(x) hdi(x, .95))

g1.df<-data.frame(t(g1.BCI))
names(g1.df) <- c("lo.g1","med.g1","hi.g1")
rs.df<-data.frame(t(rs.BCI))
names(rs.df) <- c("lo.rs","med.rs","hi.rs") 

cor.post<-c()
for(i in 1:30000){
  cor.post[i]<-cor(g1[i,],rs[i,])
}

cor.BCI <- quantile(cor.post, c(.025, .5, .975))
cor.HPDI <- hdi(cor.post, .95)

pdf(
  "~/Dropbox/clarkiaSeedBanks/products/figures/germ_rs_correlation.pdf",
  onefile=TRUE,
  paper="USr",
  height = 7.5, width = 10)

plot(rs.df$med.rs,g1.df$med.g1,xlim=c(0,20),ylim=c(0,.3),
     pch=16, 
     ylab = "Mean germination probability", 
     xlab = "Geometric SD of fitness", xaxt='n', cex.lab = 1.5, cex.axis = 1.5)
# Now, define a custom axis
axis(side = 1, at=c(0,2,4,6,8,10,12,14,16,18,20),cex.axis=1.5)
#segments(x0=rs.df$lo.rs, x1=rs.df$hi.rs, y0=g1.df$med.g1)
segments(x0=rs.df$med.rs, y0=g1.df$lo.g1, y1=g1.df$hi.g1)
text(x=15,y=.275,paste0("Pearson's r=",round(cor.BCI[2],2)),cex=1.5)

hist(cor.post,breaks = 50, main = "", xlab = "", xlim = c(-1, 1), 
     freq = FALSE, col = "azure1", cex.lab = 1.5,cex.axis=1.5)

title(xlab="Correlation of germination and geometric SD of fitness \n (Pearson's r)", line=4, cex.lab=1.5)

abline(v=cor.HPDI,lty='dashed',lwd='2')
abline(v=cor.BCI[2],lty='solid',lwd='2',col='red')

dev.off()

