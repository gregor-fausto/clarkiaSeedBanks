# -------------------------------------------------------------------
# Analysis of correlation between g and s2s3
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

cols_fun <- function(x,fun=var){
  apply(x,2,fun)
}

# -------------------------------------------------------------------
# -------------------------------------------------------------------

load("~/Dropbox/modelsF2019/output/seedbagfit")

g1<-boot::inv.logit(MCMCchains(zc,params = "alphaG1"))
s2<-boot::inv.logit(MCMCchains(zc,params = "alphaS2"))
s3<-boot::inv.logit(MCMCchains(zc,params = "alphaS3"))

nsite <- 20
nyear <- 3

survEstimates<-list()
surv<-matrix(NA,nrow=30000,ncol=nsite)
  
  for(k in 1:nsite){
    
    surv[,k]<-s2[,k]*s3[,k]
  }

cor.post<-c()
for(i in 1:30000){
  cor.post[i]<-cor(g1[i,],surv[i,])
}  

g1.BCI <- apply(g1,2,FUN = function(x) quantile(x, c(.05, .5, .95)))
g1.HPDI <- apply(g1,2,FUN = function(x) hdi(x, .95))

surv.BCI <- apply(surv,2,FUN = function(x) quantile(x, c(.05, .5, .95)))
surv.HPDI <- apply(surv,2,FUN = function(x) hdi(x, .95))

g1.df<-data.frame(t(g1.BCI))
names(g1.df) <- c("lo.g1","med.g1","hi.g1")
surv.df<-data.frame(t(surv.BCI))
names(surv.df) <- c("lo.surv","med.surv","hi.surv") 


cor.BCI <- quantile(cor.post, c(.025, .5, .975))
cor.HPDI <- hdi(cor.post, .95)


pdf(
  "~/Dropbox/clarkiaSeedBanks/products/figures/germ_surv_correlation.pdf",
  onefile=TRUE,
  paper="USr",
  height = 7.5, width = 10)

# change par 
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 

plot(surv.df$med.surv,g1.df$med.g1,xlim=c(.5,1),ylim=c(0,.4),
     pch=16, 
     ylab = "Mean germination probability [P(G)]", 
     xlab = "Probability of seed survival [P(S)]",
     cex.lab = 1.5, cex.axis = 1.5)
# Now, define a custom axis
#segments(x0=rs.df$lo.rs, x1=rs.df$hi.rs, y0=g1.df$med.g1)
segments(x0=surv.df$med.surv, y0=g1.df$lo.g1, y1=g1.df$hi.g1)
segments(y0=g1.df$med.g1, x0=surv.df$lo.surv, x1=surv.df$hi.surv)
text(x=.8,y=.375,paste0("Pearson's r=",round(cor.BCI[2],2)),cex=1.5)


hist(cor.post,breaks = 50, main = "", xlab = "", xlim = c(-1, 1), 
     freq = FALSE, col = "azure1", cex.lab = 1.5,cex.axis=1.5)

title(xlab="Correlation of germination and seed survival \n (Pearson's r)", line=4, cex.lab=1.5)


abline(v=cor.HPDI,lty='dashed',lwd='2')
abline(v=cor.BCI[2],lty='solid',lwd='2',col='red')

dev.off()



