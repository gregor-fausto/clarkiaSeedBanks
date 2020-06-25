# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Analysis of correlation between g and s2s3
# -------------------------------------------------------------------
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE,max.print=100000)
# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(MCMCvis)
library(tidyverse)
library(HDInterval)
library(bayesplot)

set.seed(10)
# -------------------------------------------------------------------
# functions to analyze data
# -------------------------------------------------------------------

cols_fun <- function(x,fun=var){
  apply(x,2,fun)
}

# -------------------------------------------------------------------
# -------------------------------------------------------------------

# read in samples from posterior distributions
load("~/Dropbox/modelsF2019/output/seedbagfit")
# readrDS(...)

# extract parameters for analysis
posterior.g1<-MCMCchains(zc,params = "alphaG1")
posterior.s2<-MCMCchains(zc,params = "alphaS2")
posterior.s3<-MCMCchains(zc,params = "alphaS3")

# place on (0,1) scale
probability.g1 <- boot::inv.logit(posterior.g1)
probability.s2 <- boot::inv.logit(posterior.s2)
probability.s3 <- boot::inv.logit(posterior.s3)

# read nSites and nYears from data file
nSites <- 20
nYears <- 3
n.iter = dim(posterior.g1)[1]

# empty matrix 
probability.survival <- matrix(NA,nrow=n.iter,ncol=nSites)

# calculate product of s2 and s3 for each site  
  for(i in 1:nSites){
    probability.survival[,i]<-probability.s2[,i]*probability.s3[,i]
  }

# empty vector for the correlation
posterior.correlation<-c()

# calculate correlation for each draw from the posterior
for(i in 1:n.iter){
  posterior.correlation[i]<-cor(probability.g1[i,],probability.survival[i,])
}  

# calculate the 95% credible interval and HPDI for g1
CI.g1 <- apply(probability.g1,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
HPDI.g1 <- apply(probability.g1,2,FUN = function(x) hdi(x, .95))

# calculate the 95% credible interval and HPDI for probability of survival (s2*s3)
CI.survival <- apply(probability.survival,2,FUN = function(x) quantile(x, c(.05, .5, .95)))
HPDI.survival <- apply(probability.survival,2,FUN = function(x) hdi(x, .95))

# put medians and credible intervals into data frame
g1PosteriorSummary <- data.frame(t(CI.g1))
names(g1PosteriorSummary) <- c("lo.g1","med.g1","hi.g1")
survivalPosteriorSummary<-data.frame(t(CI.survival))
names(survivalPosteriorSummary) <- c("lo.surv","med.surv","hi.surv") 

# calculate the 95% credible interval and HPDI for the correlation
CI.correlation <- quantile(posterior.correlation, c(.025, .5, .975))
HPDI.correlation <- hdi(posterior.correlation, .95)


pdf(
  "~/Dropbox/clarkiaSeedBanks/products/figures/germ_surv_correlation.pdf",
  onefile=TRUE,
  paper="USr",
  height = 7.5, width = 10)

# change par 
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 

# plot median of g1 vs. median of survival with CIs
plot(x = survivalPosteriorSummary$med.surv,
     y = g1PosteriorSummary$med.g1,
     xlim=c(0,1),ylim=c(0,1),
     pch=16, cex = 0.5,
     xlab = "Probability of seed survival [P(S)]",
     ylab = "Mean germination probability [P(G)]", 
     cex.lab = 1, cex.axis = 1)

segments(x0=survivalPosteriorSummary$lo.surv,x1=survivalPosteriorSummary$hi.surv,
         y0=g1PosteriorSummary$med.g1, y1=g1PosteriorSummary$med.g1)
segments(x0=survivalPosteriorSummary$med.surv,x1=survivalPosteriorSummary$med.surv,
         y0=g1PosteriorSummary$lo.g1, y1=g1PosteriorSummary$hi.g1)
text(x=.175,y=.95,
     paste0("Pearson's r=",round(CI.correlation[2],2)),
     cex=1)
abline(a=0,b=1)

# plot posterior of correlation coefficient
hist(posterior.correlation,breaks = 50, main = "", xlab = "", xlim = c(-1, 1), 
     freq = FALSE, col = "azure1", cex.lab = 1.5,cex.axis=1.5)

title(xlab="Correlation of germination and seed survival \n (Pearson's r)", line=4, cex.lab=1.5)

abline(v=CI.correlation[c(1,3)],lty='dashed',lwd='2')
abline(v=CI.correlation[2],lty='solid',lwd='2',col='red')

dev.off()



