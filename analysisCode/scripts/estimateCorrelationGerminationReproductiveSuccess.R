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

# read in samples from posterior distributions
load("~/Dropbox/modelsF2019/output/seedbagfit")
load("~/Dropbox/modelsF2019/output/rsVarPosterior")

# need to fix this part?
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

# work with g1
# extract parameters for analysis
posterior.g1<-MCMCchains(zc,params = "alphaG1")

# place on (0,1) scale
probability.g1 <- boot::inv.logit(posterior.g1)

# exclude problem sites
g1<-g1[,-probs]

# calculate the 95% credible interval and HPDI for g1
CI.g1 <- apply(probability.g1,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
HPDI.g1 <- apply(probability.g1,2,FUN = function(x) hdi(x, .95))

# RS
d<-lapply(rsEstimates,temporal_variance,fun=gsd)
reproductiveSuccess<-matrix(unlist(d), ncol = 20, byrow = FALSE)
# exclude problem sites
reproductiveSuccess<-reproductiveSuccess[,-probs]

#
CI.reproductiveSuccess <- apply(reproductiveSuccess,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
HPDI.reproductiveSuccess <- apply(reproductiveSuccess,2,FUN = function(x) hdi(x, .95))

g1PosteriorSummary<-data.frame(t(CI.g1))
names(g1PosteriorSummary) <- c("lo.g1","med.g1","hi.g1")
reproductiveSucessPosteriorSummary<-data.frame(t(CI.reproductiveSuccess))
names(reproductiveSucessPosteriorSummary) <- c("lo.rs","med.rs","hi.rs")

# empty vector for the correlation
posterior.correlation<-c()
n.iter = dim(postrior.g1)[1]

# calculate correlation for each draw from the posterior
for(i in 1:n.iter){
  posterior.correlation[i]<-cor(probability.g1[i,],probability.survival[i,])
}

# calculate the 95% credible interval and HPDI for the correlation
CI.correlation <- quantile(posterior.correlation, c(.025, .5, .975))
HPDI.correlation <- hdi(posterior.correlation, .95)

pdf(
  "~/Dropbox/clarkiaSeedBanks/products/figures/germ_rs_correlation.pdf",
  onefile=TRUE,
  paper="USr",
  height = 7.5, width = 10)

plot(x = reproductiveSucessPosteriorSummary$med.rs,
     y = g1PosteriorSummary$med.g1,
     xlim=c(0,20),ylim=c(0,.3),
     pch=16,
     ylab = "Mean germination probability",
     xlab = "Geometric SD of fitness",
     xaxt='n', cex.lab = 1.5, cex.axis = 1.5)
# Now, define a custom axis
axis(side = 1, at=c(0,2,4,6,8,10,12,14,16,18,20),cex.axis=1.5)

segments(x0=reproductiveSucessPosteriorSummary$med.rs,
         y0=g1PosteriorSummary$lo.g1,
         y1=g1PosteriorSummary$hi.g1)
text(x=15,y=.275,
     paste0("Pearson's r=",round(cor.BCI[2],2)),
     cex=1.5)

hist(posterior.correlation,breaks = 50, main = "", xlab = "", xlim = c(-1, 1),
     freq = FALSE, col = "azure1", cex.lab = 1.5,cex.axis=1.5)

title(xlab="Correlation of germination and geometric SD of fitness \n (Pearson's r)", line=4, cex.lab=1.5)

abline(v=CI.correlation[c(1,3)],lty='dashed',lwd='2')
abline(v=CI.correlation[2],lty='solid',lwd='2',col='red')

dev.off()
