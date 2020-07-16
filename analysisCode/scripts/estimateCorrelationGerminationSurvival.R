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

# -------------------------------------------------------------------
# functions to analyze data
# -------------------------------------------------------------------

cols_fun <- function(x,fun=var){
  apply(x,2,fun)
}

# -------------------------------------------------------------------
# -------------------------------------------------------------------

# read in samples from posterior distributions
# load("~/Dropbox/modelsF2019/output/seedbagfit")
zc <- readRDS("~/Dropbox/dataLibrary/posteriors/belowgroundSamplesAllYears.RDS")

# extract parameters for analysis
posterior.g1<-MCMCchains(zc,params = "g1")
posterior.s2<-MCMCchains(zc,params = "s2")
posterior.s3<-MCMCchains(zc,params = "s3")

# read nSites and nYears from data file
nSites <- 20
nYears <- 3

# get number of iterations from matrix of posterior
n.iter = dim(posterior.g1)[1]

# create empty matrix
probability.survival <- matrix(NA,nrow=n.iter,ncol=nSites)

# calculate product of s2 and s3 for each site
  for(i in 1:nSites){
    probability.survival[,i]<-posterior.s2[,i]*posterior.s3[,i]
  }

# create empty vector for the correlation
posterior.correlation<-c()

# calculate correlation for each draw from the posterior
for(i in 1:n.iter){
  posterior.correlation[i]<-cor(posterior.g1[i,],probability.survival[i,])
}

# calculate the 95% credible interval and HPDI for g1
CI.g1 <- apply(posterior.g1,2,FUN = function(x) quantile(x, c(.25, .5, .75)))
HPDI.g1 <- apply(posterior.g1,2,FUN = function(x) hdi(x, .95))

# calculate the 95% credible interval and HPDI for probability of survival (s2*s3)
CI.survival <- apply(probability.survival,2,FUN = function(x) quantile(x, c(.25, .5, .75)))
HPDI.survival <- apply(probability.survival,2,FUN = function(x) hdi(x, .95))

# put medians and credible intervals into data frame
g1PosteriorSummary <- data.frame(t(CI.g1))
names(g1PosteriorSummary) <- c("lo.g1","med.g1","hi.g1")

survivalPosteriorSummary<-data.frame(t(CI.survival))
names(survivalPosteriorSummary) <- c("lo.surv","med.surv","hi.surv")

# calculate the 95% credible interval and HPDI for the correlation
CI.correlation <- quantile(posterior.correlation, c(.25, .5, .75))
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

names<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site) 
df <- data.frame(names,survivalPosteriorSummary,g1PosteriorSummary)


library(ggrepel)

g1 <- ggplot(df,aes(x=med.surv,y=med.g1,label=site)) +
  geom_abline(intercept=0,slope=1,alpha=.5) +
  geom_point() +
  geom_text_repel(size=3,color="black") +
  theme_bw() + xlim(c(0,1)) + ylim(c(0,1)) +
  xlab("Probability of seed survival [P(S)]") +
  ylab("Mean germination probability [P(G)]")

ggsave(filename="~/Dropbox/clarkiaSeedBanks/products/figures/germ_surv_correlation-labeled.pdf",
       plot=g1,width=6,height=6)


