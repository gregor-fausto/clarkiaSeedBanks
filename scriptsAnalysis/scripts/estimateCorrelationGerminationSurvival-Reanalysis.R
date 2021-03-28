# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Analysis of correlation between g1 and s2s3
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
# read in samples from posterior distributions
# -------------------------------------------------------------------

g1 <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/g1-pop.RDS")
s2 <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/s2-pop.RDS")
s3 <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/s3-pop.RDS")


# -------------------------------------------------------------------
# Uncomment for histograms of posteriors
# -------------------------------------------------------------------
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1    )

for(i in 1:20){
  hist(g1[,i],breaks=100,main='');
  abline(v=median(g1[,i]),lwd=2,col='red');
  abline(v=quantile(g1[,i],c(.025, .5, .975)),lwd=1,col='orange',lty='dotted');
  abline(v=hdi(g1[,i],c(.95)),lwd=1,col='blue',lty='dotted');
}

for(i in 1:20){
  hist(s2[,i]*s3[,i],breaks=100,main='');
  abline(v=median(s2[,i]*s3[,i]),lwd=2,col='red');
  abline(v=quantile(s2[,i]*s3[,i],c(.025, .5, .975)),lwd=1,col='orange',lty='dotted');
  abline(v=hdi(s2[,i]*s3[,i],c(.95)),lwd=1,col='blue',lty='dotted');
}


# -------------------------------------------------------------------
# Analyze correlation
# -------------------------------------------------------------------

n.iter = dim(g1)[1]

probability.survival = s2*s3
probability.g1 = g1

# create empty vector for the correlation
posterior.correlation<-c()

# calculate correlation for each draw from the posterior
for(i in 1:n.iter){
  posterior.correlation[i]<-cor(probability.g1[i,],probability.survival[i,])
}

# calculate the 95% credible interval and HPDI for g1
CI.g1 <- apply(probability.g1,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
HPDI.g1 <- apply(probability.g1,2,FUN = function(x) hdi(x, .95))

# calculate the 95% credible interval and HPDI for probability of survival (s2*s3)
CI.survival <- apply(probability.survival,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
HPDI.survival <- apply(probability.survival,2,FUN = function(x) hdi(x, .95))

# put medians and credible intervals into data frame
g1PosteriorSummary <- data.frame(cbind(t(HPDI.g1),CI.g1[2,]))
names(g1PosteriorSummary) <- c("lo.g1","hi.g1","med.g1")

survivalPosteriorSummary<-data.frame(cbind(t(HPDI.survival),CI.survival[2,]))
names(survivalPosteriorSummary) <- c("lo.surv","hi.surv","med.surv")

# calculate the 95% credible interval and HPDI for the correlation
CI.correlation <- quantile(posterior.correlation, c(.025, .5, .975))
HPDI.correlation <- hdi(posterior.correlation, .95)
HPDI.correlation.lo <- hdi(posterior.correlation, c(.5))

pdf(
 "~/Dropbox/clarkiaSeedBanks/products/figures/analysis/correlation-germ-surv.pdf",
 height = 8, width = 6)

par(mar=c(4,4,2,1))
par(fig=c(0,10,4,10)/10)
# plot median of g1 vs. median of s2*s3 with CIs
plot(x = NA,
     y = NA,
     xlim=c(.2,.9),ylim=c(0,.55),
     pch=16, cex = 0.5,
     xlab = "",
     ylab = "",
     xaxt= "n", yaxt="n",
     cex.lab = 1, cex.axis = 1)

segments(x0=survivalPosteriorSummary$lo.surv,x1=survivalPosteriorSummary$hi.surv,
         y0=g1PosteriorSummary$med.g1, lwd=1)
segments(x0=survivalPosteriorSummary$med.surv,
         y0=g1PosteriorSummary$lo.g1, y1=g1PosteriorSummary$hi.g1,
         lwd=1)
points(survivalPosteriorSummary$med.surv,g1PosteriorSummary$med.g1,
       pch=21,col='black',bg='white',cex=1.25)

axis(1, seq(0,1,by=.1),
     labels = seq(0,1,by=.1), las = 1, line = 0,
     col = NA, col.ticks = 1, cex.axis = 1)
axis(2, seq(0,1,by=.1),
     labels = seq(0,1,by=.1), las = 1, line = 0,
     col = NA, col.ticks = 1, cex.axis = 1)
mtext("Germination probability",
      side=2,line=2.5,adj=.5,col='black',cex=1)
mtext("Survival probability",
      side=1,line=2,adj=.5,col='black',cex=1)

text(x=.1,y=.9,
     paste0("Pearson's r=",round(CI.correlation[2],2)),
     cex=1)
#abline(a=0,b=1)

par(fig=c(0,10,0,4.5)/10)
par(new=T)
# plot posterior of correlation coefficient
hist(posterior.correlation,breaks = 50, 
     main = "", xlab = "", ylab='', xaxt='n',yaxt='n', 
     xlim = c(-1, 1), ylim=c(0,2.5),
     freq = FALSE, col = "gray75", border='white',
     cex.lab = 1.25,cex.axis=1.5)

# as in Duskey dissertation
segments(x0=HPDI.correlation[1],x1=HPDI.correlation[2],y0=2.4,lwd=1.25)
#segments(x0=CI.correlation[2],y0=2.3,y1=2.5,lwd=1.5)
segments(x0=HPDI.correlation.lo[1],x1=HPDI.correlation.lo[2],y0=2.4,lwd=2.5)
points(x=median(posterior.correlation),y=2.4,pch=21,col='black',bg='white',cex=1.25)
#segments(x0=CI.correlation[2],y0=2.4,y1=0,lwd=2,lty='dotted')

axis(1, seq(-2,2,by=.2),
     labels = seq(-2,2,by=.2), las = 1, line = 0,
     col = NA, col.ticks = 1, cex.axis = 1)
axis(2, seq(0,3,by=.4),
     labels = seq(0,3,by=.4), las = 1, line = -.5,
     col = NA, col.ticks = 1, cex.axis = 1)
segments(x0=-1,x1=1,y0=-.1,lwd=2)
segments(x0=-1.05,y0=0,y1=2.4,lwd=1.5)
mtext("Density",
      side=2,line=2.5,adj=.5,col='black',cex=1)
mtext("Correlation of germination and survival",
      side=1,line=2,adj=.5,col='black',cex=1)

dev.off()

# -------------------------------------------------------------------
# Commented code below is for labeling populations
# -------------------------------------------------------------------
# 
# names<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
#   dplyr::select(site) 
# df <- data.frame(names,survivalPosteriorSummary,g1PosteriorSummary)
# ciPosteriorSummary <- data.frame(t(CI.correlation))
# 
# 
# library(ggrepel)
# 
# g1 <- ggplot(df,aes(x=med.surv,y=med.g1,label=site)) +
#   geom_point() +
#   geom_text_repel(size=3,color="black") +
#   theme_bw()+
#   xlab("Probability of seed survival [P(S)]") +
#   ylab("Mean germination probability [P(G)]") +
#   #annotate("text", label =  paste0("Pearson's r=",round(CI.correlation[2],2)), x = .15, y = .29, size = 4) +
#   theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
#  # scale_x_continuous(limits = c(0,.75), expand = c(0, 0)) +
#  # scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# 
# # ggsave(filename="~/Dropbox/clarkiaSeedBanks/products/figures/germ_surv_correlation.pdf",
# #        plot=g1,width=4,height=4)
# 
# 
# c0<-sim(rho=-.75,nsites=20)
# 
# df.sim <- data.frame(df$site,c0)
# g1.hypothesis<-ggplot(df.sim,aes(x=X,y=V2)) +
#   geom_point(color='gray') +
#   #geom_text_repel(size=3,color="black") +
#   theme_bw()+
#   xlab("Probability of seed survival [P(S)]") +
#   ylab("Mean germination probability [P(G)]") +
#   annotate("text", label =  paste0("Pearson's r=-.75"), x = .6, y = .29, size = 4,color='gray') +
#   theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
#   scale_x_continuous(limits = c(0,.75), expand = c(0, 0)) +
#   scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# # ggsave(filename="~/Dropbox/clarkiaSeedBanks/products/figures/germ_surv_correlation-hypothesis.pdf",
# #        plot=g1.hypothesis,width=4,height=4)
# 
# g1.blank<-ggplot(df.sim,aes(x=X,y=V2)) +
#   #geom_point(color='gray') +
#   #geom_text_repel(size=3,color="black") +
#   theme_bw()+
#   xlab("Probability of seed survival [P(S)]") +
#   ylab("Mean germination probability [P(G)]") +
#   annotate("text", label =  paste0("Pearson's r=-.75"), x = .6, y = .29, size = 4,color='gray') +
#   theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
#   scale_x_continuous(limits = c(0,.75), expand = c(0, 0)) +
#   scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# # ggsave(filename="~/Dropbox/clarkiaSeedBanks/products/figures/germ_surv_correlation-blank.pdf",
# #        plot=g1.blank,width=4,height=4)
