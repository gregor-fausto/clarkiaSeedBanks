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
# -------------------------------------------------------------------
# read in samples from posterior distributions
g1 <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/g1-pop.RDS")
s2 <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/s2-pop.RDS")
s3 <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/s3-pop.RDS")

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


n.iter = dim(g1)[1]
# # extract parameters for analysis
# posterior.g1<-MCMCchains(zc,params = "g1")
# posterior.s2<-MCMCchains(zc,params = "s2")
# posterior.s3<-MCMCchains(zc,params = "s3")

# read nSites and nYears from data file
# nSites <- 20
# nYears <- 3
# 
# # get number of iterations from matrix of posterior
# n.iter = dim(posterior.g1)[1]

# create empty matrix
# probability.survival <- matrix(NA,nrow=n.iter,ncol=nSites)

# calculate product of s2 and s3 for each site
  # for(i in 1:nSites){
  #   probability.survival[,i]<-posterior.s2[,i]*posterior.s3[,i]
  # }
probability.survival = s2*s3
posterior.g1 = g1
# create empty vector for the correlation
posterior.correlation<-c()

# calculate correlation for each draw from the posterior
for(i in 1:n.iter){
  posterior.correlation[i]<-cor(posterior.g1[i,],probability.survival[i,])
}

# calculate the 95% credible interval and HPDI for g1
CI.g1 <- apply(posterior.g1,2,FUN = function(x) quantile(x, c(.0255, .5, .975)))
HPDI.g1 <- apply(posterior.g1,2,FUN = function(x) hdi(x, .95))

# calculate the 95% credible interval and HPDI for probability of survival (s2*s3)
CI.survival <- apply(probability.survival,2,FUN = function(x) quantile(x, c(.0255, .5, .975)))
HPDI.survival <- apply(probability.survival,2,FUN = function(x) hdi(x, .95))

# put medians and credible intervals into data frame
g1PosteriorSummary <- data.frame(t(CI.g1))
names(g1PosteriorSummary) <- c("lo.g1","med.g1","hi.g1")

survivalPosteriorSummary<-data.frame(t(CI.survival))
names(survivalPosteriorSummary) <- c("lo.surv","med.surv","hi.surv")

# calculate the 95% credible interval and HPDI for the correlation
CI.correlation <- quantile(posterior.correlation, c(.025, .5, .975))
HPDI.correlation <- hdi(posterior.correlation, .95)


# pdf(
#  "~/Dropbox/clarkiaSeedBanks/products/figures/germ_surv_correlation.pdf",
#  onefile=TRUE,
#  paper="USr",
#  height = 7.5, width = 10)

# change par
#mar.default <- c(5,4,4,2) + 0.1
#par(mar = mar.default + c(0, 4, 0, 0))

par(mar=c(2,2,2,1))
par(fig=c(0,10,3,10)/10)
plot(1,main=1)
# par(fig=c(7,10,6,10)/10)
# par(new=T)
# plot(2,main=2)
par(fig=c(0,10,0,3)/10)
par(new=T)
plot(3,main=3)
# par(fig=c(7,10,0,6)/10)
# par(new=T)
# plot(4,main=4)

par(mar=c(4,4,2,1))
par(fig=c(0,10,4,10)/10)
# plot median of g1 vs. median of survival with CIs
plot(x = NA,
     y = NA,
     xlim=c(.1,.9),ylim=c(0,.6),
     pch=16, cex = 0.5,
     xlab = "",
     ylab = "",
     xaxt= "n", yaxt="n",
     cex.lab = 1, cex.axis = 1)

segments(x0=survivalPosteriorSummary$lo.surv,x1=survivalPosteriorSummary$hi.surv,
         y0=g1PosteriorSummary$med.g1, lwd=.5)
segments(x0=survivalPosteriorSummary$med.surv,
         y0=g1PosteriorSummary$lo.g1, y1=g1PosteriorSummary$hi.g1,
          lwd=.5)
points(survivalPosteriorSummary$med.surv,g1PosteriorSummary$med.g1,
       pch=21,col='black',bg='white',cex=.75)

axis(1, seq(0,1,by=.1),
     labels = seq(0,1,by=.1), las = 1, line = 0,
     col = NA, col.ticks = 1, cex.axis = 1)
axis(2, seq(0,1,by=.1),
     labels = seq(0,1,by=.1), las = 1, line = 0,
     col = NA, col.ticks = 1, cex.axis = 1)
mtext("Germination probability",
      side=2,line=2.5,adj=.5,col='black',cex=1)
mtext("Seed survival probability",
      side=1,line=2,adj=.5,col='black',cex=1)

text(x=.2,y=.58,
     paste0("Pearson's r=",round(CI.correlation[2],2)),
     cex=1)
#abline(a=0,b=1)

par(fig=c(0,10,0,4)/10)
par(new=T)
# plot posterior of correlation coefficient
hist(posterior.correlation,breaks = 100, main = "", xlab = "", xlim = c(-1, 1),
     freq = FALSE, col = "gray75", cex.lab = 1.25,cex.axis=1.5)

title(xlab="Correlation of germination and seed survival \n (Pearson's r)", line=3, cex.lab=1.25)

abline(v=CI.correlation[c(1,3)],lty='dashed',lwd='2')
abline(v=CI.correlation[2],lty='solid',lwd='2',col='black')

#dev.off()

names<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site) 
df <- data.frame(names,survivalPosteriorSummary,g1PosteriorSummary)
ciPosteriorSummary <- data.frame(t(CI.correlation))


library(ggrepel)

g1 <- ggplot(df,aes(x=med.surv,y=med.g1,label=site)) +
  geom_point() +
  geom_text_repel(size=3,color="black") +
  theme_bw()+
  xlab("Probability of seed survival [P(S)]") +
  ylab("Mean germination probability [P(G)]") +
  #annotate("text", label =  paste0("Pearson's r=",round(CI.correlation[2],2)), x = .15, y = .29, size = 4) +
  theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
 # scale_x_continuous(limits = c(0,.75), expand = c(0, 0)) +
 # scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# ggsave(filename="~/Dropbox/clarkiaSeedBanks/products/figures/germ_surv_correlation.pdf",
#        plot=g1,width=4,height=4)


c0<-sim(rho=-.75,nsites=20)

df.sim <- data.frame(df$site,c0)
g1.hypothesis<-ggplot(df.sim,aes(x=X,y=V2)) +
  geom_point(color='gray') +
  #geom_text_repel(size=3,color="black") +
  theme_bw()+
  xlab("Probability of seed survival [P(S)]") +
  ylab("Mean germination probability [P(G)]") +
  annotate("text", label =  paste0("Pearson's r=-.75"), x = .6, y = .29, size = 4,color='gray') +
  theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
  scale_x_continuous(limits = c(0,.75), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ggsave(filename="~/Dropbox/clarkiaSeedBanks/products/figures/germ_surv_correlation-hypothesis.pdf",
#        plot=g1.hypothesis,width=4,height=4)

g1.blank<-ggplot(df.sim,aes(x=X,y=V2)) +
  #geom_point(color='gray') +
  #geom_text_repel(size=3,color="black") +
  theme_bw()+
  xlab("Probability of seed survival [P(S)]") +
  ylab("Mean germination probability [P(G)]") +
  annotate("text", label =  paste0("Pearson's r=-.75"), x = .6, y = .29, size = 4,color='gray') +
  theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
  scale_x_continuous(limits = c(0,.75), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ggsave(filename="~/Dropbox/clarkiaSeedBanks/products/figures/germ_surv_correlation-blank.pdf",
#        plot=g1.blank,width=4,height=4)
