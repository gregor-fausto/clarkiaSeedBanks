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
library(rethinking)

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

# sample based calculation of gsd
gsd.am <- function(x){
  n = length(x)
  mu = exp(mean(log(x),na.rm=TRUE))
  y <- exp(sqrt(sum((log(x/mu))^2,na.rm=TRUE)/n))
  return(y)
}

f<-function(x="alphaS1"){
  chain<-MCMCchains(zc,params = x)
  BCI <- t(apply(p,2,FUN = function(x) quantile(x, c(.025, .5, .975))))
  return(BCI)
}

# -------------------------------------------------------------------
# -------------------------------------------------------------------

# read in samples from posterior distributions

# belowground <- readRDS("~/Dropbox/dataLibrary/posteriors/jointInferenceSamples.RDS")
# aboveground <- readRDS("~/Dropbox/clarkiaSeedBanks/products/dataFiles/rsVarFullPosterior.RDS")

g1Summary<-read.csv("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/parameterSummary/g1Summary.csv")[,-1]
rsMedians<-readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/dataFiles/rsMedianEstimates.RDS")


gsdSummary <- rsMedians %>%
  dplyr::filter(!is.na(rs)) %>%
  dplyr::select(site,year,rs) %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(var.rs=gsd.am(rs))

CI.correlation<-cor(g1Summary$med,gsdSummary$var.rs)


# pdf(
#   "~/Dropbox/clarkiaSeedBanks/products/figures/germ_rs_correlation.pdf",
#   onefile=TRUE,
#   paper="USr",
#   height = 7.5, width = 10)

plot(x = gsdSummary$var.rs,
     y = g1Summary$med,
     xlim=c(0,8),ylim=c(0,.5),
     pch=16,
     ylab = "Mean germination probability",
     xlab = "Geometric SD of fitness",
     xaxt='n', cex.lab = 1.5, cex.axis = 1.5)
# Now, define a custom axis
axis(side = 1, at=c(0,2,4,6,8),cex.axis=1.5)

# segments(x0=gsdSummary$var.rs,
#          y0=g1Summary$ci.lo95,
#          y1=g1Summary$ci.hi95)
text(x=2,y=.45,
     paste0("Pearson's r=",round(CI.correlation[1],2)),
     cex=1.5)

#dev.off()


df <- g1Summary %>% dplyr::left_join(gsdSummary,by="site")

library(ggrepel)

g1 <- ggplot(df,aes(x=var.rs,y=med,label=site)) +
  geom_point() +
  geom_text_repel(size=3,color="black") +
  theme_bw() + xlim(c(0,8)) + ylim(c(0,.5)) +
  xlab("Geometric SD of fitness") +
  ylab("Mean germination probability [P(G)]") 

# ggsave(filename=  "~/Dropbox/clarkiaSeedBanks/products/figures/germ_rs_correlation-labeled.pdf",
#        plot=g1,width=4,height=4)

# # Germination
# # extract parameters for analysis
#  posterior.g1<-MCMCchains(belowground,params = "g1")
# 
# # calculate the 95% credible interval and HPDI for g1
# CI.g1 <- apply(posterior.g1,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
# HPDI.g1 <- apply(posterior.g1,2,FUN = function(x) rethinking::HPDI(x, .95))
# 
# ## Reproductive success
# lapply(aboveground,)
# 
# d<-lapply(aboveground,temporal_variance,fun=gsd)
# reproductiveSuccess<-matrix(unlist(d), ncol = 20, byrow = FALSE)
# # # exclude problem sites
# # reproductiveSuccess<-reproductiveSuccess[,-probs]
# 
# # don't currently have the full distribution?
# CI.reproductiveSuccess <- apply(reproductiveSuccess,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
# HPDI.reproductiveSuccess <- apply(reproductiveSuccess,2,FUN = function(x) hdi(x, .95))
# 
# g1PosteriorSummary<-data.frame(t(CI.g1))
# names(g1PosteriorSummary) <- c("lo.g1","med.g1","hi.g1")
# 
# ## need to finish this section
# reproductiveSuccessPosteriorSummary<-data.frame(t(CI.reproductiveSuccess))
# names(reproductiveSuccessPosteriorSummary) <- c("lo.rs","med.rs","hi.rs")
# 
# # calculate correlation for each draw from the posterior
# n.iter=3000
# posterior.correlation<-c()
# 
# for(i in 1:n.iter){
#   posterior.correlation[i]<-cor(posterior.g1[i,],reproductiveSuccess[i,])
# }
# 
# # calculate the 95% credible interval and HPDI for the correlation
# CI.correlation <- quantile(posterior.correlation, c(.025, .5, .975))
# HPDI.correlation <- hdi(posterior.correlation, .95)
# 
# 
#   
# 
# ggsave(filename="~/Dropbox/clarkiaSeedBanks/products/figures/germ_rs_correlation-labeled.pdf",
#        plot=g1,width=6,height=6)
# 
# 
# 
# hist(posterior.correlation,breaks = 50, main = "", xlab = "", xlim = c(-1, 1),
#      freq = FALSE, col = "azure1", cex.lab = 1.5,cex.axis=1.5)
# 
# title(xlab="Correlation of germination and geometric SD of fitness \n (Pearson's r)", line=4, cex.lab=1.5)
# 
# abline(v=CI.correlation[c(1,3)],lty='dashed',lwd='2')
# abline(v=CI.correlation[2],lty='solid',lwd='2',col='red')