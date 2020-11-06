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
  n = length(x[!is.na(x)])
  mu = exp(mean(log(x),na.rm=TRUE))
  y <- exp(sqrt(sum((log(x/mu))^2,na.rm=TRUE)/(n-1)))
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
 # geom_text_repel(size=3,color="black") +
  annotate("text", label =  paste0("Pearson's r=",round(CI.correlation[1],2)), x = 2.5, y = .29, size = 4) +
 theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
  scale_x_continuous(limits = c(.99,8), expand = c(0, 0), breaks = c(1, 3, 5, 7)) +
  scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
  xlab("Geometric SD of reproductive success") +
  ylab("Mean germination probability [P(G)]") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename=  "~/Dropbox/clarkiaSeedBanks/products/figures/germ_rs_correlation-labeled.pdf",
       plot=g1,width=4,height=4)

g1.blank <- ggplot(df,aes(x=var.rs,y=med,label=site)) +
 # geom_point() +
  # geom_text_repel(size=3,color="black") +
 # annotate("text", label =  paste0("Pearson's r=",round(CI.correlation[1],2)), x = 1.5, y = .29, size = 4) +
  theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
  scale_x_continuous(limits = c(.99,8), expand = c(0, 0), breaks = c(1, 3, 5, 7)) +
  scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
  xlab("Geometric SD of reproductive success") +
  ylab("Mean germination probability [P(G)]") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename=  "~/Dropbox/clarkiaSeedBanks/products/figures/germ_rs_correlation-blank.pdf",
       plot=g1.blank,width=4,height=4)


# g1.point <- ggplot(df %>% dplyr::filter(site=="BG"),aes(x=var.rs,y=med,label=site)) +
#    geom_point() +
#   # geom_text_repel(size=3,color="black") +
#   theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
#   scale_x_continuous(limits = c(.99,8), expand = c(0, 0), breaks = c(1, 3, 5, 7)) +
#   scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
#   xlab("Geometric SD of reproductive success") +
#   ylab("Mean germination probability [P(G)]") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# ggsave(filename=  "~/Dropbox/clarkiaSeedBanks/products/figures/germ_rs_correlation-point.pdf",
#        plot=g1.point,width=4,height=4)

matplot(c0[,2],c0[,3:12],
        xlab="Mean site germination probability",
        ylab="Annual mean fitness",
        main=expression(paste("cor(",sigma,",P(G))=-.75")),
        pch=16,col='gray')
plot(c0[1:nsites,1],c0[1:nsites,2],xlab="Variance in Fitness",ylab="Probability of Germination")

set.seed(11)

# COPULA

library(MASS)

nsites = 20
nyears = 10

mu = mu = rnorm(n=nsites,mean=25,sd=0)

sim <- function(rho=0,nsites=20){
  # Defines the  sequence for stochastic trials.
  reps=1000
  
  # a and b are hyperparameters of the gamma distribution 
  # that define both the expected value and variance.   
  a = 6
  b = 1.75
  
  # alpha and beta are hyperparameters of the beta distribution that define both the expected value and
  # variance.  
  alpha =  1
  beta =  1
  
  # Defines the temporal correlation between the two parameters.
  rho = rho
  
  # Generates standard multivariate normal data with correlation structure defined by rho.
  Z <- mvrnorm(n=reps,mu=c(0,0), 
               matrix(data=c(1,rho,rho,1),
                      nrow=2,ncol=2))
  
  # Apply the Normal CDF function to Z to obtain data that is uniform on the interval [0,1], but still correlated.
  U = pnorm(Z)
  
  # x is gamma distributed
  X = qgamma(U[,1],shape=a,rate=b) 
 # X = qunif(U[,1],0,.75) 
  
  # y is beta distributed
  #X <- cbind(X,qbeta(U[,2],shape1=alpha,shape2=beta) )
  X <- cbind(X,qunif(U[,2],0,.3) )
  
  # gamma marginal of multivariate X.
  # hist(X[,1])
  # beta marginal of multivariate X
  # hist(X[,2])
  
  # plot(X[,1],X[,2])
  
  nsites = nsites
  nyears = 10
  X[1:nsites,1:2]
  
  # generate samples for each site with a different variance
  d<-lapply(X[1:nsites,1],rnorm,n=nyears,mean=rnorm(n=nsites,mean=mu,sd=0))
  # check that the variances are appropriately sampled
  cbind(unlist(lapply(d, sd)),X[1:nsites,1])
  
  d<-data.frame(d)
  names(d) <- 1:20
  
  # check variances again
  cbind(apply(d,2, sd),X[1:nsites,1])
  
  dim(d)
  
  dt<-cbind(X[1:nsites,1:2],t(d))
  return(dt)
}



c0<-sim(rho=-.75,nsites=20)

df.sim <- data.frame(df$site,c0)
g1.hypothesis<-ggplot(df.sim,aes(x=X,y=V2)) +
  geom_point(color='gray') +
  # geom_text_repel(size=3,color="black") +
   annotate("text", label =  paste0("Pearson's r=",-.75), x = 6, y = .29, size = 4,color='gray') +
  theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
  scale_x_continuous(limits = c(.99,8), expand = c(0, 0), breaks = c(1, 3, 5, 7)) +
  scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
  xlab("Geometric SD of reproductive success") +
  ylab("Mean germination probability [P(G)]") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename=  "~/Dropbox/clarkiaSeedBanks/products/figures/germ_rs_correlation-hypothesis.pdf",
       plot=g1.hypothesis,width=4,height=4)


g1.point<-ggplot(df.sim[1,],aes(x=X,y=V2)) +
  geom_segment(aes(x=X,y=0,xend=X,yend=V2),color='#E69F00') +
  geom_segment(aes(x=1,y=V2,xend=X,yend=V2),color='#CC79A7') +
  
  geom_point(color='gray',size=5) +
  # geom_text_repel(size=3,color="black") +
 # annotate("text", label =  paste0("Pearson's r=",-.75), x = 6, y = .29, size = 4,color='gray') +
  theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
  scale_x_continuous(limits = c(.99,8), expand = c(0, 0), breaks = c(1, 3, 5, 7)) +
  scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
  xlab("Geometric SD of reproductive success") +
  ylab("Mean germination probability [P(G)]") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename=  "~/Dropbox/clarkiaSeedBanks/products/figures/germ_rs_correlation-point.pdf",
       plot=g1.point,width=4,height=4)


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