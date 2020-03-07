# parameter diagnostic figures

rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)

library(MCMCvis)
library(tidybayes)
library(tidyverse)

# directory = "/Users/Gregor/Dropbox/dataLibrary/ratTumors/"
# fileList <- paste0(directory,list.files(directory))
# dirFigures = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/appendix/rat-tumors/")
# 
# rat_tumors <- read.table("~/Desktop/rat-tumors.txt",header=T)
# 
# data = rat_tumors

directory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/singleSite/")
 fileList <- paste0(directory,list.files(directory))
 dirFigures = c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/singleSite/testFigures/")
 dir.create(file.path(dirFigures), showWarnings = FALSE)
 
 data <- readRDS(fileList[[5]])

y <- data$y
n <- data$n
J <- length(y)

parameterization.MeanSamplesize <- readRDS(fileList[[1]])
parameterization.MeanVariance <- readRDS(fileList[[2]])
parameterization.LogitCentered <- readRDS(fileList[[3]])
parameterization.LogitNoncentered <- readRDS(fileList[[4]])

# MCMCsummary(parameterization.MeanSamplesize,params=c("vcount","vfreq"))
# MCMCsummary(parameterization.MeanVariance,params=c("mu","sigma"))
# MCMCsummary(parameterization.LogitCentered,params=c("mu","sigma"))
# MCMCsummary(parameterization.LogitNoncentered,params=c("mu","sigma"))

# par(mfrow=c(2,2))
# hist(MCMCsummary(parameterization.MeanSamplesize,params=c("theta"))$Rhat,xlim=c(.9,1.1),
#      breaks=seq(.9,1.1,by=.01))
# hist(MCMCsummary(parameterization.MeanVariance,params=c("theta"))$Rhat,xlim=c(.9,1.1),
#      breaks=seq(.9,1.1,by=.01))
# hist(MCMCsummary(parameterization.LogitCentered,params=c("theta"))$Rhat,xlim=c(.9,1.1),
#      breaks=seq(.9,1.1,by=.01))
# hist(MCMCsummary(parameterization.LogitNoncentered,params=c("theta"))$Rhat,xlim=c(.9,1.1),
#      breaks=seq(.9,1.1,by=.01))

colorVector = c("black","red","blue","purple")

## Plot effective sample size vs. theta
par(mfrow=c(1,1))
f<-function(x) x[order(x,decreasing=TRUE)]

nEff.MeanSamplesize = f(MCMCsummary(parameterization.MeanSamplesize,params=c("theta"))$n.eff)
nEff.MeanVariance = f(MCMCsummary(parameterization.MeanVariance,params=c("theta"))$n.eff)
nEff.LogitCentered = f(MCMCsummary(parameterization.LogitCentered,params=c("theta"))$n.eff)
nEff.LogitNoncentered = f(MCMCsummary(parameterization.LogitNoncentered,params=c("theta"))$n.eff)
nEff = c(nEff.MeanSamplesize,nEff.MeanVariance,nEff.LogitCentered,nEff.LogitNoncentered)

pdf(paste0(dirFigures,"nEff.pdf"), width=6, height=6)
plot(nEff.MeanSamplesize,
     ylim=c(min(nEff),max(nEff)),type='l',
     col=colorVector[1],
     xlab="theta_i (ordered)",
     ylab="Effective sample size")
lines(nEff.MeanVariance,
      col=colorVector[2])
lines(nEff.LogitCentered,
      col=colorVector[3])
lines(nEff.LogitNoncentered,
      col=colorVector[4])

legend(x = 0,
       y = min(nEff)+1/5*(max(nEff)-min(nEff)),
       legend=c("Mean & sample size","Mean & variance", "Centered logit", "Noncentered logit"),
       lty = c(1,1,1,1),
       col = colorVector,
       cex=.75)
dev.off()

#interesting: i = 8,1,46
# examine correlation in posterior
shapeit <- function(mu, sigma){
  a <- (mu^2 - mu^3 - mu * sigma^2)/sigma^2
  b <- (mu - 2 * mu^2 + mu^3 - sigma^2 + mu * sigma^2)/sigma^2
  shape_ps <- list(a, b)
  return(shape_ps)
}

i = sample(1:dim(data)[1],1)
samples<-sample(1:dim(vfreq<-MCMCchains(parameterization.MeanSamplesize,params="vfreq")
)[1],1000)

# Parameterization via mean and sample size 
params.MeanSamplesize <- list()
vfreq<-MCMCchains(parameterization.MeanSamplesize,params="vfreq")
vcount<-MCMCchains(parameterization.MeanSamplesize,params="vcount")
theta_i<-MCMCchains(parameterization.MeanSamplesize,params="theta")[,i]

alpha <- vfreq*vcount
beta <- (1-vfreq)*vcount

params.MeanSamplesize <- list(vfreq = vfreq,vcount = vcount,
                              theta_i = theta_i, alpha = alpha,
                              beta = beta)

# parameterization 2
params.MeanVariance <- list()

mu<-MCMCchains(parameterization.MeanVariance,params="mu")
sigma<-MCMCchains(parameterization.MeanVariance,params="sigma")
theta_i<-MCMCchains(parameterization.MeanVariance,params="theta")[,i]

betaParms<-shapeit(mu,sigma)

params.MeanVariance <- list(mu = mu,sigma = sigma,
                              theta_i = theta_i, 
                              alpha = betaParms[[1]],beta = betaParms[[2]])

# parameterization 3
params.LogitCentered <- list()

mu<-MCMCchains(parameterization.LogitCentered,params="mu")
sigma<-MCMCchains(parameterization.LogitCentered,params="sigma")
alpha_i<-MCMCchains(parameterization.LogitCentered,params="alpha")[,1]

betaParms<-shapeit(mu,sigma)

params.LogitCentered <- list(mu = mu,sigma = sigma,
                            alpha_i = alpha_i)

# parameterization 4
params.LogitNoncentered <- list()

mu<-MCMCchains(parameterization.LogitNoncentered,params="mu")
sigma<-MCMCchains(parameterization.LogitNoncentered,params="sigma")
alpha.std_i<-MCMCchains(parameterization.LogitNoncentered,params="alpha.std")[,1]

betaParms<-shapeit(mu,sigma)

params.LogitNoncentered <- list(mu = mu,sigma = sigma,
                             alpha.std_i = alpha.std_i)

## Figure for hyperparameters
pdf(paste0(dirFigures,"hyperparameterCorrelations.pdf"), width=8, height=6)

par(mfrow=c(2,2))

a = params.MeanSamplesize$alpha[samples]
b = params.MeanSamplesize$beta[samples]

plot(log(a/b),log(a+b),
     pch=16,cex=.5,
     col = colorVector[1])

a = params.MeanVariance$alpha[samples]
b = params.MeanVariance$beta[samples]

plot(log(a/b),log(a+b),
     pch=16,cex=.5,
     col = colorVector[2])

plot(params.LogitCentered$mu[samples],params.LogitCentered$sigma[samples],
     pch=16,cex=.5,
     col = colorVector[3])

plot(params.LogitNoncentered$mu[samples],params.LogitNoncentered$sigma[samples],
     pch=16,cex=.5,
     col = colorVector[4])
dev.off()

## Figure for correlation of variance and probability
pdf(paste0(dirFigures,"funnelCheck.pdf"), width=8, height=6)
par(mfrow=c(2,2))
a = params.MeanSamplesize$alpha[samples]
b = params.MeanSamplesize$beta[samples]
theta_i = params.MeanSamplesize$theta_i[samples]

plot(boot::logit(theta_i),log(a+b),
     pch=16,cex=.5,
     col = colorVector[1])

a = params.MeanVariance$alpha[samples]
b = params.MeanVariance$beta[samples]
theta_i = params.MeanVariance$theta_i[samples]

plot(boot::logit(theta_i),log(a+b),
     pch=16,cex=.5,
     col = colorVector[2])

plot(params.LogitCentered$alpha_i[samples],
     log(params.LogitCentered$sigma[samples]),
     pch=16,cex=.5,
     col = colorVector[3])

plot(params.LogitNoncentered$alpha.std_i[samples],
     log(params.LogitNoncentered$sigma[samples]),
     pch=16,cex=.5,
     col = colorVector[4])
dev.off()

# shrinkage
data$freq <- y/n

summary.MeanSamplesize<-apply(MCMCchains(parameterization.MeanSamplesize,params="theta"),2,quantile, probs = c(.05,.5,.95))
summary.MeanVariance<-apply(MCMCchains(parameterization.MeanVariance,params="theta"),2,quantile, probs = c(.05,.5,.95))
summary.LogitCentered<-apply(MCMCchains(parameterization.LogitCentered,params="theta"),2,quantile, probs = c(.05,.5,.95))
summary.LogitNoncentered<-apply(MCMCchains(parameterization.LogitNoncentered,params="theta"),2,quantile, probs = c(.05,.5,.95))

summaryList<-list(summary.MeanSamplesize,summary.MeanVariance,summary.LogitCentered,summary.LogitNoncentered)

pdf(paste0(dirFigures,"shrinkage.pdf"), width=8, height=6)

par(mfrow=c(2,2))
for(i in 1:4){
plot(x = data$freq,
     y = t(summaryList[[i]])[,2],
     xlab= "Tumor frequency in sample",
     ylab="MCMC estimate and 95% CI for tumor probability",
     pch=16,cex=.5,
     xlim=c(0,1),ylim=c(0,1))
segments(x0=data$freq,
         y0= t(summaryList[[i]])[,1],
         y1= t(summaryList[[i]])[,3])
abline(a=0,b=1)
}
dev.off()

# expand shrinkage plot
min.test<-list()
for(i in 1:4){
min.test[[i]]<-data$freq-t(summaryList[[i]])[,2]
}

pdf(paste0(dirFigures,"shrinkageSamplesize.pdf"), width=8, height=6)
par(mfrow=c(1,1))

plot(data$n,
     data$freq-t(summaryList[[i]])[,2],
     ylim=c(min(unlist(min.test)),max(unlist(min.test))),
     xlab="Number of rats in study",
     ylab="Amount of shrinkage",type='n')
abline(h=0)
for(i in 1:4){
  points(data$n,
       data$freq-t(summaryList[[i]])[,2],
       col=colorVector[i],pch=1)
}
dev.off()


# compare direct parameterization and log-odds

# the posterior means for theta are slightly more shrunk towards the population mean
# in this example for log-odds param than for direct
# this is supposedly a feature of the log-odds parameterization
# https://mc-stan.org/users/documentation/case-studies/pool-binary-trials.html
# the log-odds model is supposedly doing a bit more pooling
# For example, Roberto Clemente, the top performing player in the first 45 at bats, has an estimated 35%, 25%, and 18% chance of being the best player according to the no-pooling, partial pooling, and partial pooling of log odds models; the tenth best performing player, Ron Swaboda, has an estimated 1%, 2% and 3% chance according to the same models. That is, the probability of being best goes down for high performing players with more pooling, whereas it goes up for below-average players.
par(mfrow=c(1,1))
plot(x = data$freq,
     y = t(summaryList[[1]])[,2],
     xlab= "Tumor frequency in sample",
     ylab="MCMC estimate and 95% CI for tumor probability",
     pch=16,cex=.5,
     xlim=c(0,.35),ylim=c(0,.35),
     col = rgb( red = .5, green = 0, blue = 0.5, alpha = 0.5),type='n')
abline(a=0,b=1)
abline(h = sum(y)/sum(n),lty='dotted')

points(x = data$freq,
       y = t(summaryList[[1]])[,2],
       pch=16,cex=.5,
       col = rgb( red = .5, green = 0, blue = 0.5, alpha = 0.5))
segments(x0=data$freq,
         y0= t(summaryList[[1]])[,1],
         y1= t(summaryList[[1]])[,3],
         col = rgb( red = .5, green = 0, blue = 0.5, alpha = 0.5))

points(x = data$freq,
       y = t(summaryList[[3]])[,2],
       pch=16,cex=.5,
       col = rgb( red = 0, green = 0, blue = 0, alpha = 0.5))
segments(x0=data$freq,
         y0= t(summaryList[[3]])[,1],
         y1= t(summaryList[[3]])[,3],
         col = rgb( red = 0, green = 0, blue = 0, alpha = 0.5))


# compare population level estimates of mean probability of success

betaParms=list(params.MeanSamplesize$alpha,params.MeanSamplesize$beta)
mean.1 = betaParms[[1]]/(betaParms[[1]]+betaParms[[2]])
betaParms<-shapeit(params.MeanVariance$mu,params.MeanVariance$sigma)
mean.2 = betaParms[[1]]/(betaParms[[1]]+betaParms[[2]])
mean.3=boot::inv.logit(params.LogitCentered$mu)
mean.4=boot::inv.logit(params.LogitNoncentered$mu)

means<-cbind(mean.1,mean.2,mean.3,mean.4)
max.function<-function(x) max(density(x)$y)
density.max<-max(apply(means,2,max.function))
means<-t(apply(means,2,quantile,probs=c(0.025,.5,.975)))


pdf(paste0(dirFigures,"posteriorMean.pdf"), width=6, height=6)
plot(density(mean.1),type='n',
     xlim=c(0,1),
     ylim=c(0,density.max),
     main="Posterior distributions for mean probability",
     xlab="Mean probability")
lines(density(mean.1))
lines(density(mean.2),col="red")
lines(density(mean.3),col="blue")
lines(density(mean.4),col="purple")

points(means[,2],y=c(35,38,41,44),ylim=c(0,1),
       pch=16,col=c('black','red','blue','purple'))
segments(x0=means[,1],y0=c(35,38,41,44),x1=means[,3],
         col=c('black','red','blue','purple'))
abline(v=sum(y)/sum(n),lty='dashed')
dev.off()

##
# compare marginalized probabilities
p.1<-MCMCchains(parameterization.MeanSamplesize,params="p")
p.2<-MCMCchains(parameterization.MeanVariance,params="p")
p.3<-MCMCchains(parameterization.LogitCentered,params="p")
p.4<-MCMCchains(parameterization.LogitNoncentered,params="p")

ps<-cbind(p.1,p.2,p.3,p.4)
max.function<-function(x) max(density(x)$y)
density.max<-max(apply(ps,2,max.function))
ps<-t(apply(ps,2,quantile,probs=c(0.025,.5,.975)))
pdf(paste0(dirFigures,"posteriorProbability.pdf"), width=6, height=6)

plot(density(p.1),type='n',
     xlim=c(0,1),
     ylim=c(0,45),
     main="Posterior distributions for population-level probability",
     xlab="Probability")
lines(density(p.1))
lines(density(p.2),col="red")
lines(density(p.3),col="blue")
lines(density(p.4),col="purple")

points(ps[,2],y=c(6,7,8,9),ylim=c(0,1),
       pch=16,col=c('black','red','blue','purple'))
segments(x0=ps[,1],y0=c(6,7,8,9),x1=ps[,3],
         col=c('black','red','blue','purple'))
abline(v=sum(y)/sum(n),lty='dashed')

legend(x = .65,
       y = 9 ,
       legend=c("Mean & sample size","Mean & variance", "Centered logit", "Noncentered logit"),
       lty = c(1,1,1,1),
       col = c("black","red","blue","purple"),
       cex=.75)
dev.off()
##
sigmaBeta = function(alpha, beta) (alpha*beta)/(((alpha+beta)^2)*(alpha+beta+1))

sigma2.1 = sigmaBeta(betaParams[[1]],betaParams[[2]])
mean(sigma2.1)
p<-MCMCchains(parameterization.MeanSamplesize,params="p")
var(p)
# parameterization 2

sigma2.2 = sigmaBeta(betaParms[[1]],betaParms[[2]])
p<-MCMCchains(parameterization.MeanVariance,params="p")
var(p);mean(sigma2.2)
# parameterization 3


sigma.3=params.LogitCentered$sigma
p<-MCMCchains(parameterization.LogitCentered,params="p")
var(p)
# parameterization 4

sigma.4=params.LogitCentered$sigma
p<-MCMCchains(parameterization.LogitNoncentered,params="p")
var(p)





