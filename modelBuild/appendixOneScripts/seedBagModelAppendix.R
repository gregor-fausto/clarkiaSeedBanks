
rm(list=ls(all=TRUE)) # clear R environment

library(rjags)
library(MCMCvis)
library(tidyverse)
# pass data to list for JAGS

set.seed(sample(1:100,1))
#directory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/library/dataFromWorkflowFile/"
#load(paste0(directory,"seedBagsData.rda"))
modelDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/jagsScriptsBinomialModel/"
dirFigures = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/appendix/appendix-1/")
dir.create(file.path(dirFigures), showWarnings = FALSE)


### Appendix 1 dataset
directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/"
simFiles <- paste0(directory,list.files(directory))

simData <- readRDS(simFiles[[4]])

y = simData$fruitingPlantNumber
n = simData$seedlingNumber

# simData <- seedBags 
# 
# selectedSite<-sample(unique(simData$site),1)
# simData=simData[simData$site==selectedSite,]
# simData=simData[simData$age==2,]

# ## germinated given total intact in January
# y = simData$seedlingJan
# n = simData$totalJan

# ## intact in October given intact but not germinated in January
# y = simData$intactOct
# n = simData$intactJan

### SIMULATION STUDY
# start with 30 plots or 30 bags
# nSamples = 30
# # each bag or plot has 100 seeds at the start
# n = rnbinom(nSamples,size=20,mu=50)
# 
# # chance of success
# p = 0.25
# 
# # response 
# y = rbinom(n = nSamples, size = n, prob = p)
# 
# # length
# J <- length(y)

data = list(
  y = y,
  n = n,
  N = length(y)
)

# priors
mu.prior = 0.5
mu.prior.logit = 0
inits = list(list(vfreq = mu.prior,vcount = 20),list(vfreq = mu.prior,vcount = 20),list(vfreq = mu.prior,vcount = 20))
inits2 = list(list(mu = mu.prior,sigma = .05),list(mu = mu.prior,sigma = .05),list(mu = mu.prior,sigma = .05))
inits3 = list(list(mu = mu.prior.logit, sigma = 5),list(mu = mu.prior.logit, sigma = 10),list(mu = 0, sigma = 20))
inits4 = list(list(mu = mu.prior.logit, sigma = 5),list(mu = mu.prior.logit, sigma = 10),list(mu = mu.prior.logit, sigma = 20))


# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 1000
n.update = 5000
n.iterations = 9000
n.thin = 1

parsToMonitor = c("theta","vfreq","vcount","p")
parsToMonitor2 = c("theta","mu","sigma","p")
parsToMonitor3 = c("theta","alpha","mu","sigma","p")
parsToMonitor4 = c("theta","alpha","alpha.std","mu","sigma","p")

#set.seed(2)
# tuning (n.adapt)
jm1 = jags.model(file=paste0(modelDirectory,"hierarchicalMeanConcentration.R"), inits=inits,
                 n.chains = length(inits), data = data, n.adapt = n.adapt)
jm2 = jags.model(file=paste0(modelDirectory,"hierarchicalMomentMatched.R"), inits=inits2,
                 n.chains = length(inits2), data = data, n.adapt = n.adapt)
jm3 = jags.model(file=paste0(modelDirectory,"hierarchicalLogitCentered.R"), inits=inits3,
                 n.chains = length(inits3), data = data, n.adapt = n.adapt)
jm4 = jags.model(file=paste0(modelDirectory,"hierarchicalLogitNoncentered.R"), inits=inits4,
                 n.chains = length(inits4), data = data, n.adapt = n.adapt)
# burn-in (n.update)
update(jm1, n.iterations = n.update)
update(jm2, n.iterations = n.update)
update(jm3, n.iterations = n.update)
update(jm4, n.iterations = n.update)

# chain (n.iter)
samples.rjags1 = coda.samples(jm1, variable.names = c(parsToMonitor), n.iter = n.iterations, thin = n.thin)
samples.rjags2 = coda.samples(jm2, variable.names = c(parsToMonitor2), n.iter = n.iterations, thin = n.thin)
samples.rjags3 = coda.samples(jm3, variable.names = c(parsToMonitor3), n.iter = n.iterations, thin = n.thin)
samples.rjags4 = coda.samples(jm4, variable.names = c(parsToMonitor4), n.iter = n.iterations, thin = n.thin)

fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/appendix-1/")
dir.create(file.path(fileDirectory), showWarnings = FALSE)

saveRDS(samples.rjags1,file=paste0(fileDirectory,"samples.rjags1.rds"))
saveRDS(samples.rjags2,file=paste0(fileDirectory,"samples.rjags2.rds"))
saveRDS(samples.rjags3,file=paste0(fileDirectory,"samples.rjags3.rds"))
saveRDS(samples.rjags4,file=paste0(fileDirectory,"samples.rjags4.rds"))
saveRDS(data,file=paste0(fileDirectory,"data.rds"))

# figures
data = data.frame(cbind(y,n))

parameterization.MeanConcentration <- samples.rjags1
parameterization.MomentMatched <- samples.rjags2
parameterization.LogitCentered <- samples.rjags3
parameterization.LogitNoncentered <- samples.rjags4

colorVector = c("black","red","blue","purple")

# examine correlation in posterior
shapeit <- function(mu, sigma){
  a <- (mu^2 - mu^3 - mu * sigma^2)/sigma^2
  b <- (mu - 2 * mu^2 + mu^3 - sigma^2 + mu * sigma^2)/sigma^2
  shape_ps <- list(a, b)
  return(shape_ps)
}

i = sample(1:dim(data)[1],1)
samples<-sample(1:dim(vfreq<-MCMCchains(parameterization.MeanConcentration,params="vfreq")
)[1],1000)

# parameterization 1
params.MeanConcentration <- list()
vfreq<-MCMCchains(parameterization.MeanConcentration,params="vfreq")
vcount<-MCMCchains(parameterization.MeanConcentration,params="vcount")
theta_i<-MCMCchains(parameterization.MeanConcentration,params="theta")[,i]

alpha <- vfreq*vcount
beta <- (1-vfreq)*vcount

params.MeanConcentration <- list(vfreq = vfreq,vcount = vcount,
                              theta_i = theta_i, alpha = alpha,
                              beta = beta)

# parameterization 2
params.MomentMatched <- list()

mu<-MCMCchains(parameterization.MomentMatched,params="mu")
sigma<-MCMCchains(parameterization.MomentMatched,params="sigma")
theta_i<-MCMCchains(parameterization.MomentMatched,params="theta")[,i]

betaParms<-shapeit(mu,sigma)

params.MomentMatched <- list(mu = mu,sigma = sigma,
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


# Make a figure for the hyperparameter space.

pdf(paste0(dirFigures,"hyperparameterCorrelations.pdf"), width=8, height=6)

par(mfrow=c(2,2))

a = params.MeanConcentration$alpha[samples]
b = params.MeanConcentration$beta[samples]

plot(log(a/b),log(a+b),
     pch=16,cex=.5,
     col = colorVector[1])

a = params.MomentMatched$alpha[samples]
b = params.MomentMatched$beta[samples]

plot(log(a/b),log(a+b),
     pch=16,cex=.5,
     col = colorVector[2])

plot(params.LogitCentered$mu[samples],params.LogitCentered$sigma[samples],
     pch=16,cex=.5,
     col = colorVector[3],
     xlab="mu",ylab="log(sigma)")

plot(params.LogitNoncentered$mu[samples],params.LogitNoncentered$sigma[samples],
     pch=16,cex=.5,
     col = colorVector[4],
     xlab="mu",ylab="log(sigma)")
dev.off()

# Plot the hyperparameters against thetas.

pdf(paste0(dirFigures,"funnelCheck.pdf"), width=8, height=6)

par(mfrow=c(2,2))
a = params.MeanConcentration$alpha[samples]
b = params.MeanConcentration$beta[samples]
theta_i = params.MeanConcentration$theta_i[samples]

plot(boot::logit(theta_i),log(a+b),
     pch=16,cex=.5,
     col = colorVector[1],
     xlab="log-odds probability of success")

a = params.MomentMatched$alpha[samples]
b = params.MomentMatched$beta[samples]
theta_i = params.MomentMatched$theta_i[samples]

plot(boot::logit(theta_i),log(a+b),
     pch=16,cex=.5,
     col = colorVector[2])

plot(params.LogitCentered$alpha_i[samples],
     log(params.LogitCentered$sigma[samples]),
     pch=16,cex=.5,
     col = colorVector[3],
     ylab="log(sigma)")

plot(params.LogitNoncentered$alpha.std_i[samples],
     log(params.LogitNoncentered$sigma[samples]),
     pch=16,cex=.5,
     col = colorVector[4],
     ylab="log(sigma)")

dev.off()

# Visualize effective sample size.

pdf(paste0(dirFigures,"nEff.pdf"), width=6, height=6)
par(mfrow=c(1,1)) 
f<-function(x) x[order(x,decreasing=TRUE)] 

nEff.MeanConcentration = f(MCMCsummary(parameterization.MeanConcentration,params=c("theta"))$n.eff)
nEff.MomentMatched = f(MCMCsummary(parameterization.MomentMatched,params=c("theta"))$n.eff)
nEff.LogitCentered = f(MCMCsummary(parameterization.LogitCentered,params=c("theta"))$n.eff)
nEff.LogitNoncentered = f(MCMCsummary(parameterization.LogitNoncentered,params=c("theta"))$n.eff)
nEff = c(nEff.MeanConcentration,nEff.MomentMatched,nEff.LogitCentered,nEff.LogitNoncentered)

plot(nEff.MeanConcentration,
     ylim=c(min(nEff),max(nEff)),type='l',
     col=colorVector[1],
     xlab="theta_i (ordered)",
     ylab="Effective sample size")
lines(nEff.MomentMatched,
      col=colorVector[2])
lines(nEff.LogitCentered,
      col=colorVector[3])
lines(nEff.LogitNoncentered,
      col=colorVector[4])

legend(x = 0,
       y = min(nEff)+1/4*(max(nEff)-min(nEff)),
       legend=c("Mean & sample size","Mean & variance", "Centered logit", "Noncentered logit"),
       lty = c(1,1,1,1),
       col = colorVector,
       cex=.75)
dev.off()

# Plot shrinkage.

data <- data.frame(cbind(y=y,n=n))

data$freq <- data$y/data$n

summary.MeanConcentration<-apply(MCMCchains(parameterization.MeanConcentration,params="theta"),2,quantile, probs = c(.05,.5,.95))
summary.MomentMatched<-apply(MCMCchains(parameterization.MomentMatched,params="theta"),2,quantile, probs = c(.05,.5,.95))
summary.LogitCentered<-apply(MCMCchains(parameterization.LogitCentered,params="theta"),2,quantile, probs = c(.05,.5,.95))
summary.LogitNoncentered<-apply(MCMCchains(parameterization.LogitNoncentered,params="theta"),2,quantile, probs = c(.05,.5,.95))

summaryList<-list(summary.MeanConcentration,summary.MomentMatched,summary.LogitCentered,summary.LogitNoncentered)

min.test<-list()
max.test<-list()
for(i in 1:4){
  min.test[[i]]<-t(summaryList[[i]])[,1]
  max.test[[i]]<-t(summaryList[[i]])[,3]
}
min.val=min(unlist(min.test))-.1
max.val=max(unlist(max.test))+.1

pdf(paste0(dirFigures,"shrinkage.pdf"), width=8, height=6)

par(mfrow=c(2,2))
for(i in 1:4){
  plot(x = data$freq,
       y = t(summaryList[[i]])[,2],
       xlab= "Tumor frequency in sample",
       ylab= "Median and 95% CI",
       pch=16,cex=.5,
       xlim=c(min.val,max.val),ylim=c(min.val,max.val),
       col=colorVector[i],type='n')
  abline(a=0,b=1)
  mle=sum(y)/sum(n)
  abline(h=mle,lty='dashed')
  
  
  segments(x0=data$freq,
           y0= t(summaryList[[i]])[,1],
           y1= t(summaryList[[i]])[,3],
           col=colorVector[i])
  points(x = data$freq,
         y = t(summaryList[[i]])[,2],
         pch=16,cex=.5,
         col=colorVector[i])
}

dev.off()

# Plot shrinkage vs. sample size. 

# expand shrinkage plot
min.test<-list()
for(i in 1:4){
  min.test[[i]]<-data$freq-t(summaryList[[i]])[,2]
}


pdf(paste0(dirFigures,"shrinkageSamplesize.pdf"), width=8, height=6)
par(mfrow=c(1,2))

plot(data$n,
     data$freq-t(summaryList[[i]])[,2],
     ylim=c(min(unlist(min.test)),max(unlist(min.test))),
     xlab="Number of trials",
     ylab="Amount of shrinkage",type='n')
abline(h=0)
for(i in 1:4){
  points(data$n,
         data$freq-t(summaryList[[i]])[,2],
         col=colorVector[i],pch=1)
}

plot(rep(c(1,2),each=length(data$freq)),
     c(data$freq,t(summaryList[[i]])[,2]),
     xlab="Observed vs. estimated",
     ylab="Probability",
     ylim=c(0,1),
     pch=16,cex=0.5,
     col=colorVector[i],type='n')
for(i in 1:4){
  segments(x0=1,y0=data$freq,
           x1=2,y1=t(summaryList[[i]])[,2],
           lwd=0.5,
           col=colorVector[i])
}
abline(h=sum(y)/sum(n),lty='dashed')

dev.off()

# expand shrinkage plot
# 
# par(mfrow=c(2,2))
# 
# for(i in 1:4){
#   plot(rep(c(1,2),each=length(data$freq)),
#        c(data$freq,t(summaryList[[i]])[,2]),
#        xlab="Number of trials",
#        ylab="Amount of shrinkage",
#        ylim=c(0,1),
#        pch=16,cex=0.5,
#        col=colorVector[i])
#   segments(x0=1,y0=data$freq,
#            x1=2,y1=t(summaryList[[i]])[,2],
#            lwd=0.5,
#            col=colorVector[i])
# }



# Plot posterior density for mean probability.

betaParms=list(params.MeanConcentration$alpha,params.MeanConcentration$beta)
mean.1 = betaParms[[1]]/(betaParms[[1]]+betaParms[[2]])
betaParms<-shapeit(params.MomentMatched$mu,params.MomentMatched$sigma)
mean.2 = betaParms[[1]]/(betaParms[[1]]+betaParms[[2]])
mean.3=boot::inv.logit(MCMCchains(parameterization.LogitCentered,params="mu"))
mean.4=boot::inv.logit(MCMCchains(parameterization.LogitNoncentered,params="mu"))

means<-cbind(mean.1,mean.2,mean.3,mean.4)
max.function<-function(x) max(density(x)$y)
density.max<-max(apply(means,2,max.function))
means<-t(apply(means,2,quantile,probs=c(0.025,.5,.975)))


plot(density(mean.1),type='n',
     xlim=c(0,1),
     ylim=c(0,density.max + 20),
     main="Posterior distributions for mean probability",
     xlab="Mean probability")
lines(density(mean.1))
lines(density(mean.2),col="red")
lines(density(mean.3),col="blue")
lines(density(mean.4),col="purple")

median.vec = c(density.max+seq(4,18,length.out=4))

points(means[,2],y=median.vec,ylim=c(0,1),
       pch=16,col=c('black','red','blue','purple'))
segments(x0=means[,1],y0=median.vec,x1=means[,3],
         col=c('black','red','blue','purple'))
abline(v=sum(y)/sum(n),lty='dashed')

legend(x = .6,
       y = max(median.vec) ,
       legend=c("Mean & sample size","Mean & variance", "Centered logit", "Noncentered logit"),
       lty = c(1,1,1,1),
       col = c("black","red","blue","purple"),
       cex=.75)


# mean.deltas<-cbind(mean.2-mean.1,mean.3-mean.1,mean.4-mean.1,mean.3-mean.2,mean.4-mean.2,mean.4-mean.3)
# max.function<-function(x) max(density(x)$y)
# density.max<-max(apply(mean.deltas,2,max.function))
# mean.deltas.summary<-t(apply(mean.deltas,2,quantile,probs=c(0.025,.5,.975)))
# 
# plot(density(mean.deltas[,1]),type='n',
#      xlim=c(-.1,.1),
#      ylim=c(0,density.max + 30),
#      main="Posterior distributions for \n difference in mean probability",
#      xlab="Mean probability")
# lines(density(mean.deltas[,1]))
# lines(density(mean.deltas[,2]),col="red")
# lines(density(mean.deltas[,3]),col="blue")
# lines(density(mean.deltas[,4]),col="purple")
# lines(density(mean.deltas[,5]),col="purple")
# lines(density(mean.deltas[,6]),col="purple")
# 
# median.vec = c(density.max+seq(4,24,length.out=6))
# 
# points(mean.deltas.summary[,2],y=median.vec,ylim=c(0,1),
#        pch=16,col=c('black','red','blue','purple'))
# segments(x0=mean.deltas.summary[,1],y0=median.vec,x1=mean.deltas.summary[,3],
#          col=c('black','red','blue','purple'))
# abline(v=0)

# Plot marginal probability for population. 

pdf(paste0(dirFigures,"posteriorProbability.pdf"), width=6, height=6)

# compare marginalized probabilities
p.1<-MCMCchains(parameterization.MeanConcentration,params="p")
p.2<-MCMCchains(parameterization.MomentMatched,params="p")
p.3<-MCMCchains(parameterization.LogitCentered,params="p")
p.4<-MCMCchains(parameterization.LogitNoncentered,params="p")

ps<-cbind(p.1,p.2,p.3,p.4)
max.function<-function(x) max(density(x)$y)
density.max<-max(apply(ps,2,max.function))
ps<-t(apply(ps,2,quantile,probs=c(0.025,.5,.975)))

plot(density(p.1),type='n',
     xlim=c(0,1),
     ylim=c(0,density.max+8),
     main="Posterior distributions for population-level probability",
     xlab="Probability")
lines(density(p.1))
lines(density(p.2),col="red")
lines(density(p.3),col="blue")
lines(density(p.4),col="purple")

median.vec = c(density.max+seq(2:5))

points(ps[,2],y=median.vec,ylim=c(0,1),
       pch=16,col=c('black','red','blue','purple'))
segments(x0=ps[,1],y0=median.vec,x1=ps[,3],
         col=c('black','red','blue','purple'))
abline(v=sum(y)/sum(n),lty='dashed')

legend(x = .5,
       y = max(median.vec)+4 ,
       legend=c("Mean & sample size","Mean & variance", "Centered logit", "Noncentered logit"),
       lty = c(1,1,1,1),
       col = c("black","red","blue","purple"),
       cex=.75)
dev.off()