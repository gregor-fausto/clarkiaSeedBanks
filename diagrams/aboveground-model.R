# -------------------------------------------------------------------
# Simulate dataset
# -------------------------------------------------------------------
n.obs = 30
n.years = 5

mu0 = rnorm(1,-1.5,1)
sigma0 = runif(1,.5,1.5)

mu = rnorm(n.years,mean=mu0,sd=sigma0)
sigma = runif(n.years,min=0,max=1.5)

theta=list()
theta_p=list()
n = list()
y = list()
for(i in 1:n.years){
  theta[[i]] = rnorm(n.obs,mean=mu[i],sd=sigma[i])
  theta_p[[i]]=boot::inv.logit(theta[[i]])
  n[[i]]=sample(1000,n.obs,replace=TRUE)
  y[[i]] = rbinom(n.obs, n[[i]], theta_p[[i]])
}

y = unlist(y)
n = unlist(n)
site = rep(as.factor(1),length(n))
year = rep(as.factor(1:n.years),each=n.obs)

df=data.frame(fruitplNumber=y,seedlingNumber=n,site=site,year=year)

# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(rjags)
library(tidybayes)
library(tidyverse)
library(parallel)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------
# -------------------------------------------------------------------

data <- tidybayes::compose_data(df)

detach("package:tidyverse", unload=TRUE)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
# -------------------------------------------------------------------
# -------------------------------------------------------------------

n.adapt = 3000
n.update = 5000
n.iterations = 2000
n.thin = 1

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Describe functions
# -------------------------------------------------------------------
# -------------------------------------------------------------------

initsMu0 <- function(samps = data$n_site){
  rnorm(n = samps, mean = 0, sd = 1)
}

initsSigma0 <- function(samps = data$n_site){
  extraDistr::rhnorm(n = samps, sigma = 1)
}

initsSigma <- function(rows = data$n_site, cols = data$n_year){
  matrix(extraDistr::rhnorm(n = rows*cols, sigma = 1), rows, cols)
}

# # set inits for JAGS
inits <- list()
for(i in 1:3){
  inits[[i]] <- list(initsMu0(), initsSigma0(), initsSigma() )
  
  names(inits[[i]]) = c("mu0","sigma0","sigma")
  
}

# # Call to JAGS
#
# # tuning (n.adapt)
jm = jags.model("~/Dropbox/clarkiaSeedBanks/scriptsModelFitting/jagsScripts/jags-seedlingSurvival.R",
                data = data, inits = inits,n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

parsToMonitor = c("mu0","sigma0","mu","sigma","theta")
# parsToCheck = c(# LOG LIKELIHOODS
#   # "logLik",
#   # POSTERIOR PREDICTIVE
#   "fruitplNumber_sim",
#   "e.mu","mu",
#   # MODEL CHECKING 
#   "chi2.obs","chi2.sim")


# chain (n.iter)
samples.rjags = coda.samples(jm, 
                             variable.names = c(parsToMonitor), 
                             n.iter = n.iterations, thin = n.thin)

library(MCMCvis)
MCMCsummary(samples.rjags,params=c("mu0","sigma0","mu","sigma"))

mu0.post = MCMCchains(samples.rjags,params="mu0")
sigma0.post = MCMCchains(samples.rjags,params="sigma0")
mu.post = MCMCchains(samples.rjags,params="mu")
sigma.post = MCMCchains(samples.rjags,params="sigma")

par(mfrow=c(1,2))
hist(mu0.post,col='gray90',border=0,breaks=50)
abline(v=mu0,col='orange',lwd=2)

hist(sigma0.post,col='gray90',border=0,breaks=50)
abline(v=sigma0,col='orange',lwd=2)

par(mfcol=c(1,2))
mu.post.sum=apply(mu.post,2,quantile,c(.025,.25,.5,.75,.975))

plot(mu,mu.post.sum[3,],ylim=c(min(mu.post.sum),max(mu.post.sum)),
     pch=16,xlim=c(min(mu.post.sum),max(mu.post.sum)))
segments(x0=mu,y0=mu.post.sum[2,],y1=mu.post.sum[4,],lwd=3)
segments(x0=mu,y0=mu.post.sum[1,],y1=mu.post.sum[5,])
abline(a=0,b=1)

sigma.post.sum=apply(sigma.post,2,quantile,c(.025,.25,.5,.75,.975))

plot(sigma,sigma.post.sum[3,],ylim=c(min(sigma.post.sum),max(sigma.post.sum)),
     pch=16,xlim=c(min(sigma.post.sum),max(sigma.post.sum)))
segments(x0=sigma,y0=sigma.post.sum[2,],y1=sigma.post.sum[4,],lwd=3)
segments(x0=sigma,y0=sigma.post.sum[1,],y1=sigma.post.sum[5,])
abline(a=0,b=1)

## COMPARE MODEL ESTIMATES AND BLUPS
library(lme4)
m1=glmer(cbind(fruitplNumber,seedlingNumber)~0+(1|year),family='binomial',data=df)
# without multiple sites, can't nest random effects but the notation would be (1|site/year) or (1|site) + (1|site:year)
m2=glmer(cbind(fruitplNumber,seedlingNumber)~1+(1|year),family='binomial',data=df)

summary(m2)
MCMCsummary(samples.rjags,params=c("mu0","sigma0"))

hist(mu0.post,col='gray90',border=0,breaks=50)
abline(v=c(mu0,fixef(m2)),col='orange',lwd=2,lty=c('solid','dotted'))

mu.hat=cbind(mu.true = mu,mu.glmer = fixef(m2)+ranef(m2)[[1]],mu.hb = apply(mu.post,2,median))
colnames(mu.hat)= c("mu.true","mu.glmer","mu.hb")
cor(mu.hat)
pairs(mu.hat)
counts=df %>% 
  dplyr::mutate(p = fruitplNumber/seedlingNumber) %>%
  dplyr::group_by(year) %>% 
  dplyr::summarise(p.mean=mean(p),p.var=var(p))
mu.hat=cbind(mu.hat,counts[,2:3])
par(mfcol=c(2,3))

plot(mu.hat$p.mean,mu.hat$mu.true-mu.hat$mu.glmer)
abline(h=0)
plot(mu.hat$p.var,mu.hat$mu.true-mu.hat$mu.glmer)
abline(h=0)

plot(mu.hat$p.mean,mu.hat$mu.true-mu.hat$mu.hb)
abline(h=0)
plot(mu.hat$p.var,mu.hat$mu.true-mu.hat$mu.hb)
abline(h=0)

plot(mu.hat$p.mean,mu.hat$mu.hb-mu.hat$mu.glmer)
abline(h=0)
plot(mu.hat$p.var,mu.hat$mu.hb-mu.hat$mu.glmer)
abline(h=0)


par(mfrow=c(1,1))
mu.post.sum=apply(mu.post,2,quantile,c(.025,.25,.5,.75,.975))
mu.glmer = as.matrix(mu.hat$mu.glmer)
plot(mu.glmer[,1],mu.post.sum[3,],ylim=c(min(mu.post.sum),max(mu.post.sum)),
     pch=16,xlim=c(min(mu.post.sum),max(mu.post.sum)))
segments(x0=mu.glmer[,1],y0=mu.post.sum[2,],y1=mu.post.sum[4,],lwd=3)
segments(x0=mu.glmer[,1],y0=mu.post.sum[1,],y1=mu.post.sum[5,])
abline(a=0,b=1)


plot(mu,mu.glmer[,1],pch=16,cex=)
abline(a=0,b=1)
abline(h=fixef(m2))

par(mfrow=c(2,5))
for(i in 1:10){
  hist(mu.post[,i],col='gray90',border=0,main='')
abline(v=c(mu[i],mu.glmer[i,1]),col='orange',lwd=2,lty=c('solid','dotted'))
}

# MODEL FIGURE
f=function(x){
  tmp=density(x,from=min(x),to=max(x),adjust=2.5)
  df=cbind(tmp$x,tmp$y)
  df=rbind(c(tmp$x[1],0),df)
  df=rbind(df,c(tmp$x[length(tmp$x)],0))
  return(df)
}


## PLOT VERSION 1
par(mfrow=c(1,1))
plot(x=NA,NA,
     type='n',
     ylim=c(0,12),xlim=c(0,1),
     #axes=FALSE,
     frame=FALSE,xaxt='n',yaxt='n',
     xlab="Probability of seedling survival to fruiting",
     ylab="Year",cex=3,cex.axis=3)
year=1:5

sigma = boot::inv.logit(mu.post)

for(i in 1:5){
  df.tmp=sigma[,i]
  full.post=rnorm(length(mu.post[,i]),mean=mu.post[,i],sd=sigma.post[,i])
  full.post = boot::inv.logit(full.post)
  dat.tmp = df[df$year==i,] %>%
    dplyr::mutate(p = fruitplNumber/seedlingNumber) %>%
    dplyr::filter(!is.na(p))

  upper.limit=max(f(full.post)[,2])*.8
  polygon(y=.5+2*year[i]+f(full.post)[,2]/upper.limit,x=f(full.post)[,1],col='#1b9e77',border='#1b9e77')

  upper.limit=max(f(df.tmp)[,2])*1
  polygon(y=.5+2*year[i]+f(df.tmp)[,2]/upper.limit,x=f(df.tmp)[,1],col='#d95f02',border='#d95f02')
 
  # f.boxplot(year[i],df.tmp)
  n.obs = length(dat.tmp$seedlingNumber)
  size = dat.tmp$seedlingNumber/max(df$seedlingNumber,na.rm=TRUE)
  points(y=.5+2*rep(i-0.1,n.obs)+rnorm(n.obs,0,.025),x=dat.tmp$p,
         pch = 19, cex = .75,col=rgb(0,0,0,.5))
}

axis(1, seq(0,1,by=.2),
     labels = seq(0,1,by=.2), las = 1,
     col = 'black', col.ticks = 1, cex.axis = 1.25)
axis(2, 2*(c(2:6)-.5),
     labels = c(1:5), las = 1,
     col = 'black', col.ticks = 1, cex.axis = 1.25)

text(.4,11.91,
     labels = "Population and year level (parameters)" ,cex = .8,pos=4)
text(.75,.91,
     labels = "Population level (hyperparameters)" ,cex = .8)

df.tmp=boot::inv.logit(mu0.post)
full.post=rnorm(length(mu0.post),mean=mu0.post,sd=sigma0.post)
full.post = boot::inv.logit(full.post)

upper.limit=max(f(full.post)[,2])*.8
polygon(y=f(full.post)[,2]/upper.limit,x=f(full.post)[,1],col='#7570b3',border='#7570b3')

upper.limit=max(f(df.tmp)[,2])*1
polygon(y=f(df.tmp)[,2]/upper.limit,x=f(df.tmp)[,1],col='gray90',border='gray90')

abline(h=1.5,lty='dotted')
box()


pdf("~/Dropbox/clarkiaSeedBanks/products/manuscript/figures-overview/marginal-posterior.pdf",width=6,height=6)
## PLOT VERSION 2
par(mfrow=c(1,1),mar=c(4,.5,.5,.5),oma=c(1,3,1,1))
plot(x=NA,NA,
     type='n',
     ylim=c(0,12),xlim=c(0,1),
     frame=FALSE,xaxt='n',yaxt='n',
     xlab="Probability of seedling survival to fruiting",
     ylab="Year",cex=3,cex.axis=3,
     cex.lab = 1.5,line=2.15)
year=1:5

sigma = boot::inv.logit(mu.post)

for(i in 1:5){
  df.tmp=sigma[,i]
  full.post=rnorm(length(mu.post[,i]),mean=mu.post[,i],sd=sigma.post[,i])
  full.post = boot::inv.logit(full.post)
  dat.tmp = df[df$year==i,] %>%
    dplyr::mutate(p = fruitplNumber/seedlingNumber) %>%
    dplyr::filter(!is.na(p))
  
  upper.limit=max(f(full.post)[,2])*.8
  polygon(y=.5+2*year[i]+f(full.post)[,2]/upper.limit,x=f(full.post)[,1],col='gray80',border='gray80')

  n.obs = length(dat.tmp$seedlingNumber)
  size = dat.tmp$seedlingNumber/max(df$seedlingNumber,na.rm=TRUE)
  points(y=.5+2*rep(i-0.1,n.obs)+rnorm(n.obs,0,.025),x=dat.tmp$p,
         pch = 19, cex = .75,col=rgb(0,0,0,.5))
}

axis(1, seq(0,1,by=.2),
     labels = seq(0,1,by=.2), las = 1,
     col = 'black', col.ticks = 1, cex.axis = 1.25)
axis(2, 2*(c(2:6)-.5),
     labels = c(1:5), las = 1,
     col = 'black', col.ticks = 1, cex.axis = 1.25)
mtext("Year",side=2,adj=.6,line=2,cex=2,cex.lab=1.7)

text(.5,11.9,
     labels = "Population and year level" ,cex = 1.2, pos = 4)
text(.7,.9,
     labels = "Population level" ,cex = 1.2, pos = 4)

df.tmp=boot::inv.logit(mu0.post)
full.post=rnorm(length(mu0.post),mean=mu0.post,sd=sigma0.post)
full.post = boot::inv.logit(full.post)

upper.limit=max(f(full.post)[,2])*.8
polygon(y=f(full.post)[,2]/upper.limit,x=f(full.post)[,1],col='gray80',border='gray80')

abline(h=1.5,lty='dotted')
box()
mtext("Marginalized probabilitities",adj=0,cex=1.75)

dev.off()

pdf("~/Dropbox/clarkiaSeedBanks/products/manuscript/figures-overview/parameter.pdf",width=6,height=6)
## PLOT VERSION 3
par(mfrow=c(1,2),mar=c(4,.5,.5,.5),oma=c(1,3,1,1))
plot(x=NA,NA,
     type='n',
     ylim=c(0,12),xlim=c(-4,2),
     frame=FALSE,xaxt='n',yaxt='n',
     xlab="Means",
     ylab="Year",cex=3,cex.axis=3,
     cex.lab = 1.7,line=2.2)
year=1:5

sigma = (mu.post)

for(i in 1:5){
  
  df.tmp=sigma[,i]
  dat.tmp = df[df$year==i,] %>%
    dplyr::mutate(p = fruitplNumber/seedlingNumber) %>%
    dplyr::filter(!is.na(p))
  
  upper.limit=max(f(df.tmp)[,2])*1
  polygon(y=.5+2*year[i]+f(df.tmp)[,2]/upper.limit,x=f(df.tmp)[,1],col='#d95f02',border='#d95f02')
  
  n.obs = length(dat.tmp$seedlingNumber)
  size = dat.tmp$seedlingNumber/max(df$seedlingNumber,na.rm=TRUE)

}

axis(1, seq(-6,2,by=2),
     labels = seq(-6,2,by=2), las = 1,
     col = 'black', col.ticks = 1, cex.axis = 1)
axis(2, 2*(c(2:6)-.5),
     labels = c(1:5), las = 1,
     col = 'black', col.ticks = 1, cex.axis = 1)


df.tmp=(mu0.post)

upper.limit=max(f(df.tmp)[,2])*1
polygon(y=f(df.tmp)[,2]/upper.limit,x=f(df.tmp)[,1],col='#7570b3',border='#7570b3')

abline(h=1.5,lty='dotted')
box()
mtext("Year",side=2,adj=.6,line=2,cex=2,cex.lab=1.7)
mtext("Model parameters",adj=0,cex=1.75)

plot(x=NA,NA,
     type='n',
     ylim=c(0,12),xlim=c(0,2),
     #axes=FALSE,
     frame=FALSE,xaxt='n',yaxt='n',
     xlab="Standard deviations",
     ylab="Year",cex=3,cex.axis=3,
     cex.lab = 1.7,line=2.2)

sigma = (sigma.post)

for(i in 1:5){
  df.tmp=sigma[,i]
  dat.tmp = df[df$year==i,] %>%
    dplyr::mutate(p = fruitplNumber/seedlingNumber) %>%
    dplyr::filter(!is.na(p))
  upper.limit=max(f(df.tmp)[,2])*1
  polygon(y=.5+2*year[i]+f(df.tmp)[,2]/upper.limit,x=f(df.tmp)[,1],col='#d95f02',border='#d95f02')
  
  n.obs = length(dat.tmp$seedlingNumber)
  size = dat.tmp$seedlingNumber/max(df$seedlingNumber,na.rm=TRUE)
}

axis(1, seq(0,4,by=1),
     labels = seq(0,4,by=1), las = 1,
     col = 'black', col.ticks = 1, cex.axis = 1)


df.tmp=boot::inv.logit(mu0.post)

upper.limit=max(f(df.tmp)[,2])*1
polygon(y=f(df.tmp)[,2]/upper.limit,x=f(df.tmp)[,1],col='#7570b3',border='#7570b3')

abline(h=1.5,lty='dotted')
box()

dev.off()