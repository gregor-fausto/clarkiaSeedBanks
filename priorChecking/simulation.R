n = 100
lambda = 1.25
t = seq(0,1,by=.01)
# conditional probability of survival (all intact seeds have x% germination probability)
g = c(.25,.25,.25)

par(mfrow=c(1,2))

# f.exp = function(t,lam=lambda){
#   return(exp(-lam*t))
# }

beta=1/lambda
alpha=.75

f.exp = function(t,lam=lambda){
  return(exp(-(t/beta)^alpha))
}

# plot exponential survival function
plot(t,f.exp(t=t),xlim=c(0,1),ylim=c(0,1),type='l')


# calculate age-specific survival probability
f=function(x){
  f.exp(t.sample[x+1])/f.exp(t.sample[x])
}
t.sample = c(0,3,12,15,24,27,36)/36
# instantaneous survival is identical; differencs are due to length in time period
s.age = f(c(1,2,3,4,5,6,NA))
# calculate survivorship schedule
l.x=c()
for(i in 1:length(s.age)){
  l.x[i]=prod(s.age[1:i])
}
# plot survivorship points on to curve
points(t.sample[2:length(t.sample)],l.x[1:(length(l.x)-1)],xlim=c(0,1),ylim=c(0,1),
       pch = c(1,1,1,1,1,1))


## Add germination to the plot
t.sample = c(0,3,3,12,15,15,24,27,27,36)/36

# calculate age-specific survival probability
# full set
Oct_0 = 1
# P(S)
Janpre_1 = f.exp(t=t.sample[2])
# unconditional probability of germination: P(G|S)P(S)=P(G)
Jangerm_1 = g[1]
Janpost_1 = (1-g[1])
Oct_1 = f.exp(t=t.sample[4])/f.exp(t=t.sample[3])

Janpre_2 = f.exp(t=t.sample[5])/f.exp(t=t.sample[4])
Jangerm_2 = g[2]
Janpost_2 = (1-g[2])
Oct_2 = f.exp(t=t.sample[7])/f.exp(t=t.sample[6])

Janpre_3 = f.exp(t=t.sample[8])/f.exp(t=t.sample[7])
Jangerm_3 = g[3]
Janpost_3 = (1-g[3])
Oct_3 = f.exp(t=t.sample[10])/f.exp(t=t.sample[9])


s.age = c(Oct_0,Janpre_1,Jangerm_1,Janpost_1,Oct_1)

# calculate survivorship schedule
l.x=c()
l.x[1]=1
l.x[2]=Janpre_1
l.x[3]=Janpost_1*Janpre_1
l.x[4]=Oct_1*Janpost_1*Janpre_1
l.x[5]=Janpre_2*Oct_1*Janpost_1*Janpre_1
l.x[6]=Janpost_2*Janpre_2*Oct_1*Janpost_1*Janpre_1
l.x[7]=Oct_2*Janpost_2*Janpre_2*Oct_1*Janpost_1*Janpre_1
l.x[8]=Janpre_3*Oct_2*Janpost_2*Janpre_2*Oct_1*Janpost_1*Janpre_1
l.x[9]=Janpost_3*Janpre_3*Oct_2*Janpost_2*Janpre_2*Oct_1*Janpost_1*Janpre_1
l.x[10]=Oct_3*Janpost_3*Janpre_3*Oct_2*Janpost_2*Janpre_2*Oct_1*Janpost_1*Janpre_1
l.x[11]=NA

# points(t.sample[1:5],s.age,xlim=c(0,1),ylim=c(0,1),
#        pch = c(1,1,16,1,1),col='red',cex=.5)

points(t.sample[1:10],l.x[1:(length(l.x)-1)],xlim=c(0,1),ylim=c(0,1),
       col='red',cex=.5,
       pch=c(16,16,1,16,16,1,16,16,1,16))


# build table for model fitting formulation
t.sample = c(0,3,12,15,24,27,36)/36
counts = c(
  f.exp(t=t.sample[1]),
  f.exp(t=t.sample[2]),
  #  g[1]*exp(-lambda*t.sample[3]),
  (1-g[1])*f.exp(t=t.sample[3]),
  (1-g[1])*f.exp(t=t.sample[4]),
  # (1-g[1])*g[2]*exp(-lambda*t.sample[6]),
  (1-g[1])*(1-g[2])*f.exp(t=t.sample[5]),
  (1-g[1])*(1-g[2])*f.exp(t=t.sample[6]),
  # (1-g[1])*(1-g[2])*g[3]*exp(-lambda*t.sample[9]),
  (1-g[1])*(1-g[2])*(1-g[3])*f.exp(t=t.sample[7])
)

points(t.sample,counts,xlim=c(0,1),ylim=c(0,1),
       
       col='blue')


plot(c(1,2,3),g,ylim=c(0,1),pch=16)
#unconditional germination probabilities in red
g.uncond = c(
  Jangerm_1*Janpre_1,
  Jangerm_2*Janpre_2*Oct_1*Janpost_1*Janpre_1,
  Jangerm_3*Janpre_3*Oct_2*Janpost_2*Janpre_2*Oct_1*Janpost_1*Janpre_1)

points(c(1,2,3),g.uncond,xlim=c(0,1),ylim=c(0,1),
       col='red',cex=.5)


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Simulate experiment
# -------------------------------------------------------------------
# -------------------------------------------------------------------

nBags = 100

seedStart = rep(nBags,times=10)

t.sample = c(3,12,15,24,27,36)/36

l.x.trials=l.x[c(2,4,5,7,8,10)]

bags<-replicate(nBags,rbinom(n=6,size=100,l.x.trials))
df.list<-list()
y.list<-list()
for(i in 1:nBags){
 y = rbinom(n=6,size=100,l.x.trials)
 y.list[[i]]<-y
 df.list[[i]]<-data.frame(siteSurvival=as.character("S1"),
                                yearSurvival=as.character("2006"),
                                ageBags=c(1,1,2,2,3,3),
                                y=y,
                                seedStart=rep(100,length(y)),
                                months=t.sample,
                                gIndexSurvival=c(1,1,2,2,3,3),
                                compIndex=c(1,2,3,4,5,6))
}
survivalData = do.call(rbind,df.list)
yData = do.call(cbind,y.list)

par(mfrow=c(1,2))
plot(t.sample,bags[,1],type='n',
     xlim=c(0,1),ylim=c(0,100))
for(i in 1:nBags){
  points(t.sample,bags[,i])
}

c(l.x[2],l.x[5],l.x[8])

germinants=bags[c(1,3,5),]

for(i in 1:nBags){
  germinants[,i]=rbinom(n=3,bags[c(1,3,5),i],g)
}

df.list<-list()
for(i in 1:nBags){
  y.germinants = rbinom(n=3,yData[c(1,3,5),i],g)
  n.total = yData[c(1,3,5),i]
  df.list[[i]]<-data.frame(siteGermination=as.character("S1"),
                                 yearGermination=as.character("2006"),
                                 ageBags=c(1,2,3),
                                 seedlingJan=y.germinants,
                                 totalJan=n.total,
                                 gIndex=c(1,2,3))
}
germinationData = do.call(rbind,df.list)


plot(c(1,2,3),germinants[,1],type='n',
     ylim=c(0,100))
for(i in 1:nBags){
  points(c(1,2,3),germinants[,i])
}

#discrete.histogram(y, xlab = "y", main = "Simulated data")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Fit data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
library(tidybayes)

data <- tidybayes::compose_data(survivalData,germinationData)

data$n1 = dim(survivalData)[1]
data$n2 = dim(germinationData)[1]

refDf<-data.frame(year=data$yearGermination,gIndex=data$gIndex) %>%
  unique %>%
  dplyr::mutate(germinationReference=1:3)

data$yearRefGerm = refDf$year
data$indexRefGerm = refDf$gIndex
data$refGerm = refDf$germinationReference

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
library(magrittr)
library(tidybayes)
library(parallel)
library(stringr)
library(loo)

set.seed(10)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.chain = 3
n.adapt = 3000
n.update = 5000
n.iterations = 10000
n.thin = 1

dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/priorChecking/jagsScripts/")

# set inits for JAGS
inits <- list()

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Negative exponential, age-dependent germination
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# 
# # exponential
# jmNeAg = jags.model(paste0(dir,"binomialLikelihood-logLink-expDecay-simulation.R"), data = data, n.adapt = n.adapt,
#                     n.chains=2)
# 
# # burn-in (n.update)
# update(jmNeAg, n.iter = n.update)
# 
# samples.rjags.jmNeAg = coda.samples(jmNeAg,
#                                     variable.names = 
#                                       # GERMINATION
#                                       c("mu0_g","tau0_g","mu_g","tau_g",
#                                         "sigma0_g","sigma_g","g",
#                                         "mu0_pred","tau0_pred","mu_pred","tau_pred",
#                                         "sigma0_pred","sigma_pred","g_pred",
#                                         #PRIOR PREDICTIVE
#                                         "y_pred","y_prior.pred",
#                                         # SURVIVAL
#                                         "mu0_s","mu_s","inv.b","theta_c",
#                                         "sigma0_s","sigma_s",
#                                         "beta","beta_mod","alpha_s","theta_s",
#                                         "mu","a",
#                                         # POSTERIOR PREDICTIVE
#                                         "y_sim","seedlingJan_sim"),
#                                     n.iter = n.iterations, thin = n.thin)


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Weibull residence time, age-dependent germination
# -------------------------------------------------------------------
# -------------------------------------------------------------------

jmWbAg = jags.model(paste0(dir,"binomialLikelihood-logLink-wrDecay-simulation.R"), data = data, n.adapt = n.adapt,
                    n.chains=3)

# burn-in (n.update)
update(jmWbAg, n.iter = n.update)

samples.rjags.jmWbAg = coda.samples(jmWbAg,
                             variable.names =
                               # GERMINATION
                               c("mu0_g","tau0_g","mu_g","tau_g",
                                 "sigma0_g","sigma_g","g",
                                 "mu0_pred","tau0_pred","mu_pred","tau_pred",
                                 "sigma0_pred","sigma_pred","g_pred",
                                 #PRIOR PREDICTIVE
                                 "y_pred","y_prior.pred",
                                 # SURVIVAL
                                 "mu0_s","mu_s","inv.b","theta_c",
                                 "sigma0_s","sigma_s",
                                 "beta","beta_mod","alpha_s","theta_s",
                                 "mu","a",
                                 # POSTERIOR PREDICTIVE
                                 "y_sim","seedlingJan_sim",
                                 "p.chi2","chi2.obs","chi2.sim",
                                 "chi2.yobs","chi2.ysim",
                                 "sd.data","sd.sim","p.sd",
                                 "p.mean","p.cv"),
                             n.iter = n.iterations, thin = n.thin)


MCMCsummary(samples.rjags.jmWbAg, 
            params = c("mu0_g", "sigma0_g", "sigma_g",
                       "mu0_s","sigma0_s", "sigma_s",
                       "a"))

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Graphical Posterior Predictive Checks
# -------------------------------------------------------------------
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# Germination
# -------------------------------------------------------------------
par(mfrow=c(1,3))

hist(data$seedlingJan[data$gIndex==1], breaks = 10, 
     freq = FALSE, main = "Simulated and real data for germination", 
     xlab = expression(paste("germinant count")), cex.lab = 1.2) 
g1=MCMCchains(samples.rjags.jmWbAg, params = "seedlingJan_sim")[,data$gIndex==1]
lines(density(g1,adjust=5), col = "red")

hist(data$seedlingJan[data$gIndex==2], breaks = 10, 
     freq = FALSE, main = "Simulated and real data for germination", 
     xlab = expression(paste("germinant count")), cex.lab = 1.2) 
g2=MCMCchains(samples.rjags.jmWbAg, params = "seedlingJan_sim")[,data$gIndex==2]
lines(density(g2,adjust=5), col = "red")

hist(data$seedlingJan[data$gIndex==3], breaks = 10, 
     freq = FALSE, main = "Simulated and real data for germination", 
     xlab = expression(paste("germinant count")), cex.lab = 1.2) 
g3=MCMCchains(samples.rjags.jmWbAg, params = "seedlingJan_sim")[,data$gIndex==3]
lines(density(g3,adjust=5), col = "red")

# -------------------------------------------------------------------
# Intact seed counts
# -------------------------------------------------------------------
par(mfrow=c(2,3))

for(i in 1:6){
hist(data$y[data$compIndex==i], breaks = 10, 
     freq = FALSE, main = "Simulated and real data for germination", 
     xlab = expression(paste("germinant count")), cex.lab = 1.2) 
y_sim=MCMCchains(samples.rjags.jmWbAg, params = "y_sim")[,data$compIndex==i]
lines(density(y_sim,adjust=5), col = "red")
}

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Posterior Predictive Checks
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Binned residual plots
# see check 6: https://moodle2.units.it/pluginfile.php/290133/mod_resource/content/1/model_check_script.R
# https://discourse.mc-stan.org/t/posterior-prediction-from-logit-regression/12217
# -------------------------------------------------------------------
# 
# sims <- MCMCchains(samples.rjags.jmWbAg, params = "y_sim")
# sims.subset=sims[sample(dim(sims)[1],10000),]
# 
# sims <- MCMCchains(samples.rjags.jmWbAg, params = "seedlingJan_sim")
# 
# for(i in 1:6){
#   
#   par(mfrow=c(2,3))
#   bayesplot::ppc_error_binned(data$seedlingJan[data$gIndex==3], 
#                               sims[1:6,data$gIndex==3])
#   
# }

# -------------------------------------------------------------------
# Bayesian p-values: omnibus Chi-squared for germination counts
# entire dataset, age-specific
# -------------------------------------------------------------------
chi2.obs=MCMCchains(samples.rjags.jmWbAg,params=c("chi2.obs"))
chi2.sim=MCMCchains(samples.rjags.jmWbAg,params=c("chi2.sim"))
# calculations are rowwise
fit.obs=apply(chi2.obs,1,sum)
fit.sim=apply(chi2.sim,1,sum)
p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
mean(p.chi2.calc)

#MCMCsummary(samples.rjags.jmWbAg, params = c("p.chi2"))
# hist(MCMCchains(samples.rjags.jmWbAg, params = c("p.chi2")))
# hist(p.chi2.calc)

p.age = matrix(NA,nrow=dim(chi2.obs)[1],ncol=3)
for(i in 1:3){
  chi2.obs=MCMCchains(samples.rjags.jmWbAg,params=c("chi2.obs"))[,data$gIndex==i]
  chi2.sim=MCMCchains(samples.rjags.jmWbAg,params=c("chi2.sim"))[,data$gIndex==i]
  fit.obs=apply(chi2.obs,1,sum)
  fit.sim=apply(chi2.sim,1,sum)
  p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
  p.age[,i] = p.chi2.calc
}
apply(p.age,2,mean)

par(mfrow=c(1,1))
time.sample = 1:3
plot(time.sample,apply(p.age,2,mean),
     ylim=c(0,1),pch=16,
     xlab="Months",ylab="p-Value",
     main="Germinant counts")
abline(h=c(.1,.9),lty='dotted')

# note that for particular parameter combinations the counts
# are very low in year 3, making it challenging to obtain reasonable
# estimates of conditional germination rates
# see this by reducing lambda and comparing
# but I think this might trade off with ability to estimate survival

# -------------------------------------------------------------------
# Bayesian p-values: omnibus Chi-squared for seed counts
# entire dataset, age-specific
# -------------------------------------------------------------------
chi2.yobs=MCMCchains(samples.rjags.jmWbAg,params=c("chi2.yobs"))
chi2.ysim=MCMCchains(samples.rjags.jmWbAg,params=c("chi2.ysim"))
# calculations are rowwise
fit.obs=apply(chi2.yobs,1,sum)
fit.sim=apply(chi2.ysim,1,sum)
p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
mean(p.chi2.calc)

# hist(MCMCchains(samples.rjags.jmWbAg, params = c("p.chi2")))
# hist(p.chi2.calc)

p.index = matrix(NA,nrow=dim(chi2.obs)[1],ncol=6)
for(i in 1:6){
  chi2.yobs.index=chi2.yobs[,data$compIndex==i]
  chi2.ysim.index=chi2.ysim[,data$compIndex==i]
  fit.obs=apply(chi2.yobs.index,1,sum)
  fit.sim=apply(chi2.ysim.index,1,sum)
  p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
  p.index[,i] = p.chi2.calc
}
apply(p.index,2,mean)

par(mfrow=c(1,1))
time.sample = c(3,12,15,24,27,36)
plot(time.sample,apply(p.index,2,mean),
     ylim=c(0,1),pch=16,
     xlab="Months",ylab="p-Value",
     main="Seed counts")
abline(h=c(.1,.9),lty='dotted')

# -------------------------------------------------------------------
# Bayesian p-values: additional
# entire dataset, age-specific
# -------------------------------------------------------------------
MCMCsummary(samples.rjags.jmWbAg,params=c("p.sd","p.mean","p.cv"))

sd.data=MCMCchains(samples.rjags.jmWbAg,params=c("sd.data"))
sd(data$seedlingJan)
#sd.sim=MCMCchains(samples.rjags.jmWbAg,params=c("sd.sim"))
# plot(apply(MCMCchains(samples.rjags.jmWbAg,params=c("seedlingJan_sim")),1,sd),
#      sd.sim);abline(a=0,b=1)
sd.sim=apply(MCMCchains(samples.rjags.jmWbAg,params=c("seedlingJan_sim")),1,sd)
p.sd = ifelse(sd.sim-sd.data>=0,1,0)
mean(p.sd)

p.sd.age = matrix(NA,nrow=dim(chi2.obs)[1],ncol=3)
for(i in 1:3){
  sd.data=sd(data$seedlingJan[data$gIndex==i])
  y_samples=MCMCchains(samples.rjags.jmWbAg,params=c("seedlingJan_sim"))[,data$gIndex==i]
  sd.sim=apply(y_samples,1,sd)
  p.sd = ifelse(sd.sim-sd.data>=0,1,0)
  p.sd.age[,i] = p.sd
}
apply(p.sd.age,2,mean)

par(mfrow=c(1,1))
time.sample = 1:3
plot(time.sample,apply(p.sd.age,2,mean),
     ylim=c(0,1),pch=16,
     xlab="Months",ylab="p-Value",
     main="Germinant counts")
abline(h=c(.1,.9),lty='dotted')


# -------------------------------------------------------------------
# Bayesian p-values: additional mean
# entire dataset, age-specific
# -------------------------------------------------------------------
mean.data=mean(data$seedlingJan)
#mean.sim=MCMCchains(samples.rjags.jmWbAg,params=c("mean.sim"))
# plot(apply(MCMCchains(samples.rjags.jmWbAg,params=c("seedlingJan_sim")),1,mean),
#      mean.sim);abline(a=0,b=1)
mean.sim=apply(MCMCchains(samples.rjags.jmWbAg,params=c("seedlingJan_sim")),1,mean)
p.mean = ifelse(mean.sim-mean.data>=0,1,0)
mean(p.mean)

p.mean.age = matrix(NA,nrow=dim(chi2.obs)[1],ncol=3)
for(i in 1:3){
  mean.data=mean(data$seedlingJan[data$gIndex==i])
  y_samples=MCMCchains(samples.rjags.jmWbAg,params=c("seedlingJan_sim"))[,data$gIndex==i]
  mean.sim=apply(y_samples,1,mean)
  p.mean = ifelse(mean.sim-mean.data>=0,1,0)
  p.mean.age[,i] = p.mean
}
apply(p.mean.age,2,mean)

par(mfrow=c(1,1))
time.sample = 1:3
plot(time.sample,apply(p.mean.age,2,mean),
     ylim=c(0,1),pch=16,
     xlab="Months",ylab="p-Value",
     main="Germinant counts")
abline(h=c(.1,.9),lty='dotted')

# -------------------------------------------------------------------
# Bayesian p-values: additional median
# entire dataset, age-specific
# -------------------------------------------------------------------
median.data=median(data$seedlingJan)
#median.sim=MCMCchains(samples.rjags.jmWbAg,params=c("median.sim"))
# plot(apply(MCMCchains(samples.rjags.jmWbAg,params=c("seedlingJan_sim")),1,median),
#      median.sim);abline(a=0,b=1)
median.sim=apply(MCMCchains(samples.rjags.jmWbAg,params=c("seedlingJan_sim")),1,median)
p.median = ifelse(median.sim-median.data>=0,1,0)
mean(p.median)

p.median.age = matrix(NA,nrow=dim(chi2.obs)[1],ncol=3)
for(i in 1:3){
  median.data=median(data$seedlingJan[data$gIndex==i])
  y_samples=MCMCchains(samples.rjags.jmWbAg,params=c("seedlingJan_sim"))[,data$gIndex==i]
  median.sim=apply(y_samples,1,median)
  p.median = ifelse(median.sim-median.data>=0,1,0)
  p.median.age[,i] = p.median
}
apply(p.median.age,2,mean)

par(mfrow=c(1,1))
time.sample = 1:3
plot(time.sample,apply(p.median.age,2,mean),
     ylim=c(0,1),pch=16,
     xlab="Months",ylab="p-Value",
     main="Germinant counts")
abline(h=c(.1,.9),lty='dotted')

# -------------------------------------------------------------------
# Bayesian p-values: additional CV
# entire dataset, age-specific
# -------------------------------------------------------------------
cv.data=sd(data$seedlingJan)/mean(data$seedlingJan)
#median.sim=MCMCchains(samples.rjags.jmWbAg,params=c("median.sim"))
# plot(apply(MCMCchains(samples.rjags.jmWbAg,params=c("seedlingJan_sim")),1,median),
#      median.sim);abline(a=0,b=1)
cv.sim=sd.sim/mean.sim
p.cv = ifelse(cv.sim-cv.data>=0,1,0)
mean(p.cv)

p.cv.age = matrix(NA,nrow=dim(chi2.obs)[1],ncol=3)
for(i in 1:3){
  cv.data=sd(data$seedlingJan[data$gIndex==i])/mean(data$seedlingJan[data$gIndex==i])
  y_samples=MCMCchains(samples.rjags.jmWbAg,params=c("seedlingJan_sim"))[,data$gIndex==i]
  cv.sim = apply(y_samples,1,sd)/apply(y_samples,1,mean)
  p.cv = ifelse(cv.sim-cv.data>=0,1,0)
  p.cv.age[,i] = p.cv
}
apply(p.cv.age,2,mean)

par(mfrow=c(1,1))
time.sample = 1:3
plot(time.sample,apply(p.cv.age,2,mean),
     ylim=c(0,1),pch=16,
     xlab="Months",ylab="p-Value",
     main="Germinant counts")
abline(h=c(.1,.9),lty='dotted')


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Model comparison
# -------------------------------------------------------------------
# -------------------------------------------------------------------

LLmat.ne <- MCMCchains(samples.rjags.jmNeAg,params="logLik_y")
rel_n_eff.ne <- relative_eff(exp(LLmat.ne), chain_id = rep(1, each = 2*n.iterations))
looNe <- loo(LLmat.ne, r_eff = rel_n_eff.ne, cores = 2)

print(looNe)
plot(looNe)


LLmat.wb <- MCMCchains(samples.rjags.jmWbAg,params="logLik_y")
rel_n_eff.wb <- relative_eff(exp(LLmat.wb), chain_id = rep(1, each = 2*n.iterations))
looWb <- loo(LLmat.wb, r_eff = rel_n_eff.wb, cores = 2)

print(looWb)
plot(looWb)

loo_compare(looNe,looWb)

MCMCsummary(samples.rjags.jmWbAg,params=c("a"))
MCMCtrace(samples.rjags.jmWbAg,params=c("a"))

MCMCsummary(samples.rjags.jmNeAg,params=c("mu0_s"))
MCMCsummary(samples.rjags.jmNeAg,params=c("mu_s"))
MCMCsummary(samples.rjags.jmNeAg,params=c("mu0_g"))
MCMCsummary(samples.rjags.jmNeAg,params=c("mu_g"))

MCMCsummary(samples.rjags.jmWbAg,params=c("mu0_s"))
MCMCsummary(samples.rjags.jmWbAg,params=c("mu_s"))
MCMCsummary(samples.rjags.jmWbAg,params=c("mu0_g"))
MCMCsummary(samples.rjags.jmWbAg,params=c("mu_g"))
