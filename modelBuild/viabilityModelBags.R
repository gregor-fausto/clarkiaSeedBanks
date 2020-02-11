# -------------------------------------------------------------------
# Models for joint estimates of year 1 belowground rates
# Seed survival, germination, and viability
# Models use log-odds parameterization
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)
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
library(tidyverse)

set.seed(10)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize viability trial data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# setwd and read data files
load("/Users/Gregor/Dropbox/clarkiaSeedBanks/library/dataFromWorkflowFile/viabilityRawData.rds")
df <- viabilityRawData

df <- df %>% 
  dplyr::select(-c(germPerc,germNot,viabPerc,viabPerc2,condTest))
df$bag <-as.integer(as.numeric(df$bagNo))

## FOR NOW REMOVE MISSING DATA AND PROBLEMS
df<-subset(df,!is.na(df$germStart))
df<-subset(df,!is.na(df$germCount))
df<-subset(df,!is.na(df$viabStart))
df<-subset(df,!is.na(df$viabStain))

## data check
df %>% dplyr::filter(germStart - germCount - viabStart<0)

# filter out rows with problems
df<-df %>% 
  dplyr::filter(germStart - germCount - viabStart >= 0)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Multiple-site model (3 year)
# Log-odds parameterization for viability, trials summed across bags
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# AGE 1 SEEDS
dat<-df %>% 
  dplyr::filter(age==1)

# assign variable that combines site and bag; unique id for each bag
 dat<-dat %>% 
   tidyr::unite(col='siteBag', c(site,bagNo), sep="", remove=FALSE) %>%
   dplyr::mutate(siteBag = as.factor(siteBag))

# relevel variable
 dat$siteBag<-forcats::fct_relevel(dat$siteBag, as.vector(unique(dat$siteBag)))

# pass data to list for JAGS
data = list(
  yv = as.double(dat$germCount),
  nv = as.double(dat$germStart),
  yv2 = as.double(dat$viabStain),
  nv2 = as.double(dat$viabStart),
  N = nrow(dat),
  bag = as.double(dat$siteBag),
  nbags = length(unique(dat$siteBag))
)

# set inits for JAGS
inits = list(list(p = rep(.1,data$nbags),p2 = rep(.1,data$nbags)),
             list(p = rep(.5,data$nbags),p2 = rep(.5,data$nbags)), 
             list(p = rep(.9,data$nbags),p2 = rep(.9,data$nbags)))

sink("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/jagsScripts/viabilityModelJAGS.R")
cat("
    model { 
    
    # theta ~ dunif(0,1)
    # kappa ~ dpar(alpha=shape,c=scale)
    # kappa ~ dpar(1.5,1)
    
    for(j in 1:nbags){
        p[j] ~ dbeta(1, 1)
        p2[j] ~ dbeta(1, 1)
    }

    for(i in 1:N){

    # v viability
    yv[i] ~ dbinom( p[bag[i]] , nv[i] )
    yv2[i] ~ dbinom( p2[bag[i]] , nv2[i] )
    
    yv.sim[i] ~ dbinom( p[bag[i]], nv[i] )
    yv2.sim[i] ~ dbinom( p2[bag[i]], nv2[i] )

    # # code for deviance from Lunn 2013
    prop[i] <- yv[i]/nv[i]
    prop2[i] <- yv2[i]/nv2[i]
    
    Ds[i] <- 2*nv[i]*(prop[i])*log((prop[i]+0.00001)/p[bag[i]]) + (1-prop[i])*log((1-prop[i]+0.00001)/(1-p[bag[i]]))
    Ds2[i] <- 2*nv2[i]*(prop2[i])*log((prop2[i]+0.00001)/p2[bag[i]]) + (1-prop2[i])*log((1-prop2[i]+0.00001)/(1-p2[bag[i]]))

   # sign[i] <- 2*step(prop[i] - p[bag[i]]) - 1
   # dev.res[i] <- sign[i]*sqrt(Ds[i])
    
    }
    
    # calculate saturated deviance
    dev.sat <- sum(Ds[])
    dev.sat2 <- sum(Ds2[])
    
    # derived quantity
    for(i in 1:nbags){
        vJoint[i] = p[i] + p2[i]*(1-p[i])
    }
    
    mean.data <- mean(yv)
    mean.sim <- mean(yv.sim)
    p.mean <- step(mean.sim - mean.data)
    
    }
    ", fill = TRUE)
sink()

# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 1000
n.update = 10000
n.iter = 10000
n.thin = 10

# Call to JAGS

dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/jagsScripts/")

# tuning (n.adapt)
jm = jags.model(paste0(dir,"viabilityModelJAGS.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

viab = c("p","p2","vJoint")

# chain (n.iter)
zc = coda.samples(jm, variable.names = c(viab), n.iter = n.iter, thin = n.thin)

MCMCsummary(zc, params = c("p","p2","vJoint"))

save(zc,file="/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/output/viabilityModelFit.rds")
save(data,file="/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/output/viabilityModelData.rds")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Posterior Predictive Checks

# notes: check group = interaction(dat$site,dat$yearStart)
# for stat density grouped plots
# -------------------------------------------------------------------
# -------------------------------------------------------------------

library(bayesplot)
library(gridExtra)

postCheck = c("p.mean","yv.sim","yv2.sim","dev.sat","dev.sat2")
zc = coda.samples(jm, variable.names = c(postCheck), n.iter = n.iter, thin = n.thin)

MCMCsummary(zc, params = c("p.mean","dev.sat","dev.sat2"))

# Conn et al. 2018 suggested saturated deviance should be roughly
# equal to the sample size - but what is OK for this?
hist(MCMCchains(zc,params="dev.sat")); abline(v=data$N,col="red")
hist(MCMCchains(zc,params="dev.sat2")); abline(v=data$N,col="red")

color_scheme_set("brightblue")

iter<-dim(MCMCchains(zc,params="yv.sim"))[1]

# germination trial
ppc_dens_overlay(data$yv, MCMCchains(zc,params="yv.sim")[sample(iter,1000), ]) +
  theme_bw() + xlim(c(0,15)) + labs(title="Posterior predictive checks for seeds counted in germination trials", 
                                     caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

# these would be very large plots!
# ppc_stat_grouped(data$yv, MCMCchains(zc,params="yv.sim")[sample(iter,1000), ],group=interaction(data$bag)) +
#   theme_bw() + labs(title="Posterior predictive checks for the mean of seeds counted germination trials", 
#                     caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")  

# not all trials have >1 replicate so variance is a moot point for those

# viability trial
ppc_dens_overlay(data$yv2, MCMCchains(zc,params="yv2.sim")[sample(iter,1000), ]) +
  theme_bw() + xlim(c(0,15)) + labs(title="Posterior predictive checks for seeds counted in viability trials", 
                                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

# these would be very large plots!
# ppc_stat_grouped(data$yv2, MCMCchains(zc,params="yv2.sim")[sample(iter,1000), ],group=interaction(dat$bagNo)) +
#   theme_bw() + labs(title="Posterior predictive checks for the mean of seeds counted viability trials", 
#                     caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")  


# observed vs. predicted
plot(data$yv,apply(MCMCchains(zc,params="yv.sim"),2,mean),
     xlim=c(0,15),ylim=c(0,15),
     pch=16,cex=.5);abline(a=0,b=1)
plot(data$yv2,apply(MCMCchains(zc,params="yv2.sim"),2,mean),
     xlim=c(0,15),ylim=c(0,15),
     pch=16,cex=.5);abline(a=0,b=1)

# observed vs. predicted
plot(apply(MCMCchains(zc,params="yv.sim"),2,mean),data$yv,
     xlim=c(0,15),ylim=c(0,15),
     pch=16,cex=.25);abline(a=0,b=1)
plot(apply(MCMCchains(zc,params="yv2.sim"),2,mean),data$yv2,
     xlim=c(0,15),ylim=c(0,15),
     pch=16,cex=.25);abline(a=0,b=1)

load("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/output/viabilityModelFit.rds")

# compare observed proportion vs. predicted probabilities
d<-apply(MCMCchains(zc,params="p"),2,mean) %>% as.data.frame() 
d<-data.frame(cbind(d,siteBag=unique(dat$siteBag))) %>% dplyr::rename(est='.')
d<-dat %>%
  dplyr::left_join(d,by="siteBag")

plot(d$est,data$yv/data$nv,
     xlim=c(0,1),ylim=c(0,1),
     pch=16,cex=.25);abline(a=0,b=1)

# compare observed proportion vs. predicted probabilities
d<-apply(MCMCchains(zc,params="p2"),2,mean) %>% as.data.frame() 
d<-data.frame(cbind(d,siteBag=unique(dat$siteBag))) %>% dplyr::rename(est='.')
d<-dat %>%
  dplyr::left_join(d,by="siteBag")

plot(d$est,data$yv2/data$nv2,
     xlim=c(0,1),ylim=c(0,1),
     pch=16,cex=.25);abline(a=0,b=1)
