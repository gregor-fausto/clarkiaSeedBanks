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


# one row is coded differently so that 
# viabStart=NA and viabStain=NA
# all others have viabStart=NA and viabStain=NA
# recode
df[is.na(df$viabStart),]$viabStart = 0

## data check
df %>% dplyr::filter(germStart - germCount - viabStart<0) 

# filter out rows with problems
# these need to be corrected
df<-df %>% 
  dplyr::filter(germStart - germCount - viabStart >= 0)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Prepare data for analysis with JAGS
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
  yg = as.double(dat$germCount),
  ng = as.double(dat$germStart),
  yv = as.double(dat$viabStain),
  nv = as.double(dat$viabStart),
  N = nrow(dat),
  bag = as.double(dat$siteBag),
  nbags = length(unique(dat$siteBag))
)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 1000
n.update = 10000
n.iter = 10000
n.thin = 10

set.seed(10)
dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/jagsScripts/")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Complete pooling
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# set inits for JAGS
inits = list(list(pv = .1,pg = .1), 
             list(pv = .5,pg = .5), 
             list(pv = .9,pg = .9))

# Call to JAGS

# tuning (n.adapt)
jm = jags.model(paste0(dir,"viabilityJointCompletePoolingJAGS.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

viab = c("pv","pg","viability")

# chain (n.iter)
zc_pool = coda.samples(jm, variable.names = c(viab), n.iter = n.iter, thin = n.thin)

MCMCsummary(zc_pool, params = c("pv","pg","viability"))
save(zc_pool,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/viabilityJointCompletePoolFit.rds")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# No pooling: each trial is independent
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# set inits for JAGS
inits = list(list(pv = rep(.1,data$N), pg = rep(.1,data$N)),
             list(pv = rep(.5,data$N), pg = rep(.5,data$N)),
             list(pv = rep(.9,data$N), pg = rep(.9,data$N)))

# Call to JAGS

# tuning (n.adapt)
jm = jags.model(paste0(dir,"viabilityJointNoPoolingJAGS.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

viab = c("pv", "pg", "viability")

# chain (n.iter)
zc_nopool = coda.samples(jm, variable.names = c(viab), n.iter = n.iter, thin = n.thin)

MCMCsummary(zc_nopool, params = c("pv", "pg", "viability"))
save(zc_nopool,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/viabilityJointNoPoolFit.rds")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Partial pooling
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# set inits for JAGS
inits = list(list(pv = rep(.1,data$nbags), pg = rep(.1,data$nbags)),
             list(pv = rep(.5,data$nbags), pg = rep(.5,data$nbags)),
             list(pv = rep(.9,data$nbags), pg = rep(.9,data$nbags)))

# Call to JAGS

# tuning (n.adapt)
jm = jags.model(paste0(dir,"viabilityJointPartialPoolingJAGS.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

viab = c("pv", "pg","viability")

# chain (n.iter)
zc_partialpool = coda.samples(jm, variable.names = c(viab), n.iter = n.iter, thin = n.thin)

MCMCsummary(zc_partialpool, params = c("pv", "pg","viability"))
save(zc_partialpool,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/viabilityJointPartialPoolFit.rds")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Partial pooling: log odds
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# set inits for JAGS
inits = list(
  list( sigma.v = 50, sigma.g = 50,
        mu.v = 0, mu.g = 0) ,
  list( sigma.v = 20, sigma.g = 20, 
        mu.v = 0, mu.g = 0)  ,
  list( sigma.v = 10, sigma.g = 10,
        mu.v = 0, mu.g = 0)  
) 


# Call to JAGS

# tuning (n.adapt)
jm = jags.model(paste0(dir,"viabilityJointPartialPoolingLogsJAGS.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

viab = c("pv", "pg", "viability")

# chain (n.iter)
zc_partialpoollogs = coda.samples(jm, variable.names = c(viab), n.iter = n.iter, thin = n.thin)

MCMCsummary(zc_partialpoollogs, params = c("pv","pg","viability"))
save(zc_partialpoollogs,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/viabilityJointPartialPoolLogsFit.rds")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Partial pooling: hyperpriors
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# set inits for JAGS
inits = list(
  list( theta.v = .1, kappa.v = 1.1, 
        theta.g = .1, kappa.g = 1.1) ,
  list( theta.v = .5, kappa.v = 1.5,
        theta.g = .5, kappa.g = 1.5)  ,
  list( theta.v = .9, kappa.v = 2,
        theta.g = .9, kappa.g = 2)  
) 


# Call to JAGS

# tuning (n.adapt)
jm = jags.model(paste0(dir,"viabilityJointPartialPoolingHyperpriorsJAGS.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

viab = c("pv", "pg", "viability")

# chain (n.iter)
zc_partialpoolhyper = coda.samples(jm, variable.names = c(viab), n.iter = n.iter, thin = n.thin)

MCMCsummary(zc_partialpoolhyper, params = c("pv","pg","viability"))
save(zc_partialpoolhyper,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/viabilityJointPartialPoolHyperpriorsFit.rds")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Model comparison
# -------------------------------------------------------------------
# -------------------------------------------------------------------

x<-list.files("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/")

x<-x[sapply(x,stringr::str_detect,pattern="Fit")]
x<-x[sapply(x,stringr::str_detect,pattern="Pool")]
x<-x[sapply(x,stringr::str_detect,pattern="Joint")]

x<-paste0("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/",x)

lapply(x,load,.GlobalEnv)

library(MCMCvis)

N = dim(dat)[1]

focalParameter="viability"

ss_pool = MCMCchains(zc_pool,params="viability")
ss_nopool = MCMCchains(zc_nopool,params="viability")
ss_partialpool = MCMCchains(zc_partialpool,params="viability")
ss_partialpoollogs = MCMCchains(zc_partialpoollogs,params="viability")
ss_partialpoolhyper = MCMCchains(zc_partialpoolhyper,params="viability")

rm(zc_pool,zc_nopool,zc_partialpool,zc_partialpoollogs,zc_partialpoolhyper)

#
p_pool <- t(apply(ss_pool,2,function(x) quantile(x,probs=c(.1,.5,.9)))) %>%
  data.frame %>% 
  dplyr::rename(theta_50 = X50.,theta_10 = X10., theta_90 = X90.)

#
p_nopool <- t(apply(ss_nopool,2,function(x) quantile(x,probs=c(.1,.5,.9)))) %>%
  data.frame %>% 
  dplyr::rename(theta_50 = X50.,theta_10 = X10., theta_90 = X90.)

#
p_partialpool <- t(apply(ss_partialpool,2,function(x) quantile(x,probs=c(.1,.5,.9)))) %>%
  data.frame %>% 
  dplyr::rename(theta_50 = X50.,theta_10 = X10., theta_90 = X90.)

p_partialpool <- p_partialpool %>% 
  bind_cols(bagNoUnique = 1:dim(p_partialpool)[1])

p_partialpool<- dat %>%
  dplyr::bind_cols(bagNoUnique = data$bag) %>%
  dplyr::left_join(p_partialpool,by="bagNoUnique") 

#
p_partialpoollogs <- t(apply(ss_partialpoollogs,2,function(x) quantile(x,probs=c(.1,.5,.9)))) %>%
  data.frame %>% 
  dplyr::rename(theta_50 = X50.,theta_10 = X10., theta_90 = X90.)

p_partialpoollogs <- p_partialpoollogs %>% 
  bind_cols(bagNoUnique = 1:dim(p_partialpoollogs)[1])

p_partialpoollogs<- dat %>%
  dplyr::bind_cols(bagNoUnique = data$bag) %>%
  dplyr::left_join(p_partialpoollogs,by="bagNoUnique") 
#
p_partialpoolhyper <- t(apply(ss_partialpoolhyper,2,function(x) quantile(x,probs=c(.1,.5,.9)))) %>%
  data.frame %>% 
  dplyr::rename(theta_50 = X50.,theta_10 = X10., theta_90 = X90.)

p_partialpoolhyper <- p_partialpoolhyper %>% 
  bind_cols(bagNoUnique = 1:dim(p_partialpoolhyper)[1])

p_partialpoolhyper<- dat %>%
  dplyr::bind_cols(bagNoUnique = data$bag) %>%
  dplyr::left_join(p_partialpoolhyper,by="bagNoUnique")

lapply(list(p_pool,p_nopool,p_partialpool,p_partialpoollogs,p_partialpoolhyper),dim)

plot(p_partialpool$theta_50,p_partialpoollogs$theta_50);abline(a=0,b=1)
plot(p_partialpool$theta_50,p_partialpoolhyper$theta_50);abline(a=0,b=1)

df_plot2<-data.frame(group = rep(data$bag, 5),
                     model = c(rep("complete pooling", N),
                               rep("no pooling", N),
                               rep("partial pooling", N),
                               rep('partial pooling, logit', N),
                               rep("partial pooling, hyperpriors", N)),
                     y = c(rep(p_pool$theta_50,N),
                           p_nopool$theta_50,
                           p_partialpool$theta_50,
                           p_partialpoollogs$theta_50,
                           p_partialpoolhyper$theta_50))

dfPlot <- df_plot2 %>% 
  group_by(model) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = model, values_from = y) %>% 
  select(-row) %>%
  bind_cols(site=dat$site)

df_plotlo<-data.frame(group = rep(data$bag, 5),
                     model = c(rep("complete pooling", N),
                               rep("no pooling", N),
                               rep("partial pooling", N),
                               rep('partial pooling, logit', N),
                               rep("partial pooling, hyperpriors", N)),
                     y = c(rep(p_pool$theta_10,N),
                           p_nopool$theta_10,
                           p_partialpool$theta_10,
                           p_partialpoollogs$theta_10,
                           p_partialpoolhyper$theta_10))

dfPlotLo <- df_plotlo %>% 
  group_by(model) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = model, values_from = y) %>% 
  select(-row)


df_plothi<-data.frame(group = rep(data$bag, 5),
                      model = c(rep("complete pooling", N),
                                rep("no pooling", N),
                                rep("partial pooling", N),
                                rep('partial pooling, logit', N),
                                rep("partial pooling, hyperpriors", N)),
                      y = c(rep(p_pool$theta_90,N),
                            p_nopool$theta_90,
                            p_partialpool$theta_90,
                            p_partialpoollogs$theta_90,
                            p_partialpoolhyper$theta_90))

dfPlotHi <- df_plothi %>% 
  group_by(model) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = model, values_from = y) %>% 
  select(-row) 

library(GGally)
pairsPlot <- ggpairs(dfPlot, columns = 3:6, 
                     ggplot2::aes(alpha = 0.4,color=as.factor(site)),
                     upper = list(continuous = wrap("cor", size = 2)))
pairsPlot


ggsave(plot=pairsPlot,
       filename="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/viabilityPairsPlot.png",
       width=12,height=12)

# ggplot() + 
#   geom_point(data = dfPlot,aes(x = `no pooling`,y = `partial pooling`)) + 
#   geom_errorbar(aes(x= dfPlot$`no pooling`, ymin = dfPlotLo$`partial pooling`,ymax = dfPlotHi$`partial pooling`)) +
#   geom_errorbarh(aes(y= dfPlot$`partial pooling`, xmin = dfPlotLo$`no pooling`,xmax = dfPlotHi$`no pooling`))
# 
# ggplot() + 
#   geom_point(data = dfPlot,aes(x = `partial pooling`,y = `partial pooling, logit`)) + 
#   geom_errorbar(aes(x= dfPlot$`partial pooling`, ymin = dfPlotLo$`partial pooling, logit`,ymax = dfPlotHi$`partial pooling, logit`)) +
#   geom_errorbarh(aes(y= dfPlot$`partial pooling, logit`, xmin = dfPlotLo$`partial pooling`,xmax = dfPlotHi$`partial pooling`))
# 
# ggplot() + 
#   geom_point(data = dfPlot,aes(x = `partial pooling`,y = `partial pooling, hyperpriors`)) + 
#   geom_errorbar(aes(x= dfPlot$`partial pooling`, ymin = dfPlotLo$`partial pooling, hyperpriors`,ymax = dfPlotHi$`partial pooling, hyperpriors`)) +
#   geom_errorbarh(aes(y= dfPlot$`partial pooling, hyperpriors`, xmin = dfPlotLo$`partial pooling`,xmax = dfPlotHi$`partial pooling`))


