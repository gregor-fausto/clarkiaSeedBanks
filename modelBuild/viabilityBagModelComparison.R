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
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# AGE 1 SEEDS
dat<-df %>% 
  dplyr::filter(age==1) %>%
  dplyr::filter(site=="BG")

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
inits = list(list(p = .1), list(p = .5), list(p = .9))

# Call to JAGS

# tuning (n.adapt)
jm = jags.model(paste0(dir,"viabilityBagsCompletePoolingJAGS.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

viab = c("p")

# chain (n.iter)
zc_pool = coda.samples(jm, variable.names = c(viab), n.iter = n.iter, thin = n.thin)

MCMCsummary(zc_pool, params = c("p"))
save(zc_pool,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/viabilityCompletePoolFit.rds")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# No pooling: each trial is independent
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# set inits for JAGS
inits = list(list(p = rep(.1,data$N)),
                  list(p = rep(.5,data$N)),
                  list(p = rep(.9,data$N)))

# Call to JAGS

# tuning (n.adapt)
jm = jags.model(paste0(dir,"viabilityBagsNoPoolingJAGS.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

viab = c("p")

# chain (n.iter)
zc_nopool = coda.samples(jm, variable.names = c(viab), n.iter = n.iter, thin = n.thin)

MCMCsummary(zc_nopool, params = c("p"))
save(zc_nopool,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/viabilityNoPoolFit.rds")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Partial pooling
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# set inits for JAGS
inits = list(list(p = rep(.1,data$nbags)),
             list(p = rep(.5,data$nbags)),
             list(p = rep(.9,data$nbags)))

# Call to JAGS

# tuning (n.adapt)
jm = jags.model(paste0(dir,"viabilityBagsPartialPoolingJAGS.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

viab = c("p")

# chain (n.iter)
zc_partialpool = coda.samples(jm, variable.names = c(viab), n.iter = n.iter, thin = n.thin)

MCMCsummary(zc_partialpool, params = c("p"))
save(zc_partialpool,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/viabilityPartialPoolFit.rds")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Partial pooling: log odds
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# set inits for JAGS
inits = list(
  list( sigma = 50, 
        mu = 0) ,
  list( sigma = 20, 
        mu = 0)  ,
  list( sigma = 10, 
        mu = 0)  
) 


# Call to JAGS

# tuning (n.adapt)
jm = jags.model(paste0(dir,"viabilityBagsPartialPoolingLogsJAGS.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

viab = c("p")

# chain (n.iter)
zc_partialpoollogs = coda.samples(jm, variable.names = c(viab), n.iter = n.iter, thin = n.thin)

MCMCsummary(zc_partialpoollogs, params = c("p"))
save(zc_partialpoollogs,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/viabilityPartialPoolLogsFit.rds")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Partial pooling: hyperpriors
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# set inits for JAGS
inits = list(
  list( theta = .1, 
        kappa = 1.1) ,
  list( theta = .5, 
        kappa = 1.5)  ,
  list( theta = .9, 
        kappa = 2)  
) 


# Call to JAGS

# tuning (n.adapt)
jm = jags.model(paste0(dir,"viabilityBagsPartialPoolingHyperpriorsJAGS.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

viab = c("p")

# chain (n.iter)
zc_partialpoolhyper = coda.samples(jm, variable.names = c(viab), n.iter = n.iter, thin = n.thin)

MCMCsummary(zc_partialpoolhyper, params = c("p"))
save(zc_partialpoollogs,file="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/viabilityPartialPoolHyperpriorsFit.rds")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Model comparison
# -------------------------------------------------------------------
# -------------------------------------------------------------------

N = dim(dat)[1]

ss_pool = MCMCchains(zc_pool,params="p")
ss_nopool = MCMCchains(zc_nopool,params="p")
ss_partialpool = MCMCchains(zc_partialpool,params="p")
ss_partialpoollogs = MCMCchains(zc_partialpoollogs,params="p")
ss_partialpoolhyper = MCMCchains(zc_partialpoolhyper,params="p")

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

#
p_partialpoolhyper <- t(apply(ss_partialpoolhyper,2,function(x) quantile(x,probs=c(.1,.5,.9)))) %>%
  data.frame %>% 
  dplyr::rename(theta_50 = X50.,theta_10 = X10., theta_90 = X90.)

p_partialpoolhyper <- p_partialpoolhyper %>% 
  bind_cols(bagNoUnique = 1:dim(p_partialpoolhyper)[1])

p_partialpoolhyper<- dat %>%
  dplyr::bind_cols(bagNoUnique = data$bag) %>%
  dplyr::left_join(p_partialpoolhyper,by="bagNoUnique") 

df_plot2<-data.frame(x = rep(data$yv / data$nv, 5),
                     group = rep(data$bag, 5),
           y = c(rep(p_pool$theta_50,N),
                 p_nopool$theta_50,
                 p_partialpool$theta_50,
                 p_partialpoollogs$theta_50,
                 p_partialpoolhyper$theta_50),
           model = c(rep("complete pooling", N),
                     rep("no pooling", N),
                     rep("partial pooling", N),
                     rep('partial pooling, logit', N),
                     rep("partial pooling, hyperpriors", N)))

pop_mean <- sum(data$yv) / sum(data$nv);

# expand color palette to see pooling
library(pals)
pal.bands(alphabet, alphabet2, cols25, glasbey, kelly, polychrome, 
          stepped, tol, watlington,
          show.names=FALSE)

plot_bda3_fig_5_4 <-
  ggplot(df_plot2, aes(x=x, y=y)) +
  geom_hline(aes(yintercept=pop_mean), colour="lightpink") +
  geom_abline(intercept=0, slope=1, colour="skyblue") +
  facet_grid(. ~ model) +
  geom_errorbar(aes(ymin=c(rep(p_pool$theta_10,N),
                           p_nopool$theta_10,
                           p_partialpool$theta_10,
                           p_partialpoollogs$theta_10,
                           p_partialpoolhyper$theta_10),
                    ymax=c(rep(p_pool$theta_90,N),
                           p_nopool$theta_90,
                           p_partialpool$theta_90,
                           p_partialpoollogs$theta_90,
                           p_partialpoolhyper$theta_90)),
                width=0.005, colour="gray60") +
  geom_point(aes(color=as.factor(group)), size=0.75) +
  scale_x_continuous(breaks = seq(0,1,by=.1)) +
  xlab("observed proportion, y[n] / K[n]") +
  ylab("chance of success, theta[n]") +
  theme_bw() +
  scale_color_manual(values=as.vector(alphabet(26))) +
  ggtitle("Posterior Medians and 80% intervals\n(red line: population mean;  blue line: MLE)")
plot_bda3_fig_5_4


