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
# Summarize viability trial data 
# -------------------------------------------------------------------
# 
# df <- df %>%
#   dplyr::mutate(n_new = germCount+viabStart,
#                 y_new = germCount+viabStain) %>%
#   dplyr::group_by(site,round,age,bag) %>%
#   dplyr::summarise(y_new=sum(y_new),n_new=sum(n_new))

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Multiple-site model (3 year)
# Log-odds parameterization for viability, trials summed across bags
# -------------------------------------------------------------------
# -------------------------------------------------------------------

## AGE 1 SEEDS
dat<-df %>% 
  dplyr::filter(age==1) %>%
  dplyr::filter(site=="BG")

# remove missing data for now
# learn how to model these
# dat <- dat[!is.na(dat$n_new),]

dat$site.index <- as.integer(as.factor(dat$site))
# dat$year.index <- as.integer(as.factor(dat$yearStart))

# having problem converting site-bag number to numeric properly for JAGS call
dat<-dat %>% tidyr::unite(siteBag,c("site","bagNo"),sep="")
as.numeric(dat$siteBag)
as.numeric("a")

# pass data to list for JAGS
data = list(
  yv = as.double(dat$germCount),
  nv = as.double(dat$germStart),
  N = nrow(dat),
  bag = as.double(as.factor(paste0(dat$site,dat$bagNo))),
  nbags = length(unique(as.factor(paste0(dat$site,dat$bagNo)))),
  site = as.double(dat$site.index),
  nsites = length(unique(dat$site.index)),
  year = as.double(dat$round),
  nyears = length(unique(dat$round))
)

# right now each bag has its own prior;
# want to give nsite priors 
# inits = list(
#   list( sigma.b = matrix(rep(50,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears), 
#         mu.b = matrix(rep(0,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears)
#   ),
#   list( sigma.b = matrix(rep(10,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears), 
#         mu.b = matrix(rep(0,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears)
#   ),
#   list( sigma.b = matrix(rep(20,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears), 
#         mu.b = matrix(rep(0,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears)
#   )
# )

inits = list(list(p = rep(.1,data$nbags)),list(p = rep(.5,data$nbags)), list(p = rep(.9,data$nbags)))

#setwd("Users/Gregor/Dropbox/clarkiaSeedBags/modelBuild/jagsScripts")
sink("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/jagsScripts/viabilityModelJAGS.R")
cat("
    model { 
    
    for(j in 1:nbags){
        p[j] ~ dbeta(1, 1)
    }

    for(i in 1:N){

    # v viability
    yv[i] ~ dbinom( p[bag[i]] , nv[i])
    yv.sim[i] ~ dbinom(p[bag[i]], nv[i]) 
    
    }
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

viab = c("p","yv.sim")

# chain (n.iter)
zc = coda.samples(jm, variable.names = c(viab), n.iter = n.iter, thin = n.thin)

MCMCsummary(zc, params = c("p","yv.sim"))

bayesplot::ppc_dens_overlay(data$yv, MCMCchains(zc,params="yv.sim")[sample((n.iter*length(inits))/n.thin,1000), ]) +
  theme_bw() + xlim(c(0,15)) + labs(title="Posterior predictive checks for seeds counted in seed bags in January", 
                                     caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

bayesplot::ppc_stat_grouped(data$yv, MCMCchains(zc,params="yv.sim")[sample((n.iter*length(inits))/n.thin,1000), ],group=interaction(dat$bagNo)) +
  theme_bw() + labs(title="Posterior predictive checks for the mean of seeds counted in seed bags in January", 
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")  

save(zc,file="/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/output/viabilityModelFit.rds")
save(data,file="/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/output/viabilityModelData.rds")
