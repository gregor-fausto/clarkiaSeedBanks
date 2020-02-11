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

df <- df %>%
  dplyr::mutate(n_new = germCount+viabStart,
                y_new = germCount+viabStain) %>%
  dplyr::group_by(site,round,age,bag) %>%
  dplyr::summarise(y_new=sum(y_new),n_new=sum(n_new))

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Multiple-site model (3 year)
# Log-odds parameterization for viability, trials summed across bags
# -------------------------------------------------------------------
# -------------------------------------------------------------------

## AGE 1 SEEDS
dat<-df %>% 
  dplyr::filter(age==1)

# remove missing data for now
# learn how to model these
dat <- dat[!is.na(dat$n_new),]

dat$site.index <- as.integer(as.factor(dat$site))
# dat$year.index <- as.integer(as.factor(dat$yearStart))

# pass data to list for JAGS
data = list(
  yv = as.double(dat$y_new),
  nv = as.double(dat$n_new),
  N = nrow(dat),
  site = as.double(dat$site.index),
  nsites = length(unique(dat$site.index)),
  year = as.double(dat$round),
  nyears = length(unique(dat$round))
)

# right now each bag has its own prior;
# want to give nsite priors 
inits = list(
  list( sigma.b = matrix(rep(50,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears), 
        mu.b = matrix(rep(0,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears)
        ),
  list( sigma.b = matrix(rep(10,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears), 
        mu.b = matrix(rep(0,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears)
        ),
  list( sigma.b = matrix(rep(20,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears), 
        mu.b = matrix(rep(0,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears)
        )
)

#setwd("Users/Gregor/Dropbox/clarkiaSeedBags/modelBuild/jagsScripts")
sink("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/jagsScripts/viabilityModelJAGS.R")
cat("
    model { 
    
    ##############
    ## hyperpriors
    ##############

    # viability trials
    for(j in 1:nsites){
    for(k in 1:nyears){ 
    mu.b[j,k] ~ dnorm(0,.0001)
    
    sigma.b[j,k] ~ dunif(0,100)
    tau.b[j,k] <- 1/(sigma.b[j,k] * sigma.b[j,k])
    }
    }


    ##############
    ## priors
    ##############
    
    # viability trials 1 prior for each trial
    for(i in 1:N){
    alpha[i] ~ dnorm(mu.b[site[i],year[i]],  tau.b[site[i],year[i]])
    }

    ##############
    ## likelihoods
    ##############

    for(i in 1:N){

    # v viability
    p[i] <- ilogit(alpha[i])
    yv[i] ~ dbin( p[i] , nv[i])
    # yv.sim[i] ~ dbinom(p[i], nv[i]) 
    
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

# Call to JAGS

dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/jagsScripts/")

# tuning (n.adapt)
jm = jags.model(paste0(dir,"viabilityModelJAGS.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

viab = c("mu.b","sigma.b")

# chain (n.iter)
zc = coda.samples(jm, variable.names = c(viab), n.iter = n.iter, thin=10)

save(zc,file="/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/output/viabilityModelFit.rds")
save(data,file="/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/output/viabilityModelData.rds")
