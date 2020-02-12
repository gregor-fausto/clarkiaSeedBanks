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
# Import and organize seed bag data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# setwd and read data files
load("/Users/Gregor/Dropbox/clarkiaSeedBanks/library/dataFromWorkflowFile/seedBagsData.rda")
df <- seedBags

# -------------------------------------------------------------------
# Clean and organize seed bag data
# -------------------------------------------------------------------
df$seedStart<-100

df$seedStart <- as.integer(df$seedStart)

df$totalJan <- ifelse(df$totalJan>100,100,df$totalJan)

df <- df %>% dplyr::rename(bag=bagNo)

## FOR NOW REMOVE MISSING DATA AND PROBLEMS
df<-subset(df,!is.na(df$totalJan))
df<-subset(df,!is.na(df$intactOct))
df<-subset(df,!(intactJan<intactOct))


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize viability trial data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# setwd and read data files
setwd("~/Dropbox/modelsF2019/viability/")
df2 <- read.csv(file="viability.csv",header=TRUE)

# another check
seedBurial <- df %>% dplyr::select(site,bag,round,age) %>%
  dplyr::mutate(seedBurial = 1)
viabilityTrial <- df2 %>% dplyr::select(site,bag,round,age) %>%
  dplyr::mutate(viabilityTrial = 1)

joinedDF<-seedBurial %>% 
  dplyr::full_join(viabilityTrial) 
  
dplyr::group_by(site,round,age) %>%
  dplyr::count(seedBurial,viabilityTrial)
    
  joinedDF %>% 
    dplyr::filter(site=="BR"&round==1&age==2)

  # wy are there seed bags that show up in the seed burial experiment
  # but not in the viability trial experiment?
  
df2 <- df2 %>% dplyr::select(-c(germPerc,germNot,viabPerc,viabPerc2,condTest))
df2$bag <-as.integer(as.numeric(df2$bag))

## FOR NOW REMOVE MISSING DATA AND PROBLEMS
df2<-subset(df2,!is.na(df2$germStart))
df2<-subset(df2,!is.na(df2$germCount))
df2<-subset(df2,!is.na(df2$viabStart))
df2<-subset(df2,!is.na(df2$viabStain))

## data check
df2 %>% dplyr::filter(germStart - germCount - viabStart<0)
#
# -------------------------------------------------------------------
# Summarize viability trial data 
# -------------------------------------------------------------------

df2 <- df2 %>%
  dplyr::mutate(n_new = germCount+viabStart,
                y_new = germCount+viabStain) %>%
  dplyr::group_by(site,round,age,bag) %>%
  dplyr::summarise(y_new=sum(y_new),n_new=sum(n_new))

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Join datasets
# not sure if it's necessary to join but not sure how else to index
# -------------------------------------------------------------------
# -------------------------------------------------------------------
df_joined <- df %>%
  dplyr::left_join(df2,by=c("site","bag","round","age"))

# df_joined<-df_joined %>%
#   dplyr::group_by(site,round,age,bag) 

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Multiple-site model (3 year)
# Log-odds parameterization for viability, trials summed across bags
# -------------------------------------------------------------------
# -------------------------------------------------------------------

## AGE 1 SEEDS
dat<-df_joined %>% 
  #dplyr::filter(round=="1") %>% 
  dplyr::filter(age==1)

# remove missing data for now
# learn how to model these
dat <- dat[!is.na(dat$n_new),]

names(dat)[1] <- "site"
dat$site.index <- as.integer(as.factor(dat$site))
dat$year.index <- as.integer(as.factor(dat$yearStart))

## AGE 2 SEEDS
datTwo<-df_joined %>% 
  #dplyr::filter(round=="1") %>% 
  dplyr::filter(age==2)

# remove missing data for now
# learn how to model these
datTwo <- datTwo[!is.na(datTwo$n_new),]

names(datTwo)[1] <- "site"
datTwo$site.index <- as.integer(as.factor(datTwo$site))
datTwo$year.index <- as.integer(as.factor(datTwo$yearStart))


# pass data to list for JAGS
data = list(
  yg = as.double(dat$seedlingJan),
  yt = as.double(dat$totalJan),
  yo = as.double(dat$intactOct),
  yv = as.double(dat$y_new),
  n = as.double(dat$seedStart),
  nv = as.double(dat$n_new),
  N = nrow(dat),
  site = as.double(dat$site.index),
  nsites = length(unique(dat$site.index)),
  year = as.double(dat$year.index),
  nyears = length(unique(dat$year.index)),
  
  yt2 = as.double(datTwo$totalJan),
  n2 = as.double(datTwo$seedStart),
  N2 = nrow(datTwo),
  site2 = as.double(datTwo$site.index),
  nsites2 = length(unique(datTwo$site.index)),
  year2 = as.double(datTwo$year.index),
  nyears2 = length(unique(datTwo$year.index))
)

# right now each bag has its own prior;
# want to give nsite priors 
inits = list(
  list( sigma.b = matrix(rep(50,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears), 
        mu.b = matrix(rep(0,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears), 
        #alpha = rep(-10, data$N),
        # sigmaS1 = rep(.1, data$nsites), 
        # sigmaG1 = rep(.1, data$nsites), 
        # sigmaS2 = rep(.1, data$nsites),
        alphaS1 = rep(-10, data$nsites), 
        alphaG1 = rep(-10, data$nsites), 
        alphaS2 = rep(-10, data$nsites),
        alphaS3 = rep(0, data$nsites)
        # betaS1 = matrix(-10, nrow=data$nsites,ncol=data$nyears), 
        # betaG1 = matrix(-10, nrow=data$nsites,ncol=data$nyears), 
        # betaS2 = matrix(-10, nrow=data$nsites,ncol=data$nyears)
        ),
  list( sigma.b = matrix(rep(10,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears), 
        mu.b = matrix(rep(0,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears),
        #alpha = rep(0, data$N),
        # sigmaS1 = rep(.5, data$nsites), 
        # sigmaG1 = rep(.5, data$nsites), 
        # sigmaS2 = rep(.5, data$nsites),
        alphaS1 = rep(0, data$nsites), 
        alphaG1 = rep(0, data$nsites), 
        alphaS2 = rep(0, data$nsites),
        alphaS3 = rep(0, data$nsites)
        # betaS1 = matrix(0, nrow=data$nsites,ncol=data$nyears), 
        # betaG1 = matrix(0, nrow=data$nsites,ncol=data$nyears), 
        # betaS2 = matrix(0, nrow=data$nsites,ncol=data$nyears)
        ),
  list( sigma.b = matrix(rep(20,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears), 
        mu.b = matrix(rep(0,data$nsites*data$nyears),nrow=data$nsites,ncol=data$nyears), 
        #alpha = rep(10, data$N),
        #sigmaS1 = rep(.9, data$nsites), 
        #sigmaG1 = rep(.9, data$nsites), 
        #sigmaS2 = rep(.9, data$nsites),
        alphaS1 = rep(10, data$nsites), 
        alphaG1 = rep(10, data$nsites), 
        alphaS2 = rep(10, data$nsites),
        alphaS3 = rep(0, data$nsites)
        # betaS1 = matrix(10, nrow=data$nsites,ncol=data$nyears), 
        # betaG1 = matrix(10, nrow=data$nsites,ncol=data$nyears), 
        # betaS2 = matrix(10, nrow=data$nsites,ncol=data$nyears)
        )
)

#setwd("Users/Gregor/Dropbox/clarkiaSeedBags/modelBuild/jagsScripts")
sink("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/jagsScripts/seed_modelJAGS.R")
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

    # # seed bags
    # for(j in 1:nsites){
    # sigmaS1[j] ~ dunif(0,100)
    # tauS1[j] <- 1/(sigmaS1[j] * sigmaS1[j])
    # 
    # sigmaG1[j] ~ dunif(0,100)
    # tauG1[j] <- 1/(sigmaG1[j] * sigmaG1[j])
    # 
    # sigmaS2[j] ~ dunif(0,100)
    # tauS2[j] <- 1/(sigmaS2[j] * sigmaS2[j])
    # }

    ##############
    ## priors
    ##############
    
    # site intercepts
    for(j in 1:nsites){
    alphaS1[j] ~ dnorm(0, .001)
    alphaG1[j] ~ dnorm(0, .001)
    alphaS2[j] ~ dnorm(0, .001)
    alphaS3[j] ~ dnorm(0, .001)
    }

    # year intercepts
    # for(j in 1:nsites){
    # for(k in 1:nyears){
    # betaS1[j,k] ~ dnorm(0, tauS1[j])
    # betaG1[j,k] ~ dnorm(0, tauG1[j])
    # betaS2[j,k] ~ dnorm(0, tauS2[j])
    # }
    # }
    # 
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

    # s1 seed survival
    ps[i] = ilogit(alphaS1[site[i]])
    yt[i] ~ dbin(ps[i], n[i])
    # yt.sim[i] ~ dbinom(ps[i], n[i]) 

    # g1 seed germination
    pg[i] = ilogit(alphaG1[site[i]])
    yg[i] ~ dbin(pg[i]*(p[i])^(1/3), yt[i])
    # yg.sim[i] ~ dbinom(pg[i]*(p[i])^(1/3), yt[i]) 

    # s2 seed survival
    pr[i] = ilogit(alphaS2[site[i]])
    yo[i] ~ dbin(pr[i], yt[i]-yg[i])
    # yo.sim[i] ~ dbinom(pr[i], yt[i]-yg[i]) 

    } 

    for(i in 1:N2){
    
    ps2[i] = ilogit(alphaS1[site[i]])
    pg2[i] = ilogit(alphaG1[site[i]])
    pr2[i] = ilogit(alphaS2[site[i]])

    # s3 seed survival
    ps3[i] = ilogit(alphaS3[site2[i]])
    yt2[i] ~ dbin(ps2[i]*(1-pg2[i])*pr2[i]*ps3[i], n2[i])
    # yt2.sim[i] ~ dbinom(ps2[i]*(1-pg2[i])*pr2[i]*ps3[i], n2[i])
    
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
jm = jags.model(paste0(dir,"seed_modelJAGS.R"), data = data, inits = inits,
                n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

intercepts = c("alphaS1","alphaG1","alphaS2","alphaS3")
#slopes = c("betaS1","betaG1","betaS2")
#variances = c("sigmaS1","sigmaG1","sigmaS2")
viab = c("mu.b","sigma.b")
#sims = c("yv.sim","yt.sim","yg.sim","yo.sim","yt2.sim")



# chain (n.iter)
zc = coda.samples(jm, variable.names = c(intercepts,viab), n.iter = n.iter, thin=10)

save(zc,file="/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/output/seedbagfit.rds")
save(data,file="/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/output/seedbagdata.rds")
