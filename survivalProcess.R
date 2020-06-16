# -------------------------------------------------------------------
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

set.seed(17)

samples = 1000
n <- rpois(n = samples,8)

arm::discrete.histogram(n)

psi <- rbinom(n = samples, size = 1, prob = .5)

n2 <- rpois(n = samples, 5)

n3 <- psi*n2

z <- n+n3

arm::discrete.histogram(z)

y <- rbinom(n = samples, size = z, prob = .25)

data =data.frame(seedlingNumber = n, undercount = psi, z = z, fruitplNumber = y)
data

theta_obs = y/n
theta_true = y/z

plot(theta_true,theta_obs)

mean(theta_obs);sum(y)/sum(n);mean(y[theta_obs<=1]/n[theta_obs<=1]);sum(y[theta_obs<=1])/sum(n[theta_obs<=1])
mean(theta_true);sum(y)/sum(z)
data

data

## read in data

# # setwd and read data files
# censusSeedlingsFruitingPlants <- readRDS("~/Dropbox/dataLibrary/postProcessingData/censusSeedlingsFruitingPlants.RDS")
# 
# data = censusSeedlingsFruitingPlants %>%
#   #dplyr::filter(site=="BG") %>%
#   dplyr::mutate(n = seedlingNumber, y = fruitplNumber) %>%
#   dplyr::mutate(undercount = ifelse(y>n,1,0)) %>%
#   dplyr::group_by(site,year) %>%
#   dplyr::summarise(count = sum(undercount,na.rm=TRUE))
# 
# # par(mfrow=c(3,5)) 
# # for(i in 1:14){
# #   hist(data[data$year==(2005+i),]$theta_obs_i,main=2005+i)
# # }
# 
# # data %>%
# #   tidyr::drop_na(y, n) %>%
# #   dplyr::group_by(year) %>%
# #   dplyr::summarise(y = sum(y), n = sum(n), theta_obs = y/n)
# 
# #data[data$year==2018,] %>% View
# 
# 
# data <- data %>% dplyr::filter(year==2007)
# data=data %>% dplyr::select(seedlingNumber,fruitplNumber,undercount)

# -------------------------------------------------------------------
# Clean and organize data
# -------------------------------------------------------------------

# estimate theta using the binary undercount variable
# estimate lambda using mixture model

setwd("~/Dropbox/clarkiaSeedBanks/")
sink("test.R")
cat(
  "model{
  # priors
# mu ~ dnorm(0, 0.001)
# sigma ~ dunif(0,1.5)
# tau <- 1/(sigma*sigma)
theta ~ dbeta(1,1)
# p ~ dbeta(1, 1)
lambda ~ dgamma(0.01, 0.01)
lambda2 ~ dgamma(0.01, 0.01)

# likelihood
for(i in 1:n){

 # undercount[i] ~ dbern(p)
  seedlingNumber[i] ~ dpois(lambda)
  fruitplNumber[i] ~ dpois(lambda2)

  # u[i] =  seedlingNumber[i]# +  n_unobs[i]

 # alpha ~ dnorm(mu, tau)

  # logit(theta) <- alpha
  # fruitplNumber[i] ~ dbinom(theta, seedlingNumber[i])
  
}

  p = lambda2/lambda

}
")
sink()
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 3000
n.update = 5000
n.iterations = 1000
n.thin = 1

set.seed(10)
dir = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Complete pooling of germination and viability trials
# Partial pooling of seed burial experiment (site level)
# -------------------------------------------------------------------
# -------------------------------------------------------------------
#data = data %>% dplyr::filter(undercount==0)
data <- tidybayes::compose_data(data) 

library(rjags)
# # tuning (n.adapt)
jm = jags.model(paste0(dir,"test.R"), data = data, n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

parsToMonitor_1 = c("theta","lambda","lambda2","p")


# chain (n.iter)
samples.rjags = coda.samples(jm, 
                             variable.names = c(parsToMonitor_1), 
                             n.iter = n.iterations, thin = n.thin)
MCMCsummary(samples.rjags, n.eff = TRUE, params = parsToMonitor_1)

hist(data$fruitplNumber/data$seedlingNumber,breaks=15)
abline(v=quantile(MCMCchains(samples.rjags,params="theta"),.5),col='red')
sum(data$fruitplNumber)/sum(data$seedlingNumber)
