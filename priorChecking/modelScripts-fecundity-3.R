# list of JAGS scripts used in this file

# -------------------------------------------------------------------
# Models for 
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
library(magrittr)
library(tidybayes)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
countFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countFruitsPerPlantAllPlots.RDS")

countFruitsPerPlantAllPlots <- countFruitsPerPlantAllPlots %>%
  dplyr::rename(y_tfe = countFruitNumberPerPlant) %>%
  dplyr::select(site,year,y_tfe) 

countFruitsPerPlantAllPlots$year <- as.character(countFruitsPerPlantAllPlots$year)

countUndamagedDamagedFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countUndamagedDamagedFruitsPerPlantAllPlots.RDS")

countUndamagedDamagedFruitsPerPlantAllPlots <- countUndamagedDamagedFruitsPerPlantAllPlots %>%
  dplyr::rename(y_und = countUndamagedFruitNumberPerPlant) %>%
  dplyr::rename(y_dam = countDamagedFruitNumberPerPlant) %>%
  dplyr::rename(site2 = site) %>%
  dplyr::rename(year2 = year) %>%
  dplyr::select(site2,year2,y_und,y_dam) 

countUndamagedDamagedFruitsPerPlantAllPlots$year2 <- as.character(countUndamagedDamagedFruitsPerPlantAllPlots$year2)

countSeedPerFruit <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countSeedPerFruit.RDS")

countSeedPerUndamagedFruit <- countSeedPerFruit %>%
  dplyr::filter(demography==1) %>%
  dplyr::filter(damaged==0) %>%
  dplyr::rename(site3 = site) %>%
  dplyr::rename(year3 = year) %>%
  dplyr::select(site3,year3,sdno)

countSeedPerUndamagedFruit$year3 <- as.character(countSeedPerUndamagedFruit$year3)

countSeedPerDamagedFruit <- countSeedPerFruit %>%
  dplyr::filter(demography==1) %>%
  dplyr::filter(damaged==1) %>%
  dplyr::rename(site4 = site) %>%
  dplyr::rename(year4 = year) %>%
  dplyr::rename(sdno_dam = sdno) %>%
  dplyr::select(site4,year4,sdno_dam)

countSeedPerDamagedFruit$year4 <- as.character(countSeedPerDamagedFruit$year4)


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------
# -------------------------------------------------------------------
data <- tidybayes::compose_data(countFruitsPerPlantAllPlots,
                                countUndamagedDamagedFruitsPerPlantAllPlots,
                                countSeedPerUndamagedFruit,
                                countSeedPerDamagedFruit)

#data$n = dim(countFruitsPerPlantAllPlots)[1]
data$n = 1

data$n2 = dim(countUndamagedDamagedFruitsPerPlantAllPlots)[1]
data$n3 = dim(countSeedPerUndamagedFruit)[1]
data$n4 = dim(countSeedPerDamagedFruit)[1]
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Function to create indexes
# -------------------------------------------------------------------
# -------------------------------------------------------------------
f<-function(dataset=countUndamagedDamagedFruitsPerPlantAllPlots,response="y_und"){
  
  columnNames = colnames(dataset)
  siteCol<-grep("site",columnNames)
  yearCol<-grep("year",columnNames)
  
  df <-dataset %>% 
    dplyr::rename(site = columnNames[siteCol]) %>% 
    dplyr::rename(year = columnNames[yearCol])
  
  siteLen<-length(unique(df$site))
  yearLen<-length(unique(df$year))
  
  sites=data.frame(site=unique(df$site),siteIndex=1:siteLen)
  years=data.frame(year=unique(df$year),yearIndex=1:yearLen)
  
  
  d<-df %>%
    dplyr::left_join(sites,by="site") %>%
    dplyr::left_join(years,by="year") %>%
    # see https://stackoverflow.com/questions/48219732/pass-a-string-as-variable-name-in-dplyrfilter
    dplyr::filter(!is.na(get(response))) %>%
    dplyr::select(site,year,siteIndex,yearIndex) %>%
    unique %>%
    dplyr::select(siteIndex,yearIndex) %>%
    arrange(siteIndex) %>%
    tidyr::pivot_wider(values_from=yearIndex,names_from=siteIndex)
  
  listIndices<-lapply(d, `[[`, 1)
  
  sparse<-as.vector(unlist(listIndices))
  
  count=df %>%
    dplyr::left_join(sites,by="site") %>%
    dplyr::left_join(years,by="year") %>%
    dplyr::filter(!is.na(get(response))) %>%
    dplyr::group_by(site,year) %>%
    dplyr::summarise(count=n())
  
  counts=count$count
  item=as.vector(rep(sparse,counts))
  
  ticker=df %>%
    dplyr::left_join(sites,by="site") %>%
    dplyr::left_join(years,by="year") %>%
    dplyr::filter(!is.na(get(response))) %>%
    dplyr::mutate(number=1) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(site,year) %>%
    dplyr::summarise(ticker=cumsum(number))
  
  id = ticker[ticker$ticker==1,]
  id = data.frame(site=id$site,id=1:length(id$site)) %>%
    group_by(site) %>%
    dplyr::mutate(
      first = dplyr::first(id)
    )
  
  id=c(unique(id$first),length(sparse)+1)
  
  list(sparse=sparse,item=item,id=id)
  
}

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Set indexes
# -------------------------------------------------------------------
# -------------------------------------------------------------------
index_tfe <- f(dataset=countFruitsPerPlantAllPlots,response="y_tfe")
data$sparse_tfe = index_tfe$sparse
data$item_tfe = index_tfe$item
data$id_tfe = index_tfe$id

# index_und <- f(dataset=countUndamagedDamagedFruitsPerPlantAllPlots,response="y_und")
# data$sparse_und = index_und$sparse
# data$item_und = index_und$item
# data$id_und = index_und$id
# 
# index_seeds <- f(dataset=countSeedPerUndamagedFruit,response="sdno")
# data$sparse_seeds = index_seeds$sparse
# data$item_seeds = index_seeds$item
# data$id_seeds = index_seeds$id
# 
# index_dam_seeds <- f(dataset=countSeedPerDamagedFruit,response="sdno_dam")
# data$sparse_dam_seeds = index_dam_seeds$sparse
# data$item_dam_seeds = index_dam_seeds$item
# data$id_dam_seeds = index_dam_seeds$id
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

# # -------------------------------------------------------------------
# # -------------------------------------------------------------------
# # Complete pooling of germination and viability trials
# # Partial pooling of seed burial experiment (site level)
# # -------------------------------------------------------------------
# # -------------------------------------------------------------------
# initsMu0 <- function(samps = data$n_site2){
#   #rnorm(n = samps, mean = 0, sd = 1)
#   #truncdist::rtrunc(n = samps,"gamma",a=.001,shape=0.01,rate=0.01)
#   rgamma(n = samps, shape = 1, rate = 1)
# }
# 
# initsSigma0 <- function(samps = data$n_site2){
#   extraDistr::rhnorm(n = samps, sigma = 1)
#   #truncdist::rtrunc(n = samps,"gamma",a=.001,shape=0.001,rate=0.001)
#   # rgamma(n = samps, shape = .001, rate = .001)
#   
# }
# 
# initsR <- function(rows = data$n_site2, cols = data$n_year2){
#   samps = rows*cols
#   #matrix(truncdist::rtrunc(n = samps,"gamma",a=.001,shape=0.001,rate=0.001), rows, cols)
#   #matrix(runif(n = samps,0,1), rows, cols)
#  # matrix(truncdist::rtrunc(n = samps,"gamma",a=.001,shape=0.001,rate=0.001), rows, cols)
#   matrix(extraDistr::rhnorm(n = samps, sigma = 1), rows, cols)
# }
# 
# # set inits for JAGS
# inits <- list()
# for(i in 1:3){
#   inits[[i]] <- list(initsMu0(), initsMu0(), initsMu0(), initsMu0(), initsMu0(),
#                      initsSigma0(), initsSigma0(), initsSigma0(), initsSigma0(), initsSigma0(),
#                      initsR(rows = data$n_site2,cols=data$n_year), 
#                      initsR(rows = data$n_site2,cols=data$n_year2), 
#                      initsR(rows = data$n_site2,cols=data$n_year2), 
#                      initsR(rows = data$n_site2,cols=data$n_year3), 
#                      initsR(rows = data$n_site2,cols=data$n_year4) )
#   
#   # names(inits[[i]]) = c(paste(rep("nu",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
#   #                       paste(rep("tau0",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
#   #                       paste(rep("tau",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"))
#   
#   names(inits[[i]]) = c(paste(rep("nu",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
#                         paste(rep("sigma0",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"),
#                         paste(rep("sigma",5),c("tfe","und","dam","seeds","dam_seeds"),sep="_"))
#   
# }

# # Call to JAGS
# 
# # tuning (n.adapt)
jm = jags.model(paste0(dir,"nbLikelihood-logLink-normalHierarchical-3.R"), 
                data = data, #inits = inits,
                n.chains = n.chain, n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = n.update)

parsToMonitor = c("mu0_tfe","sigma0_tfe","tau0_tfe",
                  "mu_tfe","sigma_tfe",
                  "kappa_tfe",
                  "lambda_tfe",
                  "y_prior")

# chain (n.iter)
samples.rjags = coda.samples(jm, 
                             variable.names = c(parsToMonitor), 
                             n.iter = n.iterations, thin = n.thin)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Evaluate priors
# -------------------------------------------------------------------
# -------------------------------------------------------------------

hist.chain=function(x){
  name<-colnames(x)
  tmp = hist(x,breaks=100,plot=FALSE)
  tmp$xname = name
  plot(tmp,freq=FALSE)
}

## Priors on hyperparameters
par(mfrow=c(2,2))
hist.chain(MCMCchains(samples.rjags,params="mu0_tfe")[,1,drop=FALSE])
#hist.chain(MCMCchains(samples.rjags,params="tau0_tfe")[,1,drop=FALSE])
hist.chain(MCMCchains(samples.rjags,params="sigma0_tfe")[,1,drop=FALSE])


#par(mfrow=c(2,2))
hist.chain(MCMCchains(samples.rjags,params="mu_tfe")[,1,drop=FALSE])
hist.chain(MCMCchains(samples.rjags,params="sigma_tfe")[,1,drop=FALSE])

# hist.chain(MCMCchains(samples.rjags,params="kappa_tfe")[,1,drop=FALSE])
# 
# hist.chain(exp(MCMCchains(samples.rjags,params="mu_tfe")[,1,drop=FALSE]))
par(mfrow=c(1,2))
hist.chain(MCMCchains(samples.rjags,params="lambda_tfe")[,1,drop=FALSE])
bounds=rethinking::HPDI(MCMCchains(samples.rjags,params="lambda_tfe"),prob=.95)
segments(x0=bounds[1],y0=0,x1=bounds[2],lwd=3,col='red')

hist.chain(MCMCchains(samples.rjags,params="y_prior")[,1,drop=FALSE])
bounds=rethinking::HPDI(MCMCchains(samples.rjags,params="y_prior"),prob=.95)
segments(x0=bounds[1],y0=0,x1=bounds[2],lwd=3,col='red')



# par(mfrow=c(1,1))

# 
# plot(y_prior[order(y_prior)][1:15000])
# plot(y_prior[order(y_prior)][15001:29000])

x=seq(min(MCMCchains(samples.rjags,params="mu0_tfe")[,1,drop=FALSE]),
      max(MCMCchains(samples.rjags,params="mu0_tfe")[,1,drop=FALSE]),
      by=0.1)
plot(x,exp(x),type='l')

y_prior = MCMCchains(samples.rjags,params="y_prior")
plot((y_prior))
quantile(y_prior,c(.001,.999))
rethinking::HPDI(y_prior,prob=.999)
rethinking::HPDI(MCMCchains(samples.rjags,params="lambda_tfe"),prob=.999)
quantile(y_prior,c(.025,.975))

# 
# hist(rlnorm(1000,0,1),breaks=100)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Evaluate data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# f<-function(x){
#   plot((1:length(x))/length(x),log(x[order(x)]+.01),type='l',ylim=c(-2,8))
# }
# f2<-function(x){
#   lines((1:length(x))/length(x),log(x[order(x)]+.01),type='l')
# }
# par(mfrow=c(1,1))
# f(countFruitsPerPlantAllPlots$y_tfe)
# f2(countUndamagedDamagedFruitsPerPlantAllPlots$y_und)
# f2(countUndamagedDamagedFruitsPerPlantAllPlots$y_dam)
# f2(countSeedPerFruit$sdno)
# f2(countSeedPerUndamagedFruit$sdno)
# f2(countSeedPerDamagedFruit$sdno_dam)
# 
quantile(countFruitsPerPlantAllPlots$y_tfe,probs=c(0.01,0.5,0.999),na.rm=TRUE)
quantile(countUndamagedDamagedFruitsPerPlantAllPlots$y_und,probs=c(0.01,0.5,0.999),na.rm=TRUE)
quantile(countUndamagedDamagedFruitsPerPlantAllPlots$y_dam,probs=c(0.01,0.5,0.999),na.rm=TRUE)
quantile(countSeedPerFruit$sdno,probs=c(0.01,0.5,0.999),na.rm=TRUE)
quantile(countSeedPerUndamagedFruit$sdno,probs=c(0.01,0.5,0.999),na.rm=TRUE)
quantile(countSeedPerDamagedFruit$sdno_dam,probs=c(0.01,0.5,0.999),na.rm=TRUE)

quantile(countFruitsPerPlantAllPlots$y_tfe,probs=c(0.025,0.5,.975),na.rm=TRUE)
quantile(countUndamagedDamagedFruitsPerPlantAllPlots$y_und,probs=c(0.025,0.5,.975),na.rm=TRUE)
quantile(countUndamagedDamagedFruitsPerPlantAllPlots$y_dam,probs=c(0.025,0.5,.975),na.rm=TRUE)
quantile(countSeedPerFruit$sdno,probs=c(0.025,0.5,.975),na.rm=TRUE)
quantile(countSeedPerUndamagedFruit$sdno,probs=c(0.025,0.5,.975),na.rm=TRUE)
quantile(countSeedPerDamagedFruit$sdno_dam,probs=c(0.025,0.5,.975),na.rm=TRUE)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# plot prior
# https://bookdown.org/marklhc/notes_bookdown/one-parameter-models.html#poisson-data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# p1 <- ggplot(tibble(log_lambda = c(-10, 10)), aes(x = log_lambda)) + 
#   stat_function(fun = dnorm, args = list(mean = 0, sd = 1)) + 
#   labs(x = expression(log(lambda)), y = "Density")
# p2 <- ggplot(tibble(lambda = c(0, 20)), aes(x = lambda)) + 
#   stat_function(fun = function(x) dnorm(log(x), mean = 0, sd = 1), 
#                 n = 501) + 
#   labs(x = expression(lambda), y = "Density")
# gridExtra::grid.arrange(p1, p2, nrow = 1)
# 
# x=seq(-4,5,
#       by=0.1)
# plot(x,exp(x),type='l')
# x=c(0,4)
# c(exp(x[1]-x[2]),exp(x[1]),exp(x[1]+x[2]))
# 
# x=c(0,1.3,3.5)
# exp(x)
# 
# a=c(0,1.3,3.5)
# b=c(1,1,1)
# exp(a+(b^2)/2)
# sqrt((exp(b^2)-1)*exp(2*a+b^2))

