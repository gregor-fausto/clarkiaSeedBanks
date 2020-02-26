df<-readxl::read_excel(path="~/Dropbox/dataLibrary/Clarkia-LTREB/20_demography_sites/Survivorship & Fecundity_06-15.xls")
library(janitor)
library(tidyverse)
df<-df %>%
  janitor::clean_names(case="lower_camel")

vars <- names(df)
seedlingNo<-tidyselect::vars_select(vars,contains("seedlingNumber"))
frtPlantNo<-tidyselect::vars_select(vars,contains("fruitplNumber6"))

df<-df %>%
  dplyr::select(site,transect,position,seedlingNo,frtPlantNo)

seedlingData<-pivot_longer(data = df %>% dplyr::select(site,transect,position,seedlingNo), 
                           cols = seedlingNo,
                           names_to = "year",
                           values_to = "seedlingNumber") %>%
  tidyr::separate(year, into = c("discard", "year")) %>%
  dplyr::mutate(year=as.numeric(paste0(20,year))) %>%
  dplyr::select(-discard)
frtData<-pivot_longer(data = df %>% dplyr::select(site,transect,position,frtPlantNo), 
                      cols = frtPlantNo,
                      names_to = "year",
                      values_to = "fruitingPlantNumber") %>%
  tidyr::separate(year, into = c("discard", "year")) %>%
  dplyr::mutate(year=as.numeric(paste0(20,year))) %>%
  dplyr::select(-discard)

df<-seedlingData %>%
  dplyr::left_join(frtData,by=c("site","transect","position","year")) %>%
  dplyr::rename(plot=position)

library(janitor)

dfAnalysis<-df %>%
    dplyr::filter(seedlingNumber>0 & !is.na(fruitingPlantNumber))

out<-dfAnalysis %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(x1 = sum(seedlingNumber), x2 = sum(fruitingPlantNumber)) %>%
  dplyr::mutate(sigma = x2/x1)


dfAnalysis<-dfAnalysis %>%
  dplyr::mutate(seedlingNumber = ifelse(seedlingNumber<fruitingPlantNumber,fruitingPlantNumber,seedlingNumber))
# binomNLL1 = function(p, y, n){
#   -sum(dbinom(y, prob = p, size = n, log = TRUE))
# }
# mle2( minuslogl = binomNLL1,start=list(p=.5),data=list(n=n,y=y))

d<-split(dfAnalysis,list(dfAnalysis$site,dfAnalysis$year))
d.list <- list(d[1:20],d[21:40],d[41:60],d[61:80],d[81:100],d[101:120],d[121:140],d[141:160],d[161:180],d[181:200])
optima<-matrix(NA,20,10)
samps<-matrix(NA,20,10)
for(i in 1:10){
  tmp.list <- d.list[[i]]
  for(k in 1:20){
    tmp<-tmp.list[[k]]
    #opt<-optim(0.5, negLogLik, lower = 0, upper = 1, y = tmp$fruitingPlantNumber, n = tmp$seedlingNumber, method = "Brent")
    y = tmp$fruitingPlantNumber
    n = tmp$seedlingNumber
    negLogLik <- function(p,x=tmp$fruitingPlantNumber,size=tmp$seedlingNumber) -sum(dbinom(x=x, size=size, prob=p,log=TRUE))
    opt<-optimise(negLogLik,interval=c(0,1),maximum=FALSE)
    optima[k,i] <- opt$minimum
    samps[k,i] <- sum(tmp$seedlingNumber)
  }
}
MLEs<-data.frame(site=unique(dfAnalysis$site),optima )
#MLEs$mean = round(rowMeans(MLEs[,2:5]),3)
names(MLEs)=c("site",2006:2015)
out<-pivot_longer(data=MLEs,cols=`2006`:`2015`,names_to="year",values_to="p")
out$year<-as.numeric(out$year)

summaryMLE <- out
# ggplot(data.frame(out),aes(x=year,y=p)) +
#   geom_point() +
#   stat_summary(fun.y = "mean",col="red",geom="point")
# 
# ggplot(data.frame(out),aes(x=site,y=p)) +
#   geom_point() +
#   stat_summary(fun.y = "mean",col="red",geom="point")
# 
# ggplot() +
#   geom_point(data=dfAnalysis,aes(x=site,y=fruitingPlantNumber/seedlingNumber),cex=0.25) +
#   geom_point(data=data.frame(out),aes(x=site,y=p),color="red") +
#   facet_wrap(~year)


dfSite <- read.csv(file="/Users/Gregor/Dropbox/dataLibrary/Datafile_Site_Environment_corr.csv",header=TRUE)

dfSite<-dfSite %>%
  janitor::clean_names(case="lower_camel") %>%
  dplyr::select(site,easting) %>% unique %>%
  dplyr::filter(site %in% unique(df$site))

# ggplot() +
#   geom_point(data=dfAnalysis%>% dplyr::left_join(dfSite,by="site"),aes(x=easting,y=fruitingPlantNumber/seedlingNumber),cex=0.25) +
#   geom_point(data=data.frame(out)%>% dplyr::left_join(dfSite,by="site"),aes(x=easting,y=p),color="red") +
#   geom_smooth(data=dfAnalysis%>% dplyr::left_join(dfSite,by="site"),aes(x=easting,y=fruitingPlantNumber/seedlingNumber),method="lm") +
#   geom_smooth(data=data.frame(out)%>% dplyr::left_join(dfSite,by="site"),aes(x=easting,y=p),color="red",method="lm") +
#     facet_wrap(~year) +
#   theme_minimal()

library(tidybayes)

dfBayes<-dfAnalysis %>% dplyr::select(-c(transect,plot)) %>%
  dplyr::mutate(year = as.factor(year))
# pass data to list for JAGS
# data = list(
#   y = dfAnalysis$fruitingPlantNumber,
#   n = dfAnalysis$seedlingNumber,
#   N = dim(dfAnalysis)[1],
#   pop = as.numeric(as.factor(dfAnalysis$site)),
#   npop = length(unique(dfAnalysis$site)),
#   year = as.numeric(as.factor(dfAnalysis$year)),
#   nyear = length(unique(dfAnalysis$year))
# )
data<-tidybayes::compose_data(dfBayes)

# -------------------------------------------------------------------
# Code for JAGS model
# -------------------------------------------------------------------

modelString1 <- "model { 
  ## Priors
  
  # hyperpriors
  for(j in 1:n_site){
    for(k in 1:n_year){
      p[j,k] ~ dbeta(1,1)
    }
  }
  
  # Likelihood
  for(i in 1:n){
    fruitingPlantNumber[i] ~ dbin(p[site[i],year[i]], seedlingNumber[i])
  }
  
}"

# -------------------------------------------------------------------
# Initial values
# -------------------------------------------------------------------
# set inits for rjags
inits = list(list(p = matrix(rep(.1,data$n_site*data$n_year),nrow=data$n_site,ncol=data$n_year)),
             list(p = matrix(rep(.5,data$n_site*data$n_year),nrow=data$n_site,ncol=data$n_year)),
             list(p = matrix(rep(.9,data$n_site*data$n_year),nrow=data$n_site,ncol=data$n_year))) 

# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# -------------------------------------------------------------------
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 500
n.update = 5000
n.iterations = 1000
n.thin = 1
parsToMonitor = c("p")
library(rjags)

set.seed(2)
# tuning (n.adapt)
jm1 = jags.model(textConnection(modelString1), data = data, inits = inits,
                 n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm1, n.iterations = n.update)

# chain (n.iter)
samples.rjags1 = coda.samples(jm1, variable.names = c(parsToMonitor), n.iter = n.iterations, thin = n.thin)

MCMCvis::MCMCsummary(samples.rjags1)


samples.rjags1 %<>% recover_types(dfBayes)

summaryBayes<-samples.rjags1 %>% 
  tidybayes::recover_types(dfBayes) %>% 
  tidybayes::spread_draws(p[site,year]) %>%
  dplyr::summarise(bayes1=median(p)) 

summaryBayes$year <- as.character(summaryBayes$year)
summaryBayes$site <- as.factor(summaryBayes$site)

MCMCvis::MCMCchains(samples.rjags1) %>%
  tidybayes::spread_draws(p[pop,year]) %>%
  head(15) 
  dplyr::group_by(pop,year) %>%
  dplyr::summarise(bayes1=median(p)) %>%
  dplyr


### COMPARISON
  
ns<-  dfBayes %>%
    dplyr::group_by(site,year) %>%
    dplyr::summarise(n())

ns$site <- as.factor(ns$site)
ns$year <- as.numeric(as.character(ns$year))
summaryBayes$year <- as.numeric(as.character(summaryBayes$year))

  
comp<-  summaryMLE %>%
    dplyr::left_join(summaryBayes,by=c("site","year")) %>%
  dplyr::left_join(ns,by=c("site","year"))

comp %>% dplyr::mutate(abs(p-bayes1)) %>% View

# plot comparing MLE fit and bayes fit 

setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
pdf(file="appendix-x-mle_bayes.pdf", width=8, height=4)
par(mfrow=c(1,2))
par(mar=c(5, 6, 4, 2) + 0.1)
plot(comp$p,comp$bayes1,
     xlab="Maximum likelihood estimate",
     ylab="Median of posterior for \n beta-binomial parameterization",
     xlim=c(0,1),ylim=c(0,1),
     pch=16,cex=.5)  
abline(a=0,b=1,lwd=0.5)



d<-MCMCvis::MCMCchains(samples.rjags1,params="p")
sumMLE<-comp %>%
  dplyr::arrange(year)
tmp<-matrix(NA,ncol=200,nrow=3000)
for(i in 1:200){
  tmp[,i]<-d[,i]-sumMLE$p[i]
}
delta<-apply(tmp,2,median)
delta.conf <- apply(tmp,2,quantile,probs=c(.025,.5,.975))
par(mar=c(5, 6, 4, 2) + 0.1)
plot(sumMLE$`n()`,t(delta.conf)[,2],
     pch=16,cex=0.5,
     ylim=c(-1,1),
     xlab="Sample size (n)",
     ylab="Comparison of the maximum likelihood estimate\n and the beta-binomial parameterization\n (median and 95% CI)")
segments(x0=sumMLE$`n()`,y0=t(delta.conf)[,1],
         x1=sumMLE$`n()`,y1=t(delta.conf)[,3])
abline(h=0,lwd=1)

dev.off()


library(tidybayes)

dfBayes<-dfAnalysis %>% 
  dplyr::filter(site=="BG") %>%
  dplyr::select(-c(transect,plot)) %>%
  dplyr::mutate(year = as.factor(year))
# pass data to list for JAGS
# data = list(
#   y = dfAnalysis$fruitingPlantNumber,
#   n = dfAnalysis$seedlingNumber,
#   N = dim(dfAnalysis)[1],
#   pop = as.numeric(as.factor(dfAnalysis$site)),
#   npop = length(unique(dfAnalysis$site)),
#   year = as.numeric(as.factor(dfAnalysis$year)),
#   nyear = length(unique(dfAnalysis$year))
# )
data<-tidybayes::compose_data(dfBayes)

# --------------------------------------------
# -------------------------------------------------------------------
# Code for JAGS model
# -------------------------------------------------------------------

modelString2 <- "
  model {

    omega ~ dbeta( 1 , 1 ) # broad uniform
        kappaMinusTwo ~ dgamma( .01 , .01 )  # mode=1 , sd=10 
    kappa <- kappaMinusTwo + 2
    
        for ( i in 1:n_year ) {
      theta[i] ~ dbeta( omega*(kappa-2)+1 , (1-omega)*(kappa-2)+1 ) 
    }
    
        for ( i in 1:n ) {
      fruitingPlantNumber[i] ~ dbin( theta[year[i]], seedlingNumber[i] )
    }
  }
"

# -------------------------------------------------------------------
# Initial values
# -------------------------------------------------------------------
# set inits for rjags
# INTIALIZE THE CHAINS.
# Initial values of MCMC chains based on data:
initsList = function() {
  thetaInit = rep(0,data$n_site)
  for ( sIdx in 1:data$n_site ) { # for each subject
    includeRows = ( data$site == sIdx ) # identify rows of this subject
    yThisSubj = data$fruitingPlantNumber[includeRows]  # extract data of this subject
    resampledY = sample( yThisSubj , replace=TRUE ) # resample
    thetaInit[sIdx] = sum(resampledY)/length(resampledY) 
  }
  thetaInit = 0.001+0.998*thetaInit # keep away from 0,1
  meanThetaInit = mean( thetaInit )
  kappaInit = 100 # lazy, start high and let burn-in find better value
  return( list( list( theta = rep(.1, 10), omega=.5 , 
                kappaMinusTwo=kappaInit-2 ),
                list( theta = rep(.5, 10),  omega=.5 , 
                     kappaMinusTwo=kappaInit-2 ),
                list( theta = rep(.9, 10), omega=.5 , 
                     kappaMinusTwo=kappaInit-2 )))
}

inits = initsList()

# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# -------------------------------------------------------------------
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 500
n.update = 5000
n.iterations = 100000
n.thin = 1
parsToMonitor = c("theta","omega","kappa")
library(rjags)

set.seed(2)
# tuning (n.adapt)
jm2 = jags.model(textConnection(modelString2), data = data,
                 n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm2, n.iterations = n.update)

# chain (n.iter)
samples.rjags2 = coda.samples(jm2, variable.names = c(parsToMonitor), 
                              n.iter = n.iterations, thin = n.thin)

MCMCvis::MCMCsummary(samples.rjags2)

library(tidybayes)

samples.rjags2 %<>% recover_types(dfBayes)

summaryBayes<-samples.rjags2 %>% 
  tidybayes::recover_types(dfBayes) %>% 
  tidybayes::spread_draws(theta[site]) %>%
  dplyr::summarise(bayes2=median(theta)) 

summaryBayes$year <- as.character(summaryBayes$year)
summaryBayes$site <- as.factor(summaryBayes$site)


d<-dfBayes %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(y =sum(fruitingPlantNumber), n = sum(seedlingNumber)) %>%
  dplyr::mutate(prop = y/n)

setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
pdf(file="appendix-x-hierarchical3.pdf", width=8, height=6)

par(mfrow=c(2,5))
for(i in 1:10){
hist(MCMCvis::MCMCchains(samples.rjags2,params=c("theta"))[,i],breaks=40,
     main=c(2006:2015)[i],
     xlab=expression("theta"),
     freq=FALSE)
  abline(v=d$prop[i],col='red',lty='dashed')
  lines(density(MCMCvis::MCMCchains(samples.rjags2,params=c("omega"))),col="red")
  
}

dev.off()


# model 3
library(tidybayes)

dfBayes<-dfAnalysis %>% 
  dplyr::select(-c(transect,plot)) %>%
  dplyr::mutate(year = as.factor(year))
# pass data to list for JAGS
# data = list(
#   y = dfAnalysis$fruitingPlantNumber,
#   n = dfAnalysis$seedlingNumber,
#   N = dim(dfAnalysis)[1],
#   pop = as.numeric(as.factor(dfAnalysis$site)),
#   npop = length(unique(dfAnalysis$site)),
#   year = as.numeric(as.factor(dfAnalysis$year)),
#   nyear = length(unique(dfAnalysis$year))
# )
data<-tidybayes::compose_data(dfBayes)

# --------------------------------------------
# -------------------------------------------------------------------
# Code for JAGS model
# -------------------------------------------------------------------

modelString3 <- "
  model {

for(j in 1:n_site){
      omega[j] ~ dbeta( 1 , 1 ) # broad uniform
      kappaMinusTwo[j]  ~ dgamma( .01 , .01 )  # mode=1 , sd=10 
      kappa[j]  <- kappaMinusTwo[j]  + 2
    }
    
    for(j in 1:n_site){   
        for ( i in 1:n_year ) {
      theta[j,i] ~ dbeta( omega[j]*(kappa[j]-2)+1 , (1-omega[j])*(kappa[j]-2)+1 ) 
        }
    }
    
        for ( i in 1:n ) {
      fruitingPlantNumber[i] ~ dbin( theta[site[i],year[i]], seedlingNumber[i] )
    }
  }
"

# -------------------------------------------------------------------
# Initial values
# -------------------------------------------------------------------
kappaInit = 100 # lazy, start high and let burn-in find better value

inits<-list( list( theta = matrix(rep(.1,data$n_site*data$n_year),nrow=data$n_site,ncol=data$n_year),
            omega=rep(.5,data$n_site) ,
            kappaMinusTwo=rep(98,data$n_site) ),
      list( theta = matrix(rep(.5,data$n_site*data$n_year),nrow=data$n_site,ncol=data$n_year),
            omega=rep(.5,data$n_site) ,
            kappaMinusTwo=rep(98,data$n_site) ),
      list( theta = matrix(rep(.9,data$n_site*data$n_year),nrow=data$n_site,ncol=data$n_year),
            omega=rep(.5,data$n_site) ,
            kappaMinusTwo=rep(98,data$n_site) ))


# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# -------------------------------------------------------------------
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 500
n.update = 5000
n.iterations = 10000
n.thin = 1
parsToMonitor = c("theta","omega","kappa")
library(rjags)

set.seed(2)
# tuning (n.adapt)
jm3 = jags.model(textConnection(modelString3), data = data, inits=inits,
                 n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm3, n.iterations = n.update)

# chain (n.iter)
samples.rjags3 = coda.samples(jm3, variable.names = c(parsToMonitor), 
                              n.iter = n.iterations, thin = n.thin)

MCMCvis::MCMCsummary(samples.rjags3)

library(tidybayes)

samples.rjags3 %<>% tidybayes::recover_types(dfBayes)

summaryBayes3<-samples.rjags3 %>% 
  tidybayes::recover_types(dfBayes) %>% 
  tidybayes::spread_draws(theta[site,year]) %>%
  dplyr::summarise(bayes2=median(theta)) 

summaryBayes3$year <- as.character(summaryBayes3$year)
summaryBayes3$site <- as.factor(summaryBayes3$site)


ns<-  dfBayes %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n())

ns$site <- as.factor(ns$site)
ns$year <- as.numeric(as.character(ns$year))
summaryBayes3$year <- as.numeric(as.character(summaryBayes3$year))

comp<-  summaryMLE %>%
  dplyr::left_join(summaryBayes3,by=c("site","year")) %>%
  dplyr::left_join(ns,by=c("site","year"))

# plot comparing MLE fit and bayes fit 

setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
pdf(file="appendix-x-mle_bayeshier.pdf", width=8, height=4)
par(mfrow=c(1,2))
par(mar=c(5, 6, 4, 2) + 0.1)
plot(comp$p,comp$bayes2,
     xlab="Maximum likelihood estimate",
     ylab="Median of posterior for \n beta-binomial parameterization",
     xlim=c(0,1),ylim=c(0,1),
     pch=16,cex=.5)  
abline(a=0,b=1,lwd=0.5)

d<-MCMCvis::MCMCchains(samples.rjags3,params="theta")
sumMLE<-comp %>%
  dplyr::arrange(year)
tmp<-matrix(NA,ncol=200,nrow=30000)
for(i in 1:200){
  tmp[,i]<-d[,i]-sumMLE$p[i]
}
delta<-apply(tmp,2,median)
delta.conf <- apply(tmp,2,quantile,probs=c(.025,.5,.975))
par(mar=c(5, 6, 4, 2) + 0.1)
plot(sumMLE$`n()`,t(delta.conf)[,2],
     pch=16,cex=0.5,
     ylim=c(-1,1),
     xlab="Sample size (n)",
     ylab="Comparison of the maximum likelihood estimate\n and the beta-binomial parameterization\n (median and 95% CI)")
segments(x0=sumMLE$`n()`,y0=t(delta.conf)[,1],
         x1=sumMLE$`n()`,y1=t(delta.conf)[,3])
abline(h=0,lwd=1)

dev.off()


summaryBayes<-samples.rjags1 %>% 
  tidybayes::recover_types(dfBayes) %>% 
  tidybayes::spread_draws(p[site,year]) %>%
  dplyr::summarise(bayes1=median(p)) 

summaryBayes$year <- as.character(summaryBayes$year)
summaryBayes$site <- as.factor(summaryBayes$site)

summaryBayes3<-samples.rjags3 %>% 
  tidybayes::recover_types(dfBayes) %>% 
  tidybayes::spread_draws(theta[site,year]) %>%
  dplyr::summarise(bayes3=median(theta)) 

summaryBayes3$year <- as.character(summaryBayes3$year)
summaryBayes3$site <- as.factor(summaryBayes3$site)

sum <- summaryBayes %>%
  #dplyr::left_join(summaryBayes2,by=c("site","year")) %>%
  dplyr::left_join(summaryBayes3,by=c("site","year"))

plot(sum$bayes1,sum$bayes3)
abline(a=0,b=1)

tmp<-matrix(NA,ncol=200,nrow=3000)
for(i in 1:200){

tmp[,i]<-MCMCvis::MCMCchains(samples.rjags1,params="p")[,i] - sample(MCMCvis::MCMCchains(samples.rjags3,params="theta")[,i],3000)

}

setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
pdf(file="appendix-x-bayescomp_bayeshier.pdf", width=4, height=4)
delta.conf <- apply(tmp,2,quantile,probs=c(.025,.5,.975))
par(mar=c(5, 6, 4, 2) + 0.1)
plot(sumMLE$`n()`,t(delta.conf)[,2],
     pch=16,cex=0.5,
     ylim=c(-1,1),
     xlab="Sample size (n)",
     ylab="Comparison of the posteriors for \n complete and partial pooling parameterizations \n complete-partial (median and 95% CI)")
segments(x0=sumMLE$`n()`,y0=t(delta.conf)[,1],
         x1=sumMLE$`n()`,y1=t(delta.conf)[,3])
abline(h=0,lwd=1)
dev.off()


setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
pdf(file="appendix-x-mismatch.pdf", width=8, height=6)
v<-(1:200)[sumMLE$`n()`==1][!is.na((1:200)[sumMLE$`n()`==1])]

omegaSite<-c(11,11,13,6,7,11,12,20)

par(mfrow=c(2,4))
for(i in 1:length(v)){
plot(density(  MCMCvis::MCMCchains(samples.rjags1,params="p")[,v[i]] ),
     xlim=c(0,1),ylim=c(0,4),
     main="") 
  abline(v=median( MCMCvis::MCMCchains(samples.rjags1,params="p")[,v[i]] ))
tmp <- sample(MCMCvis::MCMCchains(samples.rjags3,params="theta")[,v[i]],3000)
lines( density(tmp),lty='dotted')
abline(v=median( tmp),lty='dotted')

lines( density(sample(MCMCvis::MCMCchains(samples.rjags3,params="omega")[,omegaSite[i]],3000)),col='red')

}
dev.off()

pdf(file="appendix-x-match.pdf", width=8, height=6)
v<-(1:200)[sumMLE$`n()`==30][!is.na((1:200)[sumMLE$`n()`==30])]


par(mfrow=c(2,4))
for(i in 1:8){
  plot(density(  MCMCvis::MCMCchains(samples.rjags1,params="p")[,v[i]] ),
       
       main="") 
  abline(v=median( MCMCvis::MCMCchains(samples.rjags1,params="p")[,v[i]] ))
  tmp <- sample(MCMCvis::MCMCchains(samples.rjags3,params="theta")[,v[i]],3000)
  lines( density(tmp),lty='dotted')
  abline(v=median( tmp),lty='dotted')

}
dev.off()


# # model 4
# library(tidybayes)
# 
# dfBayes<-dfAnalysis %>% 
#   dplyr::select(-c(transect,plot)) %>%
#   dplyr::mutate(year = as.factor(year))
# 
# data<-tidybayes::compose_data(dfBayes)
# 
# # --------------------------------------------
# # -------------------------------------------------------------------
# # Code for JAGS model
# # -------------------------------------------------------------------
# 
# modelString4 <- "
#   model {
# 
#   omega0 ~ dbeta(1,1)
#   kappaMinusTwo0 ~ dgamma(0.01,0.01)
#   kappa0 <- kappaMinusTwo0 + 2
# 
#   for(j in 1:n_site){
#       kappaMinusTwo[j]  ~ dgamma( .01 , .01 )  # mode=1 , sd=10 
#       kappa[j]  <- kappaMinusTwo[j]  + 2
#       
#       omega[j] ~ dbeta( omega0*(kappa0-2) + 1 , (1-omega0)*(kappa0-2) +1 ) # broad uniform
#     }
#     
#     for(j in 1:n_site){   
#         for ( i in 1:n_year ) {
#       theta[j,i] ~ dbeta( omega[j]*(kappa[j]-2)+1 , (1-omega[j])*(kappa[j]-2)+1 ) 
#         }
#     }
#     
#         for ( i in 1:n ) {
#       fruitingPlantNumber[i] ~ dbin( theta[site[i],year[i]], seedlingNumber[i] )
#     }
#   }
# "
# 
# # -------------------------------------------------------------------
# # Initial values
# # -------------------------------------------------------------------
# kappaInit = 100 # lazy, start high and let burn-in find better value
# 
# inits<-list( list( theta = matrix(rep(.1,data$n_site*data$n_year),nrow=data$n_site,ncol=data$n_year),
#                    omega=rep(.5,data$n_site) ,
#                    kappaMinusTwo=rep(98,data$n_site),
#                    omega0 = .5, kappaMinusTwo0 = 98),
#              list( theta = matrix(rep(.5,data$n_site*data$n_year),nrow=data$n_site,ncol=data$n_year),
#                    omega=rep(.5,data$n_site) ,
#                    kappaMinusTwo=rep(98,data$n_site),
#                    omega0 = .5, kappaMinusTwo0 = 98 ),
#              list( theta = matrix(rep(.9,data$n_site*data$n_year),nrow=data$n_site,ncol=data$n_year),
#                    omega=rep(.5,data$n_site) ,
#                    kappaMinusTwo=rep(98,data$n_site),
#                    omega0 = .5, kappaMinusTwo0 = 98 ))
# 
# 
# # -------------------------------------------------------------------
# # Set JAGS parameters and random seed
# # -------------------------------------------------------------------
# # scalars that specify the 
# # number of iterations in the chain for adaptation
# # number of iterations for burn-in
# # number of samples in the final chain
# n.adapt = 500
# n.update = 5000
# n.iterations = 10000
# n.thin = 1
# parsToMonitor = c("theta","omega","kappa","omega0","kappa0")
# library(rjags)
# 
# set.seed(2)
# # tuning (n.adapt)
# jm4 = jags.model(textConnection(modelString4), data = data, inits=inits,
#                  n.chains = length(inits), n.adapt = n.adapt)
# 
# # burn-in (n.update)
# update(jm4, n.iterations = n.update)
# 
# # chain (n.iter)
# samples.rjags4 = coda.samples(jm4, variable.names = c(parsToMonitor), 
#                               n.iter = n.iterations, thin = n.thin)
# 
# MCMCvis::MCMCsummary(samples.rjags4)
# 
# library(tidybayes)
# 
# samples.rjags4 %<>% tidybayes::recover_types(dfBayes)
# 
# summaryBayes4<-samples.rjags4 %>% 
#   tidybayes::recover_types(dfBayes) %>% 
#   tidybayes::spread_draws(theta[site,year]) %>%
#   dplyr::summarise(bayes4=median(theta)) 
# 
# summaryBayes4$year <- as.character(summaryBayes4$year)
# summaryBayes4$site <- as.factor(summaryBayes4$site)
# 
# 
# ns<-  dfBayes %>%
#   dplyr::group_by(site,year) %>%
#   dplyr::summarise(n())
# 
# ns$site <- as.factor(ns$site)
# ns$year <- as.numeric(as.character(ns$year))
# summaryBayes4$year <- as.numeric(as.character(summaryBayes4$year))
# 
# comp<-  summaryMLE %>%
#   dplyr::left_join(summaryBayes4,by=c("site","year")) %>%
#   dplyr::left_join(ns,by=c("site","year"))
# 
# # plot comparing MLE fit and bayes fit 
# 
# setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
# pdf(file="appendix-x-mle_bayeshierfull.pdf", width=8, height=4)
# par(mfrow=c(1,2))
# par(mar=c(5, 6, 4, 2) + 0.1)
# plot(comp$p,comp$bayes4,
#      xlab="Maximum likelihood estimate",
#      ylab="Median of posterior for \n beta-binomial parameterization",
#      xlim=c(0,1),ylim=c(0,1),
#      pch=16,cex=.5)  
# abline(a=0,b=1,lwd=0.5)
# 
# d<-MCMCvis::MCMCchains(samples.rjags4,params="theta")
# sumMLE<-comp %>%
#   dplyr::arrange(year)
# tmp<-matrix(NA,ncol=200,nrow=30000)
# for(i in 1:200){
#   tmp[,i]<-d[,i]-sumMLE$p[i]
# }
# delta<-apply(tmp,2,median)
# delta.conf <- apply(tmp,2,quantile,probs=c(.025,.5,.975))
# par(mar=c(5, 6, 4, 2) + 0.1)
# plot(sumMLE$`n()`,t(delta.conf)[,2],
#      pch=16,cex=0.5,
#      ylim=c(-1,1),
#      xlab="Sample size (n)",
#      ylab="Comparison of the maximum likelihood estimate\n and the beta-binomial parameterization\n (median and 95% CI)")
# segments(x0=sumMLE$`n()`,y0=t(delta.conf)[,1],
#          x1=sumMLE$`n()`,y1=t(delta.conf)[,3])
# abline(h=0,lwd=1)
# 
# dev.off()
# 
# 
# summaryBayes<-samples.rjags1 %>% 
#   tidybayes::recover_types(dfBayes) %>% 
#   tidybayes::spread_draws(p[site,year]) %>%
#   dplyr::summarise(bayes1=median(p)) 
# 
# summaryBayes$year <- as.character(summaryBayes$year)
# summaryBayes$site <- as.factor(summaryBayes$site)
# 
# summaryBayes4<-samples.rjags4 %>% 
#   tidybayes::recover_types(dfBayes) %>% 
#   tidybayes::spread_draws(theta[site,year]) %>%
#   dplyr::summarise(bayes4=median(theta)) 
# 
# summaryBayes4$year <- as.character(summaryBayes4$year)
# summaryBayes4$site <- as.factor(summaryBayes4$site)
# 
# sum <- summaryBayes %>%
#   #dplyr::left_join(summaryBayes2,by=c("site","year")) %>%
#   dplyr::left_join(summaryBayes4,by=c("site","year"))
# 
# plot(sum$bayes1,sum$bayes4)
# abline(a=0,b=1)
# 
# tmp<-matrix(NA,ncol=200,nrow=3000)
# for(i in 1:200){
#   
#   tmp[,i]<-MCMCvis::MCMCchains(samples.rjags3,params="theta")[,i] - sample(MCMCvis::MCMCchains(samples.rjags3,params="theta")[,i],3000)
#   
# }
# 
# setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
# pdf(file="appendix-x-bayescomp_bayeshierfull.pdf", width=4, height=4)
# delta.conf <- apply(tmp,2,quantile,probs=c(.025,.5,.975))
# par(mar=c(5, 6, 4, 2) + 0.1)
# plot(sumMLE$`n()`,t(delta.conf)[,2],
#      pch=16,cex=0.5,
#      ylim=c(-1,1),
#      xlab="Sample size (n)",
#      ylab="Comparison of the posteriors for \n complete and full partial pooling parameterizations \n complete-partial (median and 95% CI)")
# segments(x0=sumMLE$`n()`,y0=t(delta.conf)[,1],
#          x1=sumMLE$`n()`,y1=t(delta.conf)[,3])
# abline(h=0,lwd=1)
# dev.off()


## LOGIT PARAMETERIZATION

modelString4 <- "
model { 
  ## Priors
  ##hyperprior for intercept alpha
  
  for(i in 1:n_site){
    sigma.site[i] ~ dunif(0,100)
    tau.site[i] <- 1/(sigma.site[i]*sigma.site[i])
    mu.alpha[i] ~ dnorm(0, 0.0001)
  }
  
  for(i in 1:n_site){
  for(j in 1:n_year){
        #site intercepts
        alpha.i[i,j] ~ dnorm(mu.alpha[site[i]], tau.site[site[i]])
    }
}

  ## Likelihood
  for(i in 1:n){
    theta[i] <- ilogit(alpha.i[site[i],year[i]] )
    fruitingPlantNumber[i] ~ dbin(theta[i], seedlingNumber[i])
  }
  
}
"

library(tidybayes)

dfBayes<-dfAnalysis %>% 
  dplyr::select(-c(transect,plot)) %>%
  dplyr::mutate(year = as.factor(year))

data<-tidybayes::compose_data(dfBayes)


# -------------------------------------------------------------------
# Initial values
# -------------------------------------------------------------------
# set inits for JAGS
inits = list(list(mu.alpha = rep(rnorm(1),data$n_site), 
                  sigma.site = rep(rlnorm(1),data$n_site)),
             list(mu.alpha = rep(rnorm(1),data$n_site), 
                  sigma.site = rep(rlnorm(1),data$n_site)),
             list(mu.alpha = rep(rnorm(1),data$n_site),
                  sigma.site = rep(rlnorm(1),data$n_site))
) 


# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# -------------------------------------------------------------------


# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 500
n.update = 5000
n.iterations = 1000
n.thin = 1
parsToMonitor = c("mu.alpha","sigma.site","theta","alpha.i")
library(rjags)

set.seed(2)
# tuning (n.adapt)
jm4 = jags.model(textConnection(modelString4), data = data, inits=inits,
                 n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm4, n.iterations = n.update)

# chain (n.iter)
samples.rjags4 = coda.samples(jm4, variable.names = c(parsToMonitor), 
                              n.iter = n.iterations, thin = n.thin)

MCMCvis::MCMCsummary(samples.rjags4)

library(tidybayes)

samples.rjags4 %<>% tidybayes::recover_types(dfBayes)

summaryBayes4<-samples.rjags4 %>% 
  tidybayes::recover_types(dfBayes) %>% 
  tidybayes::spread_draws(alpha.i[site,year]) %>%
  dplyr::mutate(theta=boot::inv.logit(alpha.i)) %>%
  dplyr::summarise(bayes4=median(theta)) 

summaryBayes4$year <- as.character(summaryBayes4$year)
summaryBayes4$site <- as.factor(summaryBayes4$site)


ns<-  dfBayes %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n())

ns$site <- as.factor(ns$site)
ns$year <- as.numeric(as.character(ns$year))
summaryBayes4$year <- as.numeric(as.character(summaryBayes4$year))

comp<-  summaryMLE %>%
  dplyr::left_join(summaryBayes4,by=c("site","year")) %>%
  dplyr::left_join(ns,by=c("site","year"))

# plot comparing MLE fit and bayes fit 

setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
pdf(file="appendix-x-mle_bayeslogit.pdf", width=8, height=4)
par(mfrow=c(1,2))
par(mar=c(5, 6, 4, 2) + 0.1)
plot(comp$p,comp$bayes4,
     xlab="Maximum likelihood estimate",
     ylab="Median of posterior for \n beta-binomial parameterization",
     xlim=c(0,1),ylim=c(0,1),
     pch=16,cex=.5)  
abline(a=0,b=1,lwd=0.5)

d<-apply(MCMCvis::MCMCchains(samples.rjags4,params="alpha.i"),2,boot::inv.logit)
sumMLE<-comp %>%
  dplyr::arrange(year)
tmp<-matrix(NA,ncol=200,nrow=30000)
for(i in 1:200){
  tmp[,i]<-d[,i]-sumMLE$p[i]
}
delta<-apply(tmp,2,median)
delta.conf <- apply(tmp,2,quantile,probs=c(.025,.5,.975))
par(mar=c(5, 6, 4, 2) + 0.1)
plot(sumMLE$`n()`,t(delta.conf)[,2],
     pch=16,cex=0.5,
     ylim=c(-1,1),
     xlab="Sample size (n)",
     ylab="Comparison of the maximum likelihood estimate\n and the beta-binomial parameterization\n (median and 95% CI)")
segments(x0=sumMLE$`n()`,y0=t(delta.conf)[,1],
         x1=sumMLE$`n()`,y1=t(delta.conf)[,3])
abline(h=0,lwd=1)

dev.off()


summaryBayes<-samples.rjags1 %>% 
  tidybayes::recover_types(dfBayes) %>% 
  tidybayes::spread_draws(p[site,year]) %>%
  dplyr::summarise(bayes1=median(p)) 

summaryBayes$year <- as.character(summaryBayes$year)
summaryBayes$site <- as.factor(summaryBayes$site)

summaryBayes4<-samples.rjags4 %>% 
  tidybayes::recover_types(dfBayes) %>% 
  tidybayes::spread_draws(alpha.i[site,year]) %>%
  dplyr::mutate(theta=boot::inv.logit(alpha.i)) %>%
  dplyr::summarise(bayes4=median(theta)) 

summaryBayes4$year <- as.character(summaryBayes3$year)
summaryBayes4$site <- as.factor(summaryBayes3$site)

sum <- summaryBayes %>%
  dplyr::left_join(summaryBayes4,by=c("site","year"))

plot(sum$bayes1,sum$bayes4)
abline(a=0,b=1)

tmp<-matrix(NA,ncol=200,nrow=3000)
for(i in 1:200){
  
  tmp[,i]<-MCMCvis::MCMCchains(samples.rjags1,params="p")[,i] - boot::inv.logit(MCMCvis::MCMCchains(samples.rjags4,params="alpha.i")[,i])
  
}

setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
pdf(file="appendix-x-bayescomp_bayeslogit.pdf", width=4, height=4)
delta.conf <- apply(tmp,2,quantile,probs=c(.025,.5,.975))
par(mar=c(5, 6, 4, 2) + 0.1)
plot(sumMLE$`n()`,t(delta.conf)[,2],
     pch=16,cex=0.5,
     ylim=c(-1,1),
     xlab="Sample size (n)",
     ylab="Comparison of the posteriors for \n complete and partial pooling parameterizations \n complete-partial (median and 95% CI)")
segments(x0=sumMLE$`n()`,y0=t(delta.conf)[,1],
         x1=sumMLE$`n()`,y1=t(delta.conf)[,3])
abline(h=0,lwd=1)
dev.off()


v<-(1:200)[sumMLE$`n()`==1][!is.na((1:200)[sumMLE$`n()`==1])]

omegaSite<-c(11,11,13,6,7,11,12,20)


par(mfrow=c(2,4))
for(i in 1:length(v)){
  plot(density(  MCMCvis::MCMCchains(samples.rjags1,params="p")[,v[i]] ),
       xlim=c(0,1),ylim=c(0,4),
       main="") 
  abline(v=median( MCMCvis::MCMCchains(samples.rjags1,params="p")[,v[i]] ))
  tmp <- boot::inv.logit(MCMCvis::MCMCchains(samples.rjags4,params="alpha.i"))[,v[i]]
  lines( density(tmp),lty='dotted')
  abline(v=median( tmp),lty='dotted')
 
} 
