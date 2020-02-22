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

yd<-split(dfAnalysis,list(dfAnalysis$site,dfAnalysis$year))
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


dfSite <- read.csv(file="~/Downloads/Datafile_Site_Environment_corr.csv",header=TRUE)

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

samples.rjags1 %<>% tidybayes::recover_types(dfBayes)

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
  
comp<-  summaryMLE %>%
    dplyr::left_join(summaryBayes,by=c("site","year")) %>%
  dplyr::left_join(ns,by=c("site","year"))

comp %>% dplyr::mutate(abs(p-bayes1)) %>% View

plot(comp$p,comp$bayes1,xlim=c(0,1),ylim=c(0,1))  
abline(a=0,b=1)
