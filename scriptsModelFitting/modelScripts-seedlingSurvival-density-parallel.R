# -------------------------------------------------------------------
# Models for seedling survival to fruiting
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
# rm(list=setdiff(ls(all=TRUE),c(dataDirectory,modelDirectory,fileDirectory,n.adapt,n.update,n.iterations,n.thin))) # if using in source(script)
options(stringsAsFactors = FALSE)
# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(rjags) 
library(tidybayes)
library(tidyverse)
library(parallel)

# -------------------------------------------------------------------
# Set directories
# -------------------------------------------------------------------
dataDirectory = "/Users/Gregor/Dropbox/dataLibrary/postProcessingData/"
modelDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsModelFitting/jagsScripts/"
fileDirectory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/parallel/"

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# setwd and read data files
censusSeedlingsFruitingPlants <- readRDS(paste0(dataDirectory,"censusSeedlingsFruitingPlants.RDS"))

# -------------------------------------------------------------------
# ## Issue with number of fruiting plants greater than number of seedlings
# write.table(censusSeedlingsFruitingPlants %>% 
#               dplyr::filter((fruitplNumber>seedlingNumber)) ,
#             "~/Dropbox/clarkiaSeedBanks/dataToCheck/seedlingSurvival.txt",sep="\t",row.names=FALSE)
# write("Issue with rows of data where number of fruiting plants is greater than number of seedlings", 
#       file="~/Dropbox/clarkiaSeedBanks/dataToCheck/seedlingSurvival-metadata.txt");

# recode all data
# filter out all data with issues
censusSeedlingsFruitingPlants <- censusSeedlingsFruitingPlants %>% 
  # more fruiting plants than seedlings
  # recode s.t. number of seedlings equals number of fruiting plants
  dplyr::mutate(seedlingNumber=ifelse(fruitplNumber>seedlingNumber,fruitplNumber,seedlingNumber)) %>% 
  # NA for seedlings
  # filter these out; these are true missing data
  # recode s.t. number of seedlings equals number of fruiting plants
  #dplyr::mutate(seedlingNumber=ifelse(is.na(seedlingNumber),fruitplNumber,seedlingNumber)) %>%
  dplyr::filter(!is.na(seedlingNumber)) %>% 
  # NA for fruiting plants
  # still filter these out; missing response
  dplyr::filter(!is.na(fruitplNumber)) 
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------
# -------------------------------------------------------------------
censusSeedlingsFruitingPlants$year <- as.character(censusSeedlingsFruitingPlants$year)
censusSeedlingsFruitingPlants$density=censusSeedlingsFruitingPlants$seedlingNumber
tmp <- censusSeedlingsFruitingPlants %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(mu=mean(density),sd=sd(density))

censusSeedlingsFruitingPlants <- censusSeedlingsFruitingPlants %>% 
  dplyr::left_join(tmp,by=c('site')) %>%
  dplyr::mutate(density.std=(density-mu)/(2*sd)) %>%
  dplyr::select(-mu,sd)

cenusSeedlingsFruitingPlants <- censusSeedlingsFruitingPlants %>%
  dplyr::filter(site=="BR")

data <- tidybayes::compose_data(censusSeedlingsFruitingPlants)

detach("package:tidyverse", unload=TRUE)

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
n.iterations = 5000
n.thin = 1

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Describe functions
# -------------------------------------------------------------------
# -------------------------------------------------------------------

initsMu0 <- function(samps = data$n_site){
  rnorm(n = samps, mean = 0, sd = 1)
}

initsSigma0 <- function(samps = data$n_site){
  extraDistr::rhnorm(n = samps, sigma = 1)
}

initsSigma <- function(rows = data$n_site, cols = data$n_year){
  matrix(extraDistr::rhnorm(n = rows*cols, sigma = 1), rows, cols)
}

# # set inits for JAGS
inits <- list()
for(i in 1:3){
inits[[i]] <- list(initsMu0(), initsSigma0(), initsSigma() ,
                   initsMu0(), initsSigma0(), initsSigma())

names(inits[[i]]) = c("mu0","sigma0","sigma",
                      "mu0_b","sigma0_b","sigma_b")

}

# # Call to JAGS
#
# # tuning (n.adapt)
jm = jags.model("~/Dropbox/clarkiaSeedBanks/scriptsModelFitting/jagsScripts/jags-seedlingSurvival-density.R",
data = data, inits = inits,
n.chains = length(inits), n.adapt = n.adapt)


# burn-in (n.update)
update(jm, n.iter = n.update)

parsToMonitor = c("mu0","sigma0","mu","sigma",
                  "mu0_b","sigma0_b","mu_b","sigma_b")
parsToCheck = c(# LOG LIKELIHOODS
  # "logLik",
  # POSTERIOR PREDICTIVE
  "fruitplNumber_sim",
  "e.mu","mu",
  # MODEL CHECKING 
  "chi2.obs","chi2.sim")


# chain (n.iter)
samples.rjags = coda.samples(jm, 
                             variable.names = c(parsToMonitor), 
                             n.iter = n.iterations, thin = n.thin)

mu=MCMCchains(samples.rjags,"mu")
beta=MCMCchains(samples.rjags,"mu_b")
par(mfrow=c(4,5),mar=c(0,0,0,0))
siteNames = unique(censusSeedlingsFruitingPlants$site)
for(k in 1:20){
  data.dens.all=data$density.std[data$site==k]
  n.x=seq(min(data.dens.all),max(data.dens.all),.1)
  mat=matrix(NA,nrow=length(n.x),ncol=13)
  plot(NA,NA,type='n',xlim=c(min(n.x),max(n.x)),ylim=c(0,1),xaxt='n',yaxt='n')
  text(.9*max(n.x),.9,siteNames[k])
  for(j in 1:13){
    data.dens=data$density.std[data$year==j&data$site==k]
    index = grep(paste0("\\[",k,",",j,"\\]"),colnames(mu))
    for(i in 1:length(n.x)){
      vec=boot::inv.logit(mu[,index]+beta[,index]*n.x[i])
      mat[i,j]=t(quantile(vec,c(.5)))
    }
    lines(n.x,mat[,j],type='l')
    }
}


mu0=MCMCchains(samples.rjags,"mu0")
beta0=MCMCchains(samples.rjags,"mu0_b")

par(mfrow=c(4,5),mar=c(0,0,0,0))
for(k in 1:20){
  data.dens.all=data$density.std[data$site==k]
  n.x=seq(min(data.dens.all),max(data.dens.all),.1)
  mat=matrix(NA,nrow=length(n.x),ncol=3)
  plot(NA,NA,type='n',xlim=c(min(n.x),max(n.x)),ylim=c(0,1),xaxt='n',yaxt='n')
      data.dens=data$density.std[data$site==k]
    index = grep(paste0("\\[",k,"\\]"),colnames(mu0))
    for(i in 1:length(n.x)){
      vec=boot::inv.logit(mu0[,index]+beta0[,index]*n.x[i])
      mat[i,1:3]=t(quantile(vec,c(.025,.5,.975)))
    }
    polygon(c(n.x,rev(n.x)),c(mat[,1],rev(mat[,3])),col='gray90',border="gray90")
    lines(n.x,mat[,2],type='l',lwd=2)
    text(.9*max(n.x),.9,siteNames[k])  
    
  
}


mu0=MCMCchains(samples.rjags,"mu0")
beta0=MCMCchains(samples.rjags,"mu0_b")

dev.off()
par(mfrow=c(1,1))
plot(NA,NA,type='n',xlim=c(-.25,0),ylim=c(0,1),xaxt='n',yaxt='n')
mat=matrix(NA,nrow=3,ncol=20)

for(k in 1:20){
  data.dens.all=data$density.std[data$site==k]
  n.x=seq(min(data.dens.all),max(data.dens.all),.1)
  data.dens=data$density.std[data$site==k]
  index = grep(paste0("\\[",k,"\\]"),colnames(mu0))
  for(i in 1){
    vec=boot::inv.logit(mu0[,index]+beta0[,index]*n.x[i])
    mat[1:3,k]=t(quantile(vec,c(.025,.5,.975)))
  }
  }
mat  

df=data.frame(position,lo=(mat[1,]),mid=(mat[2,]),high=(mat[3,]))
par(mfrow=c(1,2))
plot(df$easting,df$mid,pch=19,ylim=c(0,1))
segments(x0=df$easting,y0=df$lo,y1=df$high)


stats=apply(boot::inv.logit(beta0),2,quantile,(c(.025,.5,.975)))
climate=readRDS("~/Dropbox/clarkiaSeedBanks/scriptsAnalysis/climateData.RDS")
head(climate)
position = climate %>% ungroup %>%
  dplyr::filter(intenseDemography==1) %>%
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000) %>%
  unique

df=data.frame(position,lo=(stats[1,]),mid=(stats[2,]),high=(stats[3,]))
par(mfrow=c(1,2))
plot(df$easting,df$mid,pch=19,ylim=c(0,1))
segments(x0=df$easting,y0=df$lo,y1=df$high)


stats=apply(boot::inv.logit(mu0),2,quantile,(c(.025,.5,.975)))

df=data.frame(position,lo=(stats[1,]),mid=(stats[2,]),high=(stats[3,]))
plot(df$easting,df$mid,pch=19,ylim=c(0,1))
segments(x0=df$easting,y0=df$lo,y1=df$high)
# -------------------------------------------------------------------
# Set up JAGS for parallel run
# -------------------------------------------------------------------
library(stringr)

cl <- makeCluster(3)

myWorkers <- NA
for(i in 1:3) myWorkers[i] <- stringr::word(capture.output(cl[[i]]), -1)

initsFun = function(){temp = list(initsMu0(), initsSigma0(), initsSigma(), 
                                  "base::Mersenne-Twister", runif(1, 1, 2000))

names(temp) = c("mu0",
                "sigma0",
                "sigma",
                ".RNG.name",".RNG.seed")
return(temp)
}

inits = list(initsFun(),initsFun(),initsFun())

parsToMonitor = c("mu0","sigma0","mu","sigma")
parsToCheck = c("fruitplNumber_sim","chi2.obs","chi2.sim")

parallel::clusterExport(cl, c("myWorkers","data", "inits", "n.adapt", "n.update",
                              "n.iterations", "parsToMonitor", "parsToCheck"))

out <- clusterEvalQ(cl, {
  library(rjags)
  jm = jags.model(file=paste0(modelDirectory,"jags-seedlingSurvival.R"),
                  data = data, n.chains = 1,
                  n.adapt = n.adapt, 
                  inits = inits[[which(myWorkers==Sys.getpid())]])
  update(jm, n.iter = n.update)
  zm = coda.samples(jm, variable.names = c(parsToMonitor,parsToCheck),
                    n.iter = n.iterations, thin = 1)
  return(as.mcmc(zm))
})
stopCluster(cl)
samples.rjags = mcmc.list(out)

dir.create(file.path(fileDirectory), showWarnings = FALSE)
saveRDS(samples.rjags,file=paste0(fileDirectory,"seedlingSurvivalSamples-density.rds"))
saveRDS(data,file=paste0(fileDirectory,"seedlingData-density.rds"))
