# -------------------------------------------------------------------
# Models for fitting data on fruits per plant and seeds per fruit
# -------------------------------------------------------------------
# rm(list=ls(all=TRUE)) # clear R environment
rm(list=setdiff(ls(all=TRUE),c("dataDirectory","modelDirectory","fileDirectory","n.adapt","n.update","n.iterations","n.thin"))) # if using in source(script)
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
# dataDirectory = "/Users/Gregor/Dropbox/dataLibrary/postProcessingData/"
# modelDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsModelFitting/jagsScripts/"
# fileDirectory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/parallel/"

# -------------------------------------------------------------------
# Import and organize data
# -------------------------------------------------------------------

countSeedPerFruit <-readRDS(paste0(dataDirectory,"countSeedPerFruit.RDS"))

countSeedPerUndamagedFruit <- countSeedPerFruit %>%
  dplyr::filter(demography==1) %>%
  dplyr::filter(damaged==0) %>%
  dplyr::rename(site3 = site) %>%
  dplyr::rename(year3 = year) %>%
  dplyr::select(site3,year3,sdno) %>% 
  dplyr::filter(!is.na(sdno))

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
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------
data <- tidybayes::compose_data(countSeedPerUndamagedFruit,
                                countSeedPerDamagedFruit)

data$n3 = dim(countSeedPerUndamagedFruit)[1]
data$n4 = dim(countSeedPerDamagedFruit)[1]

detach("package:tidyverse", unload=TRUE)

# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
# -------------------------------------------------------------------

# n.adapt = 3000
# n.update = 5000
# n.iterations = 5000
# n.thin = 1

# -------------------------------------------------------------------
# Describe functions
# -------------------------------------------------------------------

initsMu <- function(rows = data$n_site3, cols = data$n_year3){
  matrix(rgamma(n = rows*cols, shape = 1, rate = 1), rows, cols)
}

initsSigma <- function(rows = data$n_site3, cols = data$n_year3){
  matrix(extraDistr::rhnorm(n = rows*cols, sigma = 2), rows, cols)
}

# # set inits for JAGS
#  inits <- list()
# for(i in 1:3){
#   inits[[i]] <- list(initsMu0(), initsMu0(),
#                      initsSigma0(), initsSigma0(),
#                      initsSigma(rows = data$n_site3,cols=data$n_year3),
#                      initsSigma(rows = data$n_site3,cols=data$n_year4) )
# 
#   names(inits[[i]]) = c(paste(rep("nu",2),c("seeds","dam_seeds"),sep="_"),
#                         paste(rep("sigma0",2),c("seeds","dam_seeds"),sep="_"),
#                         paste(rep("sigma",2),c("seeds","dam_seeds"),sep="_"))
# 
# }

# Call to JAGS

# # tuning (n.adapt)
# jm = jags.model("~/Dropbox/clarkiaSeedBanks/scriptsModelFitting/jagsScripts/jags-seedsperfruit.R",
#                 data = data, inits = inits,
#                 n.chains = length(inits), n.adapt = n.adapt)
# 
# 
# # burn-in (n.update)
# update(jm, n.iter = n.update)
# 
# 
# parsToMonitor = c(paste(rep("nu",2),c("seeds","dam_seeds"),sep="_"),
#                   paste(rep("sigma0",2),c("seeds","dam_seeds"),sep="_"),
#                  # paste(rep("mu",2),c("seeds","dam_seeds"),sep="_"),
#                   paste(rep("sigma",2),c("seeds","dam_seeds"),sep="_"),
#                   paste(rep("mu.log",2),c("seeds","dam_seeds"),sep="_"))
# # parsToCheck = c("y_tfe.sim","chi2.tfe.obs","chi2.tfe.sim",
#                        # "y_und.sim","chi2.und.obs","chi2.und.sim",
#                        # "y_dam.sim","chi2.dam.obs","chi2.dam.sim",
#                        # "y_sd.sim","chi2.sd.obs","chi2.sd.sim",
#                        # "y_sd_dam.sim","chi2.sd_dam.obs","chi2.sd_dam.sim")
# 
# 
# # chain (n.iter)
# samples.rjags = coda.samples(jm, 
#                              variable.names = c(parsToMonitor), 
#                              n.iter = n.iterations, thin = n.thin)

# -------------------------------------------------------------------
# Set up JAGS for parallel run
# -------------------------------------------------------------------
library(stringr)

cl <- makeCluster(3)

myWorkers <- NA
for(i in 1:3) myWorkers[i] <- stringr::word(capture.output(cl[[i]]), -1)

initsFun = function(){temp = list(initsMu(rows = data$n_site3,cols=data$n_year3),
                                  initsMu(rows = data$n_site4,cols=data$n_year4),
                                  initsSigma(rows = data$n_site3,cols=data$n_year3),
                                  initsSigma(rows = data$n_site4,cols=data$n_year4),
                                  "base::Mersenne-Twister", runif(1, 1, 2000))

names(temp) = c(paste(rep("mu.log",2),c("seeds","dam_seeds"),sep="_"),
                paste(rep("sigma",2),c("seeds","dam_seeds"),sep="_"),
                ".RNG.name",".RNG.seed")
return(temp)
}


inits = list(initsFun(),initsFun(),initsFun())

parsToMonitor = c(paste(rep("mu.log",2),c("seeds","dam_seeds"),sep="_"),
                  paste(rep("sigma",2),c("seeds","dam_seeds"),sep="_"))
parsToCheck = c("y_sd.sim","chi2.sd.obs","chi2.sd.sim",
                "y_sd_dam.sim","chi2.sd_dam.obs","chi2.sd_dam.sim")

parallel::clusterExport(cl, c("myWorkers","data", "inits", "n.adapt", "n.update",
                              "n.iterations", "parsToMonitor", "parsToCheck", "modelDirectory"))

out <- clusterEvalQ(cl, {
  library(rjags)
  jm = jags.model(file=paste0(modelDirectory,"jags-seedsperfruit-noPool.R"),
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
saveRDS(samples.rjags,file=paste0(fileDirectory,"seedsPerFruitSamples-noPool.rds"))
saveRDS(data,file=paste0(fileDirectory,"seedsPerFruitData-noPool.rds"))



# 
# library(MCMCvis)
# 
# # nu
# MCMCsummary(samples.rjags,params="nu_seeds")
# par(mfrow=c(1,2),mar=c(2,2,2,2))
# for(i in 1:2){hist(exp(MCMCchains(samples.rjags,params="nu_seeds")[,i]),
#                    main='',breaks=50,col='gray90',border='white')}
# par(mfrow=c(1,1),mar=c(2,2,2,2))
# hist(MCMCchains(samples.rjags,params="nu_seeds")[,1],
#      main='',breaks=50,col='gray90',border='white')
# hist(rgamma(10000,1,1),
#      main='',breaks=50,col='gray90',border='white')
# 
# # sigma0
# MCMCsummary(samples.rjags,params="sigma0_seeds")
# par(mfrow=c(1,1),mar=c(2,2,2,2))
# hist(MCMCchains(samples.rjags,params="sigma0_seeds")[,2],
#      main='',breaks=50,col='gray90',border='white')
# hist(extraDistr::rhnorm(n = 10000, sigma = 1),
#      main='',breaks=50,col='gray90',border='white')
# 
# # mu vs. mu.log
# MCMCsummary(samples.rjags,params="mu.log_seeds")
# 
# hist((MCMCchains(samples.rjags,params="mu.log_seeds")[,1]),
#      main='',breaks=50,col='gray80',border='white')
# hist(log(rlnorm(10000,rgamma(10000,1,1),1)),
#      main='',breaks=50,col='gray80',border='white')
# hist(rgamma(1000,1,1))
# par(mfrow=c(4,4),mar=c(2,2,2,2))
# tmp=countSeedPerUndamagedFruit %>%
#   dplyr::group_by(site3,year3) %>%
#   dplyr::summarise(mu =median(sdno,na.rm=TRUE),
#                    n=n())
# bind_rows()
# site1=grep(paste0("\\[",2,","),colnames(MCMCchains(samples.rjags,params="mu.log_seeds")))
# for(i in 14:23){hist(exp(MCMCchains(samples.rjags,params="mu.log_seeds")[,site1[i-13]]),
#                      main='',breaks=50,col='gray80',border='white');
#   abline(v=(tmp[i,]$mu),lty='dotted')
#   text(x=(tmp[i,]$mu)*.98,y=100,tmp$n[i],col='red')
# }
# 
# # mu vs. mu.log
# MCMCsummary(samples.rjags,params="mu.log_dam_seeds")
# 
# par(mfrow=c(3,2),mar=c(2,2,2,2))
# tmp=countSeedPerDamagedFruit %>%
#   dplyr::group_by(year4) %>%
#   dplyr::summarise(mu =median(sdno_dam,na.rm=TRUE),
#                    n=n())
# for(i in 1:12){hist(exp(MCMCchains(samples.rjags,params="mu.log_dam_seeds")[,i]),
#                     main='',breaks=50,col='gray80',border='white');
#   abline(v=(tmp[i,]$mu),lty='dotted');abline(v=5,col='red',lty='dotted')
#   text(x=(tmp[i,]$mu)*.98,y=100,tmp$n[i],col='red')
# }
# 
# dev.off()
# 
# hist(exp(MCMCchains(samples.rjags,params="mu.log_dam_seeds")[,3]),
#      main='',breaks=50,col='gray80',border='white')
# plot(exp(MCMCchains(samples.rjags,params="mu.log_dam_seeds")[,3]))