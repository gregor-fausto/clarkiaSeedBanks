# -------------------------------------------------------------------
# Models for estimates of viability
# -------------------------------------------------------------------
# rm(list=ls(all=TRUE)) # clear R environment
rm(list=setdiff(ls(all=TRUE),c("dataDirectory","modelDirectory","fileDirectory","n.adapt","n.update","n.iterations","n.thin"))) # if using in source(script)
options(stringsAsFactors = FALSE)
# -------------------------------------------------------------------
# Load packages required for data cleaning and model fitting
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
# -------------------------------------------------------------------
# Import and organize viability trial data
# -------------------------------------------------------------------
viabilityRawData <- readRDS(paste0(dataDirectory,"viabilityRawData.rds"))

viabilityRawData <- viabilityRawData %>% 
  dplyr::select(-c(germPerc,germNot,viabPerc,viabPerc2,condTest))

viabilityRawData$bag <-as.integer(as.numeric(viabilityRawData$bagNo))

# -------------------------------------------------------------------
## Issue with viabStart coded as NA instead of 0
write.table(viabilityRawData[is.na(viabilityRawData$viabStart),],
            paste0(dataDirectory,"viability.txt"),sep="\t",row.names=FALSE)
write("Issue with 1 row of data; number of seeds starting viability trials is coded viabStart=NA and viabStain=NA;  all others have viabStart=0 and viabStain=NA", 
      file= paste0(dataDirectory,"viability-metadata.txt"));

# recode
viabilityRawData[is.na(viabilityRawData$viabStart),]$viabStart = 0
# -------------------------------------------------------------------

# -------------------------------------------------------------------
## Issue with 23 rows of data where more seeds started the viability trials
## than were left after the end of the germination trials; remove these data
write.table(viabilityRawData %>% dplyr::filter(germStart - germCount - viabStart<0) ,
            paste0(dataDirectory,"viability.txt"),sep="\t",
            row.names=FALSE,append=TRUE,col.names=FALSE)

#fileConn<-file();
write(c("Issue with 23 rows of data; More seeds started the viability trials than were left after the end of the germination trials; remove these data"), 
      file=paste0(dataDirectory,"viability-metadata.txt"),append=TRUE);

# filter to remove these data
viabilityRawData<-viabilityRawData %>% 
  dplyr::filter(germStart - germCount - viabStart >= 0)
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Prepare data for analysis with JAGS
# -------------------------------------------------------------------

viabilityRawData = viabilityRawData %>%
  dplyr::mutate(year = as.factor(round)) %>%
  dplyr::select(site, year, age, germStart, germCount, viabStart, viabStain, bagNo) %>%
  dplyr::rename(siteViab = site,
                yearViab = year,
                ageViab = age,
                bag = bagNo) %>%
  dplyr::mutate(bag = as.factor(bag)) %>%
  dplyr::group_by(siteViab,yearViab, ageViab, bag) %>%
  # sum observations in each bag; this is ignoring some variation
  dplyr::summarise(germStart = sum(germStart,na.rm=TRUE),
                   germCount = sum(germCount,na.rm=TRUE),
                   viabStart = sum(viabStart,na.rm=TRUE),
                   viabStain = sum(viabStain,na.rm=TRUE)) 

# indexing to avoid issue of dealing with ragged arrays in JAGS
experimentalAge = c(1,1,1,2,2,3)
experimentalRound = c(1,2,3,1,2,1)
experimentalIndex = 1:6

df.index=data.frame(ageViab=experimentalAge,yearViab=experimentalRound,indexViab=experimentalIndex)
viabilityRawData <- viabilityRawData %>%
  dplyr::mutate(yearViab=as.numeric(yearViab)) %>%
  dplyr::left_join(df.index,by=c("ageViab","yearViab")) %>%
  dplyr::mutate(yearViab=as.factor(yearViab))

data <- tidybayes::compose_data(viabilityRawData)
data$experimentalIndex = experimentalIndex

detach("package:tidyverse", unload=TRUE)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# scalars that specify the
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# n.adapt = 3000
# n.update = 5000
# n.iterations = 10000
# n.thin = 1
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Set initials
# -------------------------------------------------------------------
# -------------------------------------------------------------------

initsMu0 <- function(rows = data$n_siteViab, cols = data$n_yearViab){
  matrix( rnorm(n = rows*cols, mean = 0, sd = 1), rows, cols)
}

initsSigma0 <- function(rows = data$n_siteViab, cols = data$n_yearViab){
 # extraDistr::rhnorm(n = samps, sigma = 1)
  matrix( extraDistr::rhnorm(n = rows*cols, sigma = 1), rows, cols)
}

# initsMu <- function(samps = data$n_siteViab){
#   rnorm(n = samps, mean = 0, sd = 1)
# }

initsSigma <- function(rows = data$n_siteViab, cols = 6){
  matrix( extraDistr::rhnorm(n = rows*cols, sigma = 1), rows, cols)
}

# set inits for JAGS
# inits <- list()
# for(i in 1:3){
#   inits[[i]] <- list(initsMu0(), initsMu0(), initsSigma0(), initsSigma0(), initsSigma(), initsSigma())
# 
#   names(inits[[i]]) = c("mu0_g","mu0_v",
#                         "sigma0_g","sigma0_v",
#                         "sigma_g","sigma_v")
# 
# }

# jm = jags.model("~/Dropbox/clarkiaSeedBanks/scriptsModelFitting/jagsScripts/jags-viability.R", data = data,
#                 n.chains = 3, n.adapt = n.adapt, inits = inits)

# -------------------------------------------------------------------
# Set up JAGS for parallel run
# -------------------------------------------------------------------
library(stringr)

cl <- makeCluster(3)

myWorkers <- NA
for(i in 1:3) myWorkers[i] <- stringr::word(capture.output(cl[[i]]), -1)

initsFun = function(){temp = list(initsMu0(),initsMu0(), initsSigma0(), initsSigma0(), initsSigma(), initsSigma(),
                                  "base::Mersenne-Twister", runif(1, 1, 2000))

names(temp) = c("mu0_g","mu0_v",
                "sigma0_g","sigma0_v",
                "sigma_g","sigma_v",
                ".RNG.name",".RNG.seed")
return(temp)
}

inits = list(initsFun(),initsFun(),initsFun())

 parsToMonitor_g = c("mu0_g","sigma0_g","mu_g","sigma_g","germCount_sim","chi2.germCount.obs","chi2.germCount.sim")
 parsToMonitor_v = c("mu0_v","sigma0_v","mu_v","sigma_v","viabStain_sim","chi2.viabStain.obs","chi2.viabStain.sim")

parallel::clusterExport(cl, c("myWorkers","data", "inits", "n.adapt", "n.update",
                              "n.iterations", "parsToMonitor_g", "parsToMonitor_v","modelDirectory"))


out <- clusterEvalQ(cl, {
  library(rjags)
  jm = jags.model(file=paste0(modelDirectory,"jags-viability.R"),
                  data = data, n.chains = 1,
                  n.adapt = n.adapt, 
                  inits = inits[[which(myWorkers==Sys.getpid())]])
  update(jm, n.iter = n.update)
  zm = coda.samples(jm, variable.names = c(parsToMonitor_g,parsToMonitor_v),
                    n.iter = n.iterations, thin = 1)
  return(as.mcmc(zm))
})
stopCluster(cl)
samples.rjags = mcmc.list(out)

dir.create(file.path(fileDirectory), showWarnings = FALSE)
saveRDS(samples.rjags,file=paste0(fileDirectory,"viabilityTrialSamples.rds"))
saveRDS(data,file=paste0(fileDirectory,"viabilityData.rds"))
