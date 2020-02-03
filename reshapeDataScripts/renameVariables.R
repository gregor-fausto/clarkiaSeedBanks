# -------------------------------------------------------------------
# Rename Variable Names
# -------------------------------------------------------------------
# rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)
# -------------------------------------------------------------------
# Set appropriate directories
# -------------------------------------------------------------------
setwd("~/Dropbox/projects/clarkiaScripts/stanAnalysis/cleanData")
# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(rstan) # stan interface
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
load(file="sigmaDF.RData")
load(file="fecDF.RData")
load(file="phiIndDF.RData")
load(file="s1g1DF.RData")

reassign<-function(df){
  df$site<-ifelse(df$site=="BG","BGU",as.character(df$site))
  df$site<-ifelse(df$site=="BR","BRD",as.character(df$site))
  df$site<-ifelse(df$site=="CF","CFL",as.character(df$site))
  df$site<-ifelse(df$site=="EC","ECR",as.character(df$site))
  df$site<-ifelse(df$site=="FR","FWR",as.character(df$site))
  df$site<-ifelse(df$site=="LO","LOA",as.character(df$site))
  df$site<-ifelse(df$site=="MC","MCR",as.character(df$site))
  df$site<-ifelse(df$site=="SM","SMT",as.character(df$site))
  df$site<-ifelse(df$site=="OKRE","OKE",as.character(df$site))
  df$site<-ifelse(df$site=="OKRW","OKW",as.character(df$site))
  return(df)
}

sigmaDF$site<-as.factor(reassign(sigmaDF)$site)
fecDF$site<-as.factor(reassign(fecDF)$site)
phiIndDF$site<-as.factor(reassign(phiIndDF)$site)
s1g1DF$site<-as.factor(reassign(s1g1DF)$site)

# -------------------------------------------------------------------
# Save data object
# -------------------------------------------------------------------
setwd("~/Dropbox/projects/clarkiaScripts/stanAnalysis/cleanData")
save(sigmaDF, file = "sigmaDF.RData")
save(fecDF, file = "fecDF.RData")
save(phiIndDF, file = "phiIndDF.RData")
save(s1g1DF, file = "s1g1DF.RData")
