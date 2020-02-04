# -------------------------------------------------------------------
# SEEDS PER FRUIT (PHI)
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)
# -------------------------------------------------------------------
# Set appropriate directories
# -------------------------------------------------------------------
setwd("~/Dropbox/projects/clarkiaScripts/data/reshapeData/")
# -------------------------------------------------------------------
# Load packages
# -------------------------------------------------------------------
library(dplyr) # manipulate data frames
library(tidyr) # manipulate data frames
library(reshape2) # manipulate data frames
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize data
# note; all these data come from permanet plots
# -------------------------------------------------------------------
# -------------------------------------------------------------------
seeds <- read.csv(file = "seeds.csv", header=TRUE)
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# visualize
# -------------------------------------------------------------------
# -------------------------------------------------------------------
par(mfrow=c(2,1))
hist(seeds$noSeedsPerFruit)
hist(seeds$noSeedsPerUndFruit,add=TRUE,col="red")
hist(seeds$noSeedsPerDamFruit,add=TRUE,col="blue")
hist(summarySeeds$meanSeedsPerFruitEq)
hist(summarySeeds$meanSeedsPerUndFruit,add=TRUE,col="red")
hist(summarySeeds$meanSeedsPerDamFruit,add=TRUE,col="blue")

plot(density(seeds$noSeedsPerFruit,na.rm=TRUE),ylim=c(0,0.08))
lines(density(seeds$noSeedsPerUndFruit,na.rm=TRUE),col="red")
lines(density(seeds$noSeedsPerDamFruit,na.rm=TRUE),col="blue")
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Gather
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# question whether or not to just include from permanent plot
df<-seeds %>% dplyr::select(-permanentPlot) %>% tidyr::gather(variable, response, -c(year,site))
df<-df[complete.cases(df), ]

# -------------------------------------------------------------------
# Create regression data frame 
# -------------------------------------------------------------------
regressionDf <- dplyr::select(df, c(site, year,  variable, response)) 

# -------------------------------------------------------------------
# Save data object
# -------------------------------------------------------------------
setwd("~/Dropbox/clarkiaSeedBanks/library/dataForAnalysis")
phiIndDF <- regressionDf
save(phiIndDF, file = "phiIndDF.RData")
