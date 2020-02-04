# -------------------------------------------------------------------
# Fruits per plant (Fec)
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

# custom paste (for NAs)
# https://stackoverflow.com/questions/13673894/suppress-nas-in-paste
paste3 <- function(...,sep=", ") {
  L <- list(...)
  L <- lapply(L,function(x) {x[is.na(x)] <- ""; x})
  ret <-gsub(paste0("(^",sep,"|",sep,"$)"),"",
             gsub(paste0(sep,sep),sep,
                  do.call(paste,c(L,list(sep=sep)))))
  is.na(ret) <- ret==""
  ret
}

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize data
# note; all these data come from permanet plots
# -------------------------------------------------------------------
# -------------------------------------------------------------------
seeds <- read.csv(file = "seeds.csv", header=TRUE)
# -------------------------------------------------------------------
# Calculate seeds per fruit in 2013-2015
# weights to convert damaged to undamaged fruits
# -------------------------------------------------------------------
# mean seeds/fruit in undamaged fruits, by site
summary <- seeds %>%
  dplyr::filter(year > 2012) %>%
  dplyr::group_by(year,site) %>%
  dplyr::summarize(count=n(),
                   meanSeedPerUndFruit = mean(noSeedsPerUndFruit, na.rm = T),
                   meanSeedPerDamFruit = mean(noSeedsPerDamFruit, na.rm = T)) %>%
  dplyr::mutate(weight = meanSeedPerDamFruit/meanSeedPerUndFruit)

# set seed so that rnorm always generates same distribution
# this could be changed for simulations?
set.seed(117)

# randomly sample other data to estimate mean seed per damaged fruit 
# when data isn't available
summary$weight[is.na(summary$weight)] <- rnorm(1,
                                               mean=mean(summary$weight,na.rm=TRUE),
                                               sd=sd(summary$weight,na.rm=TRUE))

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
fruits <- read.csv(file = "fruitsPerPlantFromPlots.csv", header=TRUE)
fruits$position<-as.character(fruits$position)

fruits <- fruits %>%
  dplyr::left_join(summary, by = c("year","site"))
fruits$totalFruitEquivalents <- fruits$noFruitsPerPlant

fruitsConvert <- fruits %>%
  dplyr::filter(year > 2012)
fruitsConvert<-fruitsConvert %>%
  dplyr::mutate(TFE=noUndFruitsPerPlant+noDamFruitsPerPlant*weight)
fruitsConvert$totalFruitEquivalents<-floor(fruitsConvert$TFE)
fruitsConvert <- fruitsConvert %>%
  dplyr::select(-TFE)

fruits <- rbind(filter(fruits,year<2013),fruitsConvert)

# -------------------------------------------------------------------
# Create regression data frame 
# -------------------------------------------------------------------
fecDF <- dplyr::select(fruits, c(site, year, transect, position,
                                 noFruitingPlants,
                                 totalFruitEquivalents)) 
# -------------------------------------------------------------------
# Clean up data for regression
# -------------------------------------------------------------------
# what data to exclude?
# recode noSeedlings==0 & ??
# throwing out 320/9000+ data points
#regressionDf$totalFruitEquivalents[regressionDf$noSeedlings==0] = NA

# unique transect/site variables
fecDF$uniqueTransect<-paste3(fecDF$site,
                            fecDF$transect,sep="")
fecDF$uniquePosition<-paste3(fecDF$site,
                            fecDF$transect,
                            fecDF$position,sep="")

# site, transect and position as factors
fecDF$uniqueTransect <- as.factor(fecDF$uniqueTransect)
fecDF$uniquePosition <- as.factor(fecDF$uniquePosition)
fecDF$otherClarkia <- as.factor(fecDF$otherClarkia)

# -------------------------------------------------------------------
# Save data object
# -------------------------------------------------------------------
setwd("~/Dropbox/clarkiaSeedBanks/library/dataForAnalysis")
save(fecDF, file = "fecDF.RData")