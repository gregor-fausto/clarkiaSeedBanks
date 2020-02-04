# -------------------------------------------------------------------
# Checking individual models
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
library(bayesplot)
library(stringr)

set.seed(10)
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
setwd("~/Dropbox/projects/clarkiaScripts/data/cleanData")
load(file="sigmaDF.RData") 

setwd("~/Dropbox/projects/clarkiaScripts/data/cleanData")
load(file="fecDF.RData")

setwd("~/Dropbox/projects/clarkiaScripts/data/cleanData")
load(file="phiIndDF.RData")

### MULTIPLE SITE, MULTIPLE YEAR MODEL (POOLED)

siteFilter = as.character(unique(sigmaDF$site))[1:20]
yearFilter = 2007:2012

# need to exlude 2 observations if using all data points without latent states
sigmaDat <- sigmaDF %>% 
  dplyr::filter(year %in% yearFilter & site %in% siteFilter) %>%
  dplyr::select(c(year,site,transect,position,noSeedlings,noFruitingPlants)) %>%
  dplyr::mutate( p = noFruitingPlants/noSeedlings ) %>%
  dplyr::mutate( test = ifelse(p>1, 1, 0)) %>%
  dplyr::filter( test < 1 | is.na(test)) %>% 
  dplyr::filter( !is.na(noSeedlings)) %>%
  dplyr::select(c(site, year, transect, position, noSeedlings, noFruitingPlants)) %>%
  unique

# remove plots without plants (NA)
fecDat <- fecDF %>% 
  dplyr::filter(year %in% yearFilter & site %in% siteFilter) %>%
  dplyr::select(c(year,site,transect,position,totalFruitEquivalents)) %>%
  dplyr::filter( !is.na(totalFruitEquivalents) )

# remove plots without plants (NA)
phiDat <- phiIndDF %>% 
  dplyr::filter(year %in% yearFilter & site %in% siteFilter) %>%
  dplyr::select(c(year,site,response))

# index 
sigmaDat$site.index <- as.integer(as.factor(sigmaDat$site))
sigmaDat$year.index <- as.integer(as.factor(sigmaDat$year))

fecDat$site.index <- as.integer(as.factor(fecDat$site))
fecDat$year.index <- as.integer(as.factor(fecDat$year))

phiDat$site.index <- as.integer(as.factor(phiDat$site))
phiDat$year.index <- as.integer(as.factor(phiDat$year))

if(length(unique(sigmaDat$year.index)) != length(unique(fecDat$year.index))) stop("unequal years")
if(length(unique(sigmaDat$year.index)) != length(unique(phiDat$year.index))) stop("unequal years")

if(length(unique(sigmaDat$site.index)) != length(unique(fecDat$site.index))) stop("unequal sites")
if(length(unique(sigmaDat$site.index)) != length(unique(phiDat$site.index))) stop("unequal sites")


data = list(
  n = as.double(sigmaDat$noSeedlings),
  y = as.double(sigmaDat$noFruitingPlants),
  N = nrow(sigmaDat),
  siteSigma = as.double(sigmaDat$site.index),
  yearSigma = as.double(sigmaDat$year.index),
  y1 = as.double(fecDat$totalFruitEquivalents),
  N1 = nrow(fecDat),
  siteFec = as.double(fecDat$site.index),
  yearFec = as.double(fecDat$year.index),
  y2 = as.double(phiDat$response),
  N2 = nrow(phiDat),
  sitePhi = as.double(phiDat$site.index),
  yearPhi = as.double(phiDat$year.index),
  nyears = length(unique(sigmaDat$year.index)),
  nsites = length(unique(sigmaDat$site.index))
)


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Evaluate convergence
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# intercepts = c("alphaS1","alphaG1","alphaS2")
# slopes = c("betaS1","betaG1","betaS2")
# variances = c("sigmaS1","sigmaG1","sigmaS2")
# viab = c("mu.b","sigma.b","alpha")
# sims = c("yv.sim","yt.sim","yg.sim","yo.sim")

# summary table
# MCMCsummary(zc, n.eff = TRUE, params = c(intercepts,slopes,variances))
# not all Rhat statistics are under 1.05; need more iterations
# what is a reasonable sample size? right now in the 100s

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Posterior Predictive Checks
# -------------------------------------------------------------------
# -------------------------------------------------------------------
color_scheme_set("brightblue")

library(gridExtra)


pdf(
  "~/Dropbox/modelsF2019/writing/fitnessPosteriorPredictiveChecks.pdf",
  onefile=TRUE,
  paper="USr",
  height = 5.5, width = 10)

load("~/Dropbox/modelsF2019/output/sigmaFit")

iter<-dim(MCMCchains(zmP,params=
                       "y.sim"))[1]

# Seedling survival to fruiting
ppc_dens_overlay(data$y, MCMCchains(zmP,params="y.sim")[sample(iter,1000), ]) +
  theme_bw() + labs(title="Posterior predictive checks for seedling survival", 
                                     caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$y, MCMCchains(zmP,params="y.sim")[sample(iter,1000), ],group=sigmaDat$site) +
  theme_bw() + labs(title="Posterior predictive checks for mean of seedling survival", 
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")  

ppc_stat_grouped(data$y, MCMCchains(zmP,params="y.sim")[sample(iter,1000), ],group=sigmaDat$site, stat="var") +
  theme_bw() + labs(title="Posterior predictive checks for variance of seedling survival", 
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

# Fruits per plant

load("~/Dropbox/modelsF2019/output/fecFit")

iter<-dim(MCMCchains(zmP,params=
                       "y1.sim"))[1]

ppc_dens_overlay(data$y1, MCMCchains(zmP,params="y1.sim")[sample(iter,1000), ]) +
  theme_bw() + labs(title="Posterior predictive checks for fruits per plant", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$y1, MCMCchains(zmP,params="y1.sim")[sample(iter,1000), ],group=fecDat$site) +
  theme_bw() + labs(title="Posterior predictive checks for the mean of fruits per plant", 
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior") 

ppc_stat_grouped(data$y1, MCMCchains(zmP,params="y1.sim")[sample(iter,1000), ],group=fecDat$site, stat="var") +
  theme_bw() + labs(title="Posterior predictive checks for the variance of fruits per plant", 
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior") 

#Seeds per fruit

load("~/Dropbox/modelsF2019/output/phiFit") 

iter<-dim(MCMCchains(zmP,params=
                       "y2.sim"))[1]

ppc_dens_overlay(data$y2[1:2927], MCMCchains(zmP,params="y2.sim")[sample(iter,1000), ]) +
  theme_bw() + labs(title="Posterior predictive checks for seeds per fruit", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$y2[1:2927], MCMCchains(zmP,params="y2.sim")[sample(iter,1000), ],group=phiDat$site[1:2927]) +
  theme_bw() + labs(title="Posterior predictive checks for mean of seeds per fruit", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$y2[1:2927], MCMCchains(zmP,params="y2.sim")[sample(iter,1000), ],group=phiDat$site, stat="var") +
  theme_bw()+ labs(title="Posterior predictive checks for variance of seeds per fruit", 
                   caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

dev.off()
