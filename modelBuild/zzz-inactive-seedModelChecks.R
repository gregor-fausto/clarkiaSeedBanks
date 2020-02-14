# -------------------------------------------------------------------
# Models for joint estimates of year 1 belowground rates
# Seed survival, germination, and viability
# Models use log-odds parameterization
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

set.seed(10)
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize seed bag data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# setwd and read data files
setwd("~/Dropbox/modelsF2019/seedBags")
df <- read.csv(file="seedBags.csv",header=TRUE)

# -------------------------------------------------------------------
# Clean and organize seed bag data
# -------------------------------------------------------------------
df$seedStart<-100

df$seedStart <- as.integer(df$seedStart)

df$totalJan <- ifelse(df$totalJan>100,100,df$totalJan)

df <- df %>% dplyr::rename(bag=bagNo)

## FOR NOW REMOVE MISSING DATA AND PROBLEMS
df<-subset(df,!is.na(df$totalJan))
df<-subset(df,!is.na(df$intactOct))
df<-subset(df,!(intactJan<intactOct))

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize viability trial data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# setwd and read data files
setwd("~/Dropbox/modelsF2019/viability/")
df2 <- read.csv(file="viability.csv",header=TRUE)

df2 <- df2 %>% dplyr::select(-c(germPerc,germNot,viabPerc,viabPerc2,condTest))
df2$bag <-as.integer(as.numeric(df2$bag))

## FOR NOW REMOVE MISSING DATA AND PROBLEMS
df2<-subset(df2,!is.na(df2$germStart))
df2<-subset(df2,!is.na(df2$germCount))
df2<-subset(df2,!is.na(df2$viabStart))
df2<-subset(df2,!is.na(df2$viabStain))

## data check
df2 %>% dplyr::filter(germStart - germCount - viabStart<0)
#
# -------------------------------------------------------------------
# Summarize viability trial data 
# -------------------------------------------------------------------

df2 <- df2 %>%
  dplyr::mutate(n_new = germCount+viabStart,
                y_new = germCount+viabStain) %>%
  dplyr::group_by(site,round,age,bag) %>%
  dplyr::summarise(y_new=sum(y_new),n_new=sum(n_new))

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Join datasets
# not sure if it's necessary to join but not sure how else to index
# -------------------------------------------------------------------
# -------------------------------------------------------------------
df_joined <- df %>%
  dplyr::left_join(df2,by=c("site","bag","round","age"))

# df_joined<-df_joined %>%
#   dplyr::group_by(site,round,age,bag) 

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Multiple-site model (3 year)
# Log-odds parameterization for viability, trials summed across bags
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# -------------------------------------------------------------------

## AGE 1 SEEDS
dat<-df_joined %>% 
  #dplyr::filter(round=="1") %>% 
  dplyr::filter(age==1)

# remove missing data for now
# learn how to model these
dat <- dat[!is.na(dat$n_new),]

names(dat)[1] <- "site"
dat$site.index <- as.integer(as.factor(dat$site))
dat$year.index <- as.integer(as.factor(dat$yearStart))

## AGE 2 SEEDS
datTwo<-df_joined %>% 
  #dplyr::filter(round=="1") %>% 
  dplyr::filter(age==2)

# remove missing data for now
# learn how to model these
datTwo <- datTwo[!is.na(datTwo$n_new),]

names(datTwo)[1] <- "site"
datTwo$site.index <- as.integer(as.factor(datTwo$site))
datTwo$year.index <- as.integer(as.factor(datTwo$yearStart))


# pass data to list for JAGS
data = list(
  yg = as.double(dat$seedling),
  yt = as.double(dat$totalJan),
  yo = as.double(dat$intactOct),
  yv = as.double(dat$y_new),
  n = as.double(dat$seedStart),
  nv = as.double(dat$n_new),
  N = nrow(dat),
  site = as.double(dat$site.index),
  nsites = length(unique(dat$site.index)),
  year = as.double(dat$year.index),
  nyears = length(unique(dat$year.index)),
  
  yt2 = as.double(datTwo$totalJan),
  n2 = as.double(datTwo$seedStart),
  N2 = nrow(datTwo),
  site2 = as.double(datTwo$site.index),
  nsites2 = length(unique(datTwo$site.index)),
  year2 = as.double(datTwo$year.index),
  nyears2 = length(unique(datTwo$year.index))
)

load("~/Dropbox/modelsF2019/output/seedbagfit")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Evaluate convergence
# -------------------------------------------------------------------
# -------------------------------------------------------------------

intercepts = c("alphaS1","alphaG1","alphaS2")
slopes = c("betaS1","betaG1","betaS2")
variances = c("sigmaS1","sigmaG1","sigmaS2")
viab = c("mu.b","sigma.b","alpha")
sims = c("yv.sim","yt.sim","yg.sim","yo.sim")

# summary table
# MCMCsummary(zc, n.eff = TRUE, params = c(intercepts,slopes,variances))
# not all Rhat statistics are under 1.05; need more iterations
# what is a reasonable sample size? right now in the 100s

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Posterior Predictive Checks

# notes: check group = interaction(dat$site,dat$yearStart)
# for stat density grouped plots
# -------------------------------------------------------------------
# -------------------------------------------------------------------
color_scheme_set("brightblue")

library(gridExtra)

 iter<-dim(MCMCchains(zc,params=
                        "yt.sim"))[1]


pdf(
  "~/Dropbox/modelsF2019/writing/seedPosteriorPredictiveChecks.pdf",
  onefile=TRUE,
  paper="USr",
  height = 5.5, width = 10)

# Seed survival (s1)
ppc_dens_overlay(data$yt, MCMCchains(zc,params="yt.sim")[sample(iter,1000), ]) +
  theme_bw() + xlim(c(0,100)) + labs(title="Posterior predictive checks for seeds counted in seed bags in January", 
                                     caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$yt, MCMCchains(zc,params="yt.sim")[sample(iter,1000), ],group=interaction(dat$site,dat$yearStart)) +
  theme_bw() + labs(title="Posterior predictive checks for the mean of seeds counted in seed bags in January", 
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")  

ppc_stat_grouped(data$yt, MCMCchains(zc,params="yt.sim")[sample(iter,1000), ],group=dat$site, stat="var") +
  theme_bw() + labs(title="Posterior predictive checks for the variance of seeds counted in seed bags in January", 
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

# Seed germination (g1)
ppc_dens_overlay(data$yg, MCMCchains(zc,params="yg.sim")[sample(iter,1000), ]) +
  theme_bw() + labs(title="Posterior predictive checks for germinated seeds counted in seed bags in January", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$yg, MCMCchains(zc,params="yg.sim")[sample(iter,1000), ],group=dat$site) +
  theme_bw() + labs(title="Posterior predictive checks for the mean of germinated seeds counted in seed bags in January", 
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

ppc_stat_grouped(data$yg, MCMCchains(zc,params="yg.sim")[sample(iter,1000), ],group=dat$site, stat="var") +
  theme_bw() + labs(title="Posterior predictive checks for the variance of germinated seeds counted in seed bags in January", 
                     caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

# Seed survival (s2)
ppc_dens_overlay(data$yo, MCMCchains(zc,params="yo.sim")[sample(iter,1000), ]) +
  theme_bw() + labs(title="Posterior predictive checks for seeds counted in seed bags in October", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$yo, MCMCchains(zc,params="yo.sim")[sample(iter,1000), ],group=dat$site) +
  theme_bw() + labs(title="Posterior predictive checks for mean of seeds counted in seed bags in October", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$yo, MCMCchains(zc,params="yo.sim")[sample(iter,1000), ],group=dat$site, stat="var") +
  theme_bw()+ labs(title="Posterior predictive checks for variance of seeds counted in seed bags in October", 
                     caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

# Seed survival (s3)
ppc_dens_overlay(data$yt2, MCMCchains(zc,params="yt2.sim")[sample(iter,1000), ]) +
  theme_bw() + labs(title="Posterior predictive checks for seeds counted in seed bags in January, year two", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$yt2, MCMCchains(zc,params="yt2.sim")[sample(iter,1000), ],group=datTwo$site) +
  theme_bw() + labs(title="Posterior predictive checks for mean of seeds counted in seed bags in January, year two", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$yt2, MCMCchains(zc,params="yt2.sim")[sample(iter,1000), ],group=datTwo$site, stat="var") +
  theme_bw()+ labs(title="Posterior predictive checks for variance of seeds counted in seed bags in January, year two", 
                   caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")


dev.off()
