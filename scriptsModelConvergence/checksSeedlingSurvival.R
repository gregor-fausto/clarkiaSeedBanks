#################################################################################
################################################################################
################################################################################
# Code for figures to compare the following modeling approaches for the seedling survivorship data
# 
# Scripts by Gregor Siegmund
# fausto.siegmund@gmail.com
# last updated 04-22-2020
#################################################################################
#################################################################################
#################################################################################
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)

library(MCMCvis)
library(tidybayes)
library(tidyverse)
library(magrittr)
library(bayesplot)
library(rethinking)
library(gridExtra)


directory = "/Users/Gregor/Dropbox/dataLibrary/posteriors/"
simFiles <- paste0(directory,list.files(directory))

mcmcSamples <- readRDS(simFiles[[11]])

################################################################################
# Data directory
#################################################################################

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedlingSurvival/"
dataFiles <- paste0(directory,list.files(directory))

data <- readRDS(dataFiles[[2]])
censusSeedlingsFruitingPlants <- readRDS(dataFiles[[1]])

################################################################################
# Seedling survival
#################################################################################

color_scheme_set("brightblue")

zc <- MCMCchains(mcmcSamples,params="fruitplNumber.sim")

iter<-dim(zc)[1]

ppc_hist(data$fruitplNumber, zc[sample(iter,5), ]) + theme_bw()
# 
# # Seedling survival to fruiting
# ppc_dens_overlay(data$fruitplNumber, zc[sample(iter,100), ]) +
#   theme_bw() + labs(title="Posterior predictive checks for seedling survival", 
#                     caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$fruitplNumber, zc[sample(iter,1000), ],group=censusSeedlingsFruitingPlants$site) +
  theme_bw() + labs(title="Posterior predictive checks for mean of seedling survival", 
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")  

ppc_stat_grouped(data$fruitplNumber, zc[sample(iter,1000), ],group=censusSeedlingsFruitingPlants$site, stat="var") +
  theme_bw() + labs(title="Posterior predictive checks for variance of seedling survival", 
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")


ppc_stat_grouped(data$fruitplNumber, zc[sample(iter,1000), ],group=interaction(censusSeedlingsFruitingPlants$site,censusSeedlingsFruitingPlants$year)) +
  theme_bw() + labs(title="Posterior predictive checks for mean of seedling survival", 
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")  
