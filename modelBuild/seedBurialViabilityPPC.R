# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Posterior Predictive Checks

# notes: check group = interaction(dat$site,dat$yearStart)
# for stat density grouped plots
# -------------------------------------------------------------------
# -------------------------------------------------------------------
library(gridExtra)
library(bayesplot)
library(dplyr)

color_scheme_set("brightblue")

load("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBagsCompletePoolingFit.rds")
load("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBagsModelData.rds")

iter<-dim(MCMCchains(zc_pool,params="yvSim"))[1]


pdf(
  "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/products/seedBagsCompletePoolingPPC.pdf",
  onefile=TRUE,
  paper="USr",
  height = 5.5, width = 10)

# Seed survival (s1)
ppc_dens_overlay(data$y_total, MCMCchains(zc_pool,params="yTotalSim")[sample(iter,1000), ]) +
  theme_bw() + xlim(c(0,100)) + labs(title="Posterior predictive checks for seeds counted in seed bags in January", 
                                     caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

# ppc_stat_grouped(data$y_total, MCMCchains(zc_pool,params="yTotalSim")[sample(iter,1000), ],group=data$bag_burial) +
#   theme_bw() + labs(title="Posterior predictive checks for the mean of seeds counted in seed bags in January", 
#                     caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")  

# Seed germination (g1)
ppc_dens_overlay(data$y_seedlings, MCMCchains(zc_pool,params="ySeedlingsSim")[sample(iter,1000), ]) +
  theme_bw() + labs(title="Posterior predictive checks for germinated seeds counted in seed bags in January", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

# ppc_stat_grouped(data$y_seedlings, MCMCchains(zc_pool,params="ySeedlingsSim")[sample(iter,1000), ],group=data$bag_burial) +
#   theme_bw() + labs(title="Posterior predictive checks for the mean of germinated seeds counted in seed bags in January", 
#                     caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

# Viability trials
# need to restrict data to responses that have data
# see github issue here https://github.com/stan-dev/bayesplot/issues/151

ppc_dens_overlay(data$yv[!is.na(data$yv)], MCMCchains(zc_pool,params="yvSim")[sample(iter,1000), !is.na(data$yv)]) +
  theme_bw() + labs(title="Posterior predictive checks for viability trials", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

# ppc_stat_grouped(data$yv[!is.na(data$yv)], MCMCchains(zc_pool,params="yvSim")[sample(iter,1000), !is.na(data$yv)],group=data$bag[!is.na(data$yv)]) +
#   theme_bw() + labs(title="Posterior predictive checks for the mean of viability trials", 
#                     caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

# Germination experiment
ppc_dens_overlay(data$yg[!is.na(data$yg)], MCMCchains(zc_pool,params="ygSim")[sample(iter,1000), !is.na(data$yg)]) +
  theme_bw() + labs(title="Posterior predictive checks for germination trials", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

# ppc_stat_grouped(data$yg[!is.na(data$yg)], MCMCchains(zc_pool,params="ygSim")[sample(iter,1000), !is.na(data$yg)],group=data$bag[!is.na(data$yg)]) +
#   theme_bw() + labs(title="Posterior predictive checks for the mean of germination trials", 
#                     caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

dev.off()


## no pool
library(gridExtra)
library(bayesplot)
library(dplyr)

color_scheme_set("brightblue")

load("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBagsNoPoolingFit.rds")
load("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBagsModelData.rds")

iter<-dim(MCMCchains(zc_nopool,params="yvSim"))[1]



pdf(
  "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/products/seedBagsNoPoolingPPC.pdf",
  onefile=TRUE,
  paper="USr",
  height = 5.5, width = 10)

# Seed survival (s1)
ppc_dens_overlay(data$y_total, MCMCchains(zc_nopool,params="yTotalSim")[sample(iter,1000), ]) +
  theme_bw() + xlim(c(0,100)) + labs(title="Posterior predictive checks for seeds counted in seed bags in January", 
                                     caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$y_total, MCMCchains(zc_nopool,params="yTotalSim")[sample(iter,1000), ],group=data$bag_burial) +
  theme_bw() + labs(title="Posterior predictive checks for the mean of seeds counted in seed bags in January",
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

# Seed germination (g1)
ppc_dens_overlay(data$y_seedlings, MCMCchains(zc_nopool,params="ySeedlingsSim")[sample(iter,1000), ]) +
  theme_bw() + labs(title="Posterior predictive checks for germinated seeds counted in seed bags in January", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$y_seedlings, MCMCchains(zc_nopool,params="ySeedlingsSim")[sample(iter,1000), ],group=data$bag_burial) +
  theme_bw() + labs(title="Posterior predictive checks for the mean of germinated seeds counted in seed bags in January",
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

# Viability trials
# need to restrict data to responses that have data
# see github issue here https://github.com/stan-dev/bayesplot/issues/151

ppc_dens_overlay(data$yv[!is.na(data$yv)], MCMCchains(zc_nopool,params="yvSim")[sample(iter,1000), !is.na(data$yv)]) +
  theme_bw() + labs(title="Posterior predictive checks for viability trials", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$yv[!is.na(data$yv)], MCMCchains(zc_nopool,params="yvSim")[sample(iter,1000), !is.na(data$yv)],group=data$bag[!is.na(data$yv)]) +
  theme_bw() + labs(title="Posterior predictive checks for the mean of viability trials",
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

# Germination experiment
ppc_dens_overlay(data$yg[!is.na(data$yg)], MCMCchains(zc_nopool,params="ygSim")[sample(iter,1000), !is.na(data$yg)]) +
  theme_bw() + labs(title="Posterior predictive checks for germination trials", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$yg[!is.na(data$yg)], MCMCchains(zc_nopool,params="ygSim")[sample(iter,1000), !is.na(data$yg)],group=data$bag[!is.na(data$yg)]) +
  theme_bw() + labs(title="Posterior predictive checks for the mean of germination trials",
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

dev.off()


## partial pooling
library(gridExtra)
library(bayesplot)
library(dplyr)

color_scheme_set("brightblue")

load("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBagsPartialPoolingFit.rds")
load("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBagsModelData.rds")

iter<-dim(MCMCchains(zc_partialpool,params="yvSim"))[1]


pdf(
  "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/products/seedBagsPartialPoolingPPC.pdf",
  onefile=TRUE,
  paper="USr",
  height = 5.5, width = 10)

# Seed survival (s1)
ppc_dens_overlay(data$y_total, MCMCchains(zc_partialpool,params="yTotalSim")[sample(iter,1000), ]) +
  theme_bw() + xlim(c(0,100)) + labs(title="Posterior predictive checks for seeds counted in seed bags in January", 
                                     caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$y_total, MCMCchains(zc_partialpool,params="yTotalSim")[sample(iter,1000), ],group=data$bag_burial) +
  theme_bw() + labs(title="Posterior predictive checks for the mean of seeds counted in seed bags in January",
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

# Seed germination (g1)
ppc_dens_overlay(data$y_seedlings, MCMCchains(zc_partialpool,params="ySeedlingsSim")[sample(iter,1000), ]) +
  theme_bw() + labs(title="Posterior predictive checks for germinated seeds counted in seed bags in January", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$y_seedlings, MCMCchains(zc_partialpool,params="ySeedlingsSim")[sample(iter,1000), ],group=data$bag_burial) +
  theme_bw() + labs(title="Posterior predictive checks for the mean of germinated seeds counted in seed bags in January",
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

# Viability trials
# need to restrict data to responses that have data
# see github issue here https://github.com/stan-dev/bayesplot/issues/151

ppc_dens_overlay(data$yv[!is.na(data$yv)], MCMCchains(zc_partialpool,params="yvSim")[sample(iter,1000), !is.na(data$yv)]) +
  theme_bw() + labs(title="Posterior predictive checks for viability trials", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$yv[!is.na(data$yv)], MCMCchains(zc_partialpool,params="yvSim")[sample(iter,1000), !is.na(data$yv)],group=data$bag[!is.na(data$yv)]) +
  theme_bw() + labs(title="Posterior predictive checks for the mean of viability trials",
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

# Germination experiment
ppc_dens_overlay(data$yg[!is.na(data$yg)], MCMCchains(zc_partialpool,params="ygSim")[sample(iter,1000), !is.na(data$yg)]) +
  theme_bw() + labs(title="Posterior predictive checks for germination trials", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$yg[!is.na(data$yg)], MCMCchains(zc_partialpool,params="ygSim")[sample(iter,1000), !is.na(data$yg)],group=data$bag[!is.na(data$yg)]) +
  theme_bw() + labs(title="Posterior predictive checks for the mean of germination trials",
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

dev.off()


## partial pooling
library(gridExtra)
library(bayesplot)
library(dplyr)

color_scheme_set("brightblue")

load("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBagsPartialPoolingViabilityPartialPoolingBurialFit.rds")
load("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBagIndexSiteModelData.rds")

iter<-dim(MCMCchains(zc_partialpool,params="yvSim"))[1]


pdf(
  "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/products/seedBagsPartialPoolingViabilityBurialPPC.pdf",
  onefile=TRUE,
  paper="USr",
  height = 5.5, width = 10)

# Seed survival (s1)
ppc_dens_overlay(data$y_total, MCMCchains(zc_partialpool,params="yTotalSim")[sample(iter,1000), ]) +
  theme_bw() + xlim(c(0,100)) + labs(title="Posterior predictive checks for seeds counted in seed bags in January", 
                                     caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$y_total, MCMCchains(zc_partialpool,params="yTotalSim")[sample(iter,1000), ],group=data$site) +
  theme_bw() + labs(title="Posterior predictive checks for the mean of seeds counted in seed bags in January",
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

# variance is very poor
ppc_stat_grouped(data$y_total, MCMCchains(zc_partialpool,params="yTotalSim")[sample(iter,1000), ],group=data$site, stat="var") +
  theme_bw() + labs(title="Posterior predictive checks for the mean of seeds counted in seed bags in January",
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")


# Seed germination (g1)
ppc_dens_overlay(data$y_seedlings, MCMCchains(zc_partialpool,params="ySeedlingsSim")[sample(iter,1000), ]) +
  theme_bw() + labs(title="Posterior predictive checks for germinated seeds counted in seed bags in January", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$y_seedlings, MCMCchains(zc_partialpool,params="ySeedlingsSim")[sample(iter,1000), ],group=data$site) +
  theme_bw() + labs(title="Posterior predictive checks for the mean of germinated seeds counted in seed bags in January",
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

ppc_stat_grouped(data$y_seedlings, MCMCchains(zc_partialpool,params="ySeedlingsSim")[sample(iter,1000), ],group=data$site, stat = "var") +
  theme_bw() + labs(title="Posterior predictive checks for the mean of germinated seeds counted in seed bags in January",
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")


# Viability trials
# need to restrict data to responses that have data
# see github issue here https://github.com/stan-dev/bayesplot/issues/151

ppc_dens_overlay(data$yv[!is.na(data$yv)], MCMCchains(zc_partialpool,params="yvSim")[sample(iter,1000), !is.na(data$yv)]) +
  theme_bw() + labs(title="Posterior predictive checks for viability trials", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$yv[!is.na(data$yv)], MCMCchains(zc_partialpool,params="yvSim")[sample(iter,1000), !is.na(data$yv)],group=data$bag[!is.na(data$yv)]) +
  theme_bw() + labs(title="Posterior predictive checks for the mean of viability trials",
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

# Germination experiment
ppc_dens_overlay(data$yg[!is.na(data$yg)], MCMCchains(zc_partialpool,params="ygSim")[sample(iter,1000), !is.na(data$yg)]) +
  theme_bw() + labs(title="Posterior predictive checks for germination trials", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$yg[!is.na(data$yg)], MCMCchains(zc_partialpool,params="ygSim")[sample(iter,1000), !is.na(data$yg)],group=data$bag[!is.na(data$yg)]) +
  theme_bw() + labs(title="Posterior predictive checks for the mean of germination trials",
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

dev.off()
