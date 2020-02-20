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

## partial pooling, hyperpriors
library(gridExtra)
library(bayesplot)
library(dplyr)

color_scheme_set("brightblue")

load("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBagsPartialPoolingViabilityPartialPoolingHyperpriorsBurialFit.rds")
load("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBagIndexSiteModelData.rds")

iter<-dim(MCMCchains(zc_partialpoolhyperpriors,params="yvSim"))[1]


pdf(
  "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/products/seedBagsPartialPoolingViabilityBurialPPC.pdf",
  onefile=TRUE,
  paper="USr",
  height = 5.5, width = 10)

# Seed survival (s1)
ppc_dens_overlay(data$y_total, MCMCchains(zc_partialpoolhyperpriors,params="yTotalSim")[sample(iter,1000), ]) +
  theme_bw() + xlim(c(0,100)) + labs(title="Posterior predictive checks for seeds counted in seed bags in January", 
                                     caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$y_total, MCMCchains(zc_partialpoolhyperpriors,params="yTotalSim")[sample(iter,1000), ],group=data$site) +
  theme_bw() + labs(title="Posterior predictive checks for the mean of seeds counted in seed bags in January",
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

# variance is MUCH better in this model
ppc_stat_grouped(data$y_total, MCMCchains(zc_partialpoolhyperpriors,params="yTotalSim")[sample(iter,1000), ],group=data$site, stat="var") +
  theme_bw() + labs(title="Posterior predictive checks for the variance of seeds counted in seed bags in January",
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")


# Seed germination (g1)
ppc_dens_overlay(data$y_seedlings, MCMCchains(zc_partialpoolhyperpriors,params="ySeedlingsSim")[sample(iter,1000), ]) +
  theme_bw() + labs(title="Posterior predictive checks for germinated seeds counted in seed bags in January", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$y_seedlings, MCMCchains(zc_partialpoolhyperpriors,params="ySeedlingsSim")[sample(iter,1000), ],group=data$site) +
  theme_bw() + labs(title="Posterior predictive checks for the mean of germinated seeds counted in seed bags in January",
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

ppc_stat_grouped(data$y_seedlings, MCMCchains(zc_partialpoolhyperpriors,params="ySeedlingsSim")[sample(iter,1000), ],group=data$site, stat = "var") +
  theme_bw() + labs(title="Posterior predictive checks for the variance of germinated seeds counted in seed bags in January",
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

# Viability trials
# need to restrict data to responses that have data
# see github issue here https://github.com/stan-dev/bayesplot/issues/151

ppc_dens_overlay(data$yv[!is.na(data$yv)], MCMCchains(zc_partialpoolhyperpriors,params="yvSim")[sample(iter,1000), !is.na(data$yv)]) +
  theme_bw() + labs(title="Posterior predictive checks for viability trials", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

# ppc_stat_grouped(data$yv[!is.na(data$yv)], MCMCchains(zc_partialpoolhyperpriors,params="yvSim")[sample(iter,1000), !is.na(data$yv)],group=data$bag[!is.na(data$yv)]) +
#   theme_bw() + labs(title="Posterior predictive checks for the mean of viability trials",
#                     caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

# Germination experiment
ppc_dens_overlay(data$yg[!is.na(data$yg)], MCMCchains(zc_partialpoolhyperpriors,params="ygSim")[sample(iter,1000), !is.na(data$yg)]) +
  theme_bw() + labs(title="Posterior predictive checks for germination trials", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

# ppc_stat_grouped(data$yg[!is.na(data$yg)], MCMCchains(zc_partialpoolhyperpriors,params="ygSim")[sample(iter,1000), !is.na(data$yg)],group=data$bag[!is.na(data$yg)]) +
#   theme_bw() + labs(title="Posterior predictive checks for the mean of germination trials",
#                     caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

dev.off()

plot(data$y_total,apply(MCMCchains(zc_partialpoolhyperpriors,params="yTotalSim"),2,median))

plot(data$y_total/data$n_buried,apply(MCMCchains(zc_partialpoolhyperpriors,params="pi"),2,median))

p = apply(MCMCchains(zc_partialpoolhyperpriors,params="ps"),2,median)*apply(MCMCchains(zc_partialpoolhyperpriors,params="pi"),2,median)*((apply(MCMCchains(zc_partialpoolhyperpriors,params="viability"),2,median))^(1/3))
plot(data$y_seedlings/data$n_buried,p)

theta.i<-apply(MCMCchains(zc_partialpoolhyperpriors,params="theta.i"),2,median)
kappa.i<-apply(MCMCchains(zc_partialpoolhyperpriors,params="kappa.i"),2,median)

# pretty close to frequentist estimates 
alpha <- function(kappa,theta){
  alpha = kappa*theta
  return(alpha)
}

beta <- function(kappa,theta){
  kappa*(1-theta)
}

p=c()
for(i in 1:2){
a=alpha(kappa.i[i],theta.i[i])
b=beta(kappa.i[i],theta.i[i])
p[i]=a/(a+b)
}
p

seedBagExperiment %>%

  dplyr::mutate(p=totalJan/seedStart) %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(mean(p))

## partial pooling, site-year
library(gridExtra)
library(bayesplot)
library(dplyr)

color_scheme_set("brightblue")

load("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBagsPartialPoolingLogitSiteYearFit.rds")
load("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBagIndexSiteModelData.rds")

iter<-dim(MCMCchains(zc_partialpoollogit,params="yvSim"))[1]

data$year=seedBagExperiment$yearStart

pdf(
  "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/products/seedBagsPartialPoolingViabilityBurialPPC.pdf",
  onefile=TRUE,
  paper="USr",
  height = 5.5, width = 10)

# Seed survival (s1)
ppc_dens_overlay(data$y_total, MCMCchains(zc_partialpoollogit,params="yTotalSim")[sample(iter,1000), ]) +
  theme_bw() + xlim(c(0,100)) + labs(title="Posterior predictive checks for seeds counted in seed bags in January", 
                                     caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$y_total, MCMCchains(zc_partialpoollogit,params="yTotalSim")[sample(iter,1000), ],group=interaction(data$site,data$siteyear)) +
  theme_bw() + labs(title="Posterior predictive checks for the mean of seeds counted in seed bags in January",
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

# variance is MUCH better in this model
ppc_stat_grouped(data$y_total, MCMCchains(zc_partialpoollogit,params="yTotalSim")[sample(iter,1000), ],group=interaction(data$site,data$siteyear), stat="var") +
  theme_bw() + labs(title="Posterior predictive checks for the variance of seeds counted in seed bags in January",
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")


# Seed germination (g1)
ppc_dens_overlay(data$y_seedlings, MCMCchains(zc_partialpoollogit,params="ySeedlingsSim")[sample(iter,1000), ]) +
  theme_bw() + labs(title="Posterior predictive checks for germinated seeds counted in seed bags in January", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$y_seedlings, MCMCchains(zc_partialpoollogit,params="ySeedlingsSim")[sample(iter,1000), ],group=interaction(data$site,data$year)) +
  theme_bw() + labs(title="Posterior predictive checks for the mean of germinated seeds counted in seed bags in January",
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

ppc_stat_grouped(data$y_seedlings, MCMCchains(zc_partialpoollogit,params="ySeedlingsSim")[sample(iter,1000), ],group=interaction(data$site,data$siteyear), stat = "var") +
  theme_bw() + labs(title="Posterior predictive checks for zc_partialpoollogit variance of germinated seeds counted in seed bags in January",
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

# Viability trials
# need to restrict data to responses that have data
# see github issue here https://github.com/stan-dev/bayesplot/issues/151

ppc_dens_overlay(data$yv[!is.na(data$yv)], MCMCchains(zc_partialpoolhyperpriors,params="yvSim")[sample(iter,1000), !is.na(data$yv)]) +
  theme_bw() + labs(title="Posterior predictive checks for viability trials", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

# ppc_stat_grouped(data$yv[!is.na(data$yv)], MCMCchains(zc_partialpoolhyperpriors,params="yvSim")[sample(iter,1000), !is.na(data$yv)],group=data$bag[!is.na(data$yv)]) +
#   theme_bw() + labs(title="Posterior predictive checks for the mean of viability trials",
#                     caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

# Germination experiment
ppc_dens_overlay(data$yg[!is.na(data$yg)], MCMCchains(zc_partialpoolhyperpriors,params="ygSim")[sample(iter,1000), !is.na(data$yg)]) +
  theme_bw() + labs(title="Posterior predictive checks for germination trials", 
                    caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

# ppc_stat_grouped(data$yg[!is.na(data$yg)], MCMCchains(zc_partialpoolhyperpriors,params="ygSim")[sample(iter,1000), !is.na(data$yg)],group=data$bag[!is.na(data$yg)]) +
#   theme_bw() + labs(title="Posterior predictive checks for the mean of germination trials",
#                     caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

dev.off()

plot(data$y_total,apply(MCMCchains(zc_partialpoolhyperpriors,params="yTotalSim"),2,median))

plot(data$y_total/data$n_buried,apply(MCMCchains(zc_partialpoolhyperpriors,params="pi"),2,median))

p = apply(MCMCchains(zc_partialpoolhyperpriors,params="ps"),2,median)*apply(MCMCchains(zc_partialpoolhyperpriors,params="pi"),2,median)*((apply(MCMCchains(zc_partialpoolhyperpriors,params="viability"),2,median))^(1/3))
plot(data$y_seedlings/data$n_buried,p)

theta.i<-apply(MCMCchains(zc_partialpoolhyperpriors,params="theta.i"),2,median)
kappa.i<-apply(MCMCchains(zc_partialpoolhyperpriors,params="kappa.i"),2,median)

# pretty close to frequentist estimates 
alpha <- function(kappa,theta){
  alpha = kappa*theta
  return(alpha)
}

beta <- function(kappa,theta){
  kappa*(1-theta)
}

p=c()
for(i in 1:2){
  a=alpha(kappa.i[i],theta.i[i])
  b=beta(kappa.i[i],theta.i[i])
  p[i]=a/(a+b)
}
p

seedBagExperiment %>%
  
  dplyr::mutate(p=totalJan/seedStart) %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(mean(p))
