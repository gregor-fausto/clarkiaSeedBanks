## model checking for the viability model

load("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/output/viabilityModelFit.rds")

library(MCMCvis)

MCMCsummary(zc, params = c("mu.b","sigma.b"))

load("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/output/viabilityModelData.rds")

library(ggplot2)
library(bayesplot)

iter=3000

ppc_dens_overlay(data$yv, MCMCchains(zc,params="yv.sim")[sample(iter,1000), ]) +
  theme_bw() + xlim(c(0,100)) + labs(title="Posterior predictive checks for viable seeds in viability trials", 
                                     caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(data$yv, MCMCchains(zc,params="yv.sim")[sample(iter,1000), ],group=interaction(dat$site)) +
  theme_bw() + labs(title="Posterior predictive checks for the mean of seeds counted in viability trials", 
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")  

ppc_stat_grouped(data$yv, MCMCchains(zc,params="yv.sim")[sample(iter,1000), ],group=dat$site, stat="var") +
  theme_bw() + labs(title="Posterior predictive checks for the variance of seeds counted in viability trials", 
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

plot(data$yv,apply(MCMCchains(zc,params="yv.sim")[sample(iter,1000), ],2,median),pch=16,cex=0.5)
abline(a=0,b=1)
