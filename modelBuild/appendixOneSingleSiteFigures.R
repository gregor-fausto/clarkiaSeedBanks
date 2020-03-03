#################################################################################
################################################################################
################################################################################
# Code for figures to compare the following modeling approaches for the seedling survivorship data
# 1. binomial likelihood
# 2. binomial likelihood with beta prior, complete pooling
# 3. binomial likelihood with beta prior, partial pooling, parameterize mode
# 4. binomial likelihood with beta prior, partial pooling, parameterize mean
# 4. binomial likelihood with logit parameterization
#
# Scripts by Gregor Siegmund
# fausto.siegmund@gmail.com
# last updated 02-29-2020
#################################################################################
#################################################################################
#################################################################################
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)

library(MCMCvis)
library(tidybayes)

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/survivorship/"
simFiles <- paste0(directory,list.files(directory))
dirFigures = c("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/")

samplesBinLinkBetaPriorComplete <- readRDS(simFiles[[5]])
samplesBinLinkBetaPriorPartialMean <- readRDS(simFiles[[7]])
samplesBinLinkBetaPriorPartialMeanHier <- readRDS(simFiles[[6]])
simData <- readRDS(simFiles[[4]])
yearNames = levels(simData$year)

nSites = length(unique(simData$site))
nYears = length(unique(simData$year))
nPlots = length(unique(simData$plot))

# MCMCsummary(samplesBinLinkBetaPriorComplete)
# MCMCsummary(samplesBinLinkBetaPriorPartialMean,params="phi")
# MCMCsummary(samplesBinLinkBetaPriorPartialMode,params="omega")

nh.p <- MCMCchains(samplesBinLinkBetaPriorComplete,params="p") 
hyear.phi<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params="phi")

hyear.p<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params="p")
hyear.theta <- MCMCchains(samplesBinLinkBetaPriorPartialMean,params="theta")

hyearpop.phi0<-MCMCchains(samplesBinLinkBetaPriorPartialMeanHier,params="phi0")
hyearpop.phi<-MCMCchains(samplesBinLinkBetaPriorPartialMeanHier,params="phi")
hyearpop.theta <- MCMCchains(samplesBinLinkBetaPriorPartialMeanHier,params="theta")

hyearpop.p <- MCMCchains(samplesBinLinkBetaPriorPartialMeanHier,params="p")
hyearpop.p_site <- MCMCchains(samplesBinLinkBetaPriorPartialMeanHier,params="p_site")

par(mfrow=c(1,1))
plot(density(hyearpop.phi0),
     xlim=c(0,1),ylim=c(0,20),
     xlab="Probability of survivorship",
     main="Posterior distribution for hierarchical model parameters")
for(i in 1:nYears){
  lines(density(hyearpop.phi[,i]),lty=1,lwd=1,col='gray')
}
points(simData$fruitingPlantNumber/simData$seedlingNumber,
       rep(0,dim(simData)[1]),
       pch=16,
       col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))


simData$year <- as.numeric(simData$year)

pdf(paste0(dirFigures,"appendix-x-hierarchyPosteriors_nh_hyear.pdf"), width=8, height=6)
par(mfrow=c(2,5))
# compare the hierarchical model to non-hierarchical
# note the posterior is narrower for the non-hierarchical
par(mfrow=c(2,5))
for(i in 1:nYears){
  # plot hierarchical population-level estimate
  plot(density(hyearpop.phi0),
       xlim=c(0,1),ylim=c(0,30), type = 'n',
       xlab="Probability of survivorship",
       main=yearNames[i])
  
  # plot hierarchical year-level estimates
  lines(density(hyear.phi[,i]),lty=1,lwd=1,col='red')
  
  # plot the non-hierarchical year-level estimates
  lines(density(nh.p[,i]),lty=1,lwd=1,col='blue')
  
  tmp = simData %>% dplyr::filter(year==i)
  points(tmp$fruitingPlantNumber/tmp$seedlingNumber,
         rep(0,dim(tmp)[1]),
         pch=16,
         col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
  
}
dev.off()
pdf(paste0(dirFigures,"appendix-x-hierarchyPosteriors_nh_hyearpop.pdf"), width=8, height=6)
par(mfrow=c(2,5))
# compare hierarchical to non-hierarchical
# note the posterior is narrower for the non-hierarchical
for(i in 1:nYears){
  # plot hierarchical population-level estimate
  plot(density(hyearpop.phi0),
       xlim=c(0,1),ylim=c(0,30),
       xlab="Probability of survivorship",
       main=yearNames[i])
  
  # plot hierarchical year-level estimates
  lines(density(hyearpop.phi[,i]),lty=1,lwd=1,col='red')
  
  lines(density(nh.p[,i]),lty=1,lwd=1,col='blue')
  
  tmp = simData %>% dplyr::filter(year==i)
  points(tmp$fruitingPlantNumber/tmp$seedlingNumber,
         rep(0,dim(tmp)[1]),
         pch=16,
         col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
}
dev.off()


nh.p.sum<-apply(nh.p,2,quantile,probs=c(.025,.5,.975))
hyear.phi.sum<-apply(hyear.phi,2,quantile,probs=c(.025,.5,.975))
hyearpop.phi.sum<-apply(hyearpop.phi,2,quantile,probs=c(.025,.5,.975))

# show summary of posteriors of means with 95% credible intervals
pdf(paste0(dirFigures,"appendix-x-hierarchyPosteriorsSummary.pdf"), width=4, height=4)
par(mfrow=c(1,1))
plot(1:10-0.25,t(nh.p.sum)[,2],pch=16,
     ylim=c(0,1),xlim=c(0,11),
     main="",
     xlab="Year",
     ylab="Probability of survivorship")
segments(x0=1:10-0.25,y0=t(nh.p.sum)[,1],
         x1=1:10-0.25,y1=t(nh.p.sum)[,3])
points(1:10,t(hyear.phi.sum)[,2],pch=17)
segments(x0=1:10,y0=t(hyear.phi.sum)[,1],
         x1=1:10,y1=t(hyear.phi.sum)[,3])
points(1:10+0.25,t(hyearpop.phi.sum)[,2],pch=18)
segments(x0=1:10+0.25,y0=t(hyearpop.phi.sum)[,1],
         x1=1:10+0.25,y1=t(hyearpop.phi.sum)[,3])
legend(x=0,y=1,c("Non-hierarchical model",
                  "Hierarchical model with year-level parameters",
                  "Hierarchical model with year- and population-level parameters"),
       pch=c(16,17,18),
       cex=0.5)
dev.off()


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

iter<-dim(MCMCchains(samplesBinLinkBetaPriorPartialMeanHier,params="fruitingPlantNumberSim"))[1]


# pdf(
#   "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/products/seedBagsCompletePoolingPPC.pdf",
#   onefile=TRUE,
#   paper="USr",
#   height = 5.5, width = 10)

# ppc_dens_overlay(simData$fruitingPlantNumber, MCMCchains(samplesBinLinkBetaPriorPartialMeanHier,params="fruitingPlantNumberSim")[sample(iter,1000), ]) +
#   theme_bw() +
#   labs(title="Posterior predictive checks for number of fruiting plants", 
#   caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")

y <- simData$fruitingPlantNumber
yrep.nh <- MCMCchains(samplesBinLinkBetaPriorComplete,params="fruitingPlantNumberSim")
group <- as.factor(simData$year)

nh.mean<-ppc_stat_grouped(y, yrep.nh,
                   group=group,facet_args=list(nrow=2)) +
  theme_bw() +
  labs(title="Posterior predictive checks for number of fruiting plants", 
       caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")

nh.var<-ppc_stat_grouped(y, yrep.nh, group=group,stat="var",facet_args=list(nrow=2)) +
  theme_bw() +
  labs(title="Posterior predictive checks for number of fruiting plants", 
       caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")


yrep.hyear = MCMCchains(samplesBinLinkBetaPriorPartialMean,params="fruitingPlantNumberSim")

hieryear.mean<-ppc_stat_grouped(y, yrep.hyear,group=group) +
  theme_bw() +
  labs(title="Posterior predictive checks for number of fruiting plants", 
       caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")

hieryear.var<-ppc_stat_grouped(y, yrep.hyear,group=group,stat='var') +
  theme_bw()
  labs(title="Posterior predictive checks for number of fruiting plants", 
       caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")


 yrep.hyearpop= MCMCchains(samplesBinLinkBetaPriorPartialMeanHier,params="fruitingPlantNumberSim")
hieryearpop.mean<-ppc_stat_grouped(y, yrep.hyearpop,group=group) +
  theme_bw() +
  labs(title="Posterior predictive checks for number of fruiting plants", 
       caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")

hieryearpop.var<-ppc_stat_grouped(y, yrep.hyearpop,group=group,stat='var') +
  theme_bw() +
  labs(title="Posterior predictive checks for number of fruiting plants", 
       caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")

l = list(nh.mean,nh.var,hieryear.mean,hieryear.var,hieryearpop.mean,hieryearpop.var)

ggsave(paste0(dirFigures,"appendix-x-ppc.pdf"), marrangeGrob(grobs = l,nrow=1,ncol=1))





par(mfrow=c(2,5))
# compare hierarchical to non-hierarchical
# for posterior of p
for(i in 1:nYears){
  # plot hierarchical population-level estimate
  plot(density(hyearpop.p_site),
       xlim=c(0,1),ylim=c(0,30),
       xlab="Probability of survivorship",
       main=yearNames[i])
  
  # plot hierarchical year-level estimates
  lines(density(hyearpop.p[,i]),lty=1,lwd=1,col='red')
  
  lines(density(hyear.p[,i]),lty=1,lwd=1,col='blue')
  
  tmp = simData %>% dplyr::filter(year==i)
  points(tmp$fruitingPlantNumber/tmp$seedlingNumber,
         rep(0,dim(tmp)[1]),
         pch=16,
         col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
}

hyear.p.sum<-apply(hyear.p,2,quantile,probs=c(.025,.5,.975))
hyearpop.p.sum<-apply(hyearpop.p,2,quantile,probs=c(.025,.5,.975))

# show summary of posteriors with 95% credible intervals
par(mfrow=c(1,1))
plot(1:9-0.25,t(hyear.p.sum)[,2],pch=1,ylim=c(0,1),xlim=c(0,10))
segments(x0=1:9-0.25,y0=t(hyear.p.sum)[,1],
         x1=1:9-0.25,y1=t(hyear.p.sum)[,3])
points(1:9+0.25,t(hyearpop.p.sum)[,2],pch=3)
segments(x0=1:9+0.25,y0=t(hyearpop.p.sum)[,1],
         x1=1:9+0.25,y1=t(hyearpop.p.sum)[,3])

yearNames = levels(simData$year)

# par(mfrow=c(2,5))
# for(i in 1:nYears){
#   plot(density(p[,i]),xlim=c(-.1,1.1),main=yearNames[i])
#   #lines(density(omega[,i]),lty='dotted')
#   lines(density(phi[,i]),lty='dashed')
# }
# 
# simData$year <- as.numeric(simData$year)
# ggplot() +
#   stat_density(data=samplesBinLinkBetaPriorComplete %>%
#                  tidybayes::spread_draws(p[year]),aes(p),
#                geom='line') +
#   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
#                  tidybayes::spread_draws(phi[year]),aes(phi),color="red",
#                geom='line') +
#   facet_wrap(~year,scales="free",nrow=2) +
#   geom_jitter(data=simData ,aes(x=fruitingPlantNumber/seedlingNumber,y=0)) +
#   theme_minimal()

# posteriors for the complete pooling model, the mean of the beta, the mode of the beta
# it shows that the complete pooling model (one p per year, black)
# has the smallest credible intervals which generally
# overlap those from the model parameterized via the mean (red)
# the models parameterized via the mean are slightly broader credible intervals

# ggplot() +
#   stat_density(data=samplesBinLinkBetaPriorComplete %>%
#                  tidybayes::spread_draws(p[year]),aes(p),geom='line') +
#   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
#                  tidybayes::spread_draws(phi[year]),aes(phi),color="red",geom='line') +
#   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
#                  tidybayes::spread_draws(p[year]),aes(p),color="red",geom='line',linetype='dotted') +
#   stat_density(data=samplesBinLinkBetaPriorPartialMode %>%
#                  tidybayes::spread_draws(omega[year]),aes(omega),color="blue",geom='line') +
#   stat_density(data=samplesBinLinkBetaPriorPartialMode %>%
#                  tidybayes::spread_draws(p[year]),aes(p),color="blue",geom='line',linetype='dotted') +
#   facet_wrap(~year,scales="free",nrow=2) +
#   geom_jitter(data=simData ,aes(x=fruitingPlantNumber/seedlingNumber,y=0)) +
#   theme_minimal()

# the posterior distribution for probability of success in year p (marginalized over all plots)
# is broader than the posterior for the mean probability of success in year p
# parameterization via the mean and mode give similar results once marginalized
# though the mode parameterization seems to give broader posteriors

simData$year <- as.numeric(as.factor(simData$year))
simData$plot <- as.numeric(as.factor(simData$plot))
# ggplot() +
#   stat_density(data=samplesBinLinkBetaPriorComplete %>%
#                  tidybayes::spread_draws(p[year]) %>%
#                  dplyr::filter(year==5),aes(p),lwd=1, geom='line') +
#   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
#                  tidybayes::spread_draws(phi[year]) %>% 
#                  dplyr::filter(year==5),aes(phi),color="red",lwd=1, geom='line') +
#   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
#                  tidybayes::spread_draws(p[year]) %>% 
#                  dplyr::filter(year==5),aes(p),
#                color="red",linetype='dashed',lwd=0.5, geom="line") +
#   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
#                  tidybayes::spread_draws(theta[year,plot]) %>% 
#                  dplyr::filter(year==5) ,aes(theta,group=plot),
#                color="black",lwd=0.25, geom="line") +
#   geom_point(data=simData  %>% 
#                dplyr::filter(year==5),
#              aes(x=fruitingPlantNumber/seedlingNumber,y=0),position=ggstance::position_dodgev(height=0.3)) +
#   theme_minimal() +
#   facet_wrap(~plot, scales = 'free')


samplesBinLinkBetaPriorPartialMeanHier %<>% recover_types(simData)
# plot 2
ggplot() +
  stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
                 tidybayes::spread_draws(phi0),aes(phi0),color="red",geom='line',size=1) +
  stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
                 tidybayes::spread_draws(phi[site,year]),aes(phi,group=year),geom='line',size=0.25, position='identity') +
  geom_point(data=simData,
             aes(x=fruitingPlantNumber/seedlingNumber,y=0),
             position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
  # facet_wrap(~year,scales="free_y",nrow=2) +
  theme_minimal()

ggplot() +
  stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
                 tidybayes::spread_draws(p_site),aes(p_site),color="red",geom='line',size=1) +
  stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
                 tidybayes::spread_draws(p[site,year]),aes(p,group=year),geom='line',size=0.25, position='identity') +
  geom_point(data=simData,
             aes(x=fruitingPlantNumber/seedlingNumber,y=0),
             position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
  facet_wrap(~year,scales="free_y",nrow=2) +
  theme_minimal()

# compare levels of the hierarchy
ggplot() +
  # posterior distribution of the mean probability of success at the site
  stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
                 tidybayes::spread_draws(phi0),aes(phi0),color="red",geom='line',size=.5) +
  # posterior distribution of the mean probability of success in each year
  stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
                 tidybayes::spread_draws(phi[site,year]),aes(phi,group=year),geom='line',size=0.25, position='identity') +
  # posterior distribution for the probability of success at the site
  stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
                 tidybayes::spread_draws(p_site),aes(p_site),color="red",geom='line',size=1) +
  # posterior distribution for the probability of success in each year 
  stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
                 tidybayes::spread_draws(p[site,year]),aes(p,group=year),geom='line',size=1, position='identity') +
  geom_point(data=simData,
             aes(x=fruitingPlantNumber/seedlingNumber,y=0),
             position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
  facet_wrap(~year,scales="free_y",nrow=2) +
  theme_minimal()

# key takeaway from the plot above: 
# the posterior for probability of success at each level (site, year-within-site) is broader than the mean probability of success at that level (compare thin to solid lines)

ggplot() +
  # posterior distribution of the mean probability of success at the site
  stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
                 tidybayes::spread_draws(phi0),aes(phi0),color="red",geom='line',size=1) +
  # posterior distribution of the mean probability of success in each year
  # stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
  #                tidybayes::spread_draws(phi[site,year]),aes(phi,group=as.factor(year)),geom='line',size=0.25, position='identity') +
  # posterior distribution for the probability of success at the site
  stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
                 tidybayes::spread_draws(p_site),aes(p_site),color="red",geom='line',size=.5) +
  # posterior distribution for the probability of success in each year 
  # stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
  #                tidybayes::spread_draws(p[site,year]),aes(p,group=year),geom='line',size=1, position='identity') +
  geom_point(data=simData,
             aes(x=fruitingPlantNumber/seedlingNumber,y=0),
             position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
  #facet_wrap(~year,scales="free_y",nrow=2) +
  theme_minimal()

ggplot() +
  # posterior distribution of the mean probability of success at the site
  # stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
  #                tidybayes::spread_draws(phi0),aes(phi0),color="red",geom='line',size=1) +
  # posterior distribution of the mean probability of success in each year
  stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
                 tidybayes::spread_draws(phi[site,year]),aes(phi,group=as.factor(year)),geom='line',size=0.25, position='identity') +
  # posterior distribution for the probability of success at the site
  # stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
  #                tidybayes::spread_draws(p_site),aes(p_site),color="red",geom='line',size=.5) +
  # posterior distribution for the probability of success in each year 
  stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
                 tidybayes::spread_draws(p[site,year]),aes(p,group=year),geom='line',size=1, position='identity') +
  geom_point(data=simData,
             aes(x=fruitingPlantNumber/seedlingNumber,y=0),
             position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
  facet_wrap(~year,scales="free_y",nrow=2) +
  theme_minimal()

# effect of hierarchy on estimate of means - estimates pool towards hierarchy with little data
ggplot() +
  # posterior distribution of the mean probability of success at the site
  stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
                 tidybayes::spread_draws(phi0),aes(phi0),color="red",geom='line',size=1) +
  # posterior distribution of the mean probability of success in each year
  stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
                 tidybayes::spread_draws(phi[site,year]),aes(phi,group=as.factor(year)),geom='line',size=0.25, position='identity',color='red') +
  # posterior distribution of the mean probability of success in each year
  stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
                 tidybayes::spread_draws(phi[year]),aes(phi,group=as.factor(year)),geom='line',size=0.25, position='identity') +
  geom_point(data=simData,
             aes(x=fruitingPlantNumber/seedlingNumber,y=0),
             position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
  facet_wrap(~year,scales="free_y",nrow=2) +
  theme_minimal()





# compare the three different fits (no hierarchy, hierarchy at 1 level, hierarchy at 2 levels)
# difference is greatest with small datasets (year = 9 only has 1 data point)
ggplot() +
  stat_density(data=samplesBinLinkBetaPriorComplete %>%
                 tidybayes::spread_draws(p[year]),aes(p),color="black",geom='line',size=1) +
  stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
                 tidybayes::spread_draws(phi[year]),aes(phi),color="red",geom='line') +
  stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
                 tidybayes::spread_draws(phi[site,year]),aes(phi),color="blue",geom='line') +
  geom_point(data=simData,
             aes(x=fruitingPlantNumber/seedlingNumber,y=0),
             position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
  facet_wrap(~year,scales="free_y",nrow=2) +
  theme_minimal()

# compare distribution of p
ggplot() +
  stat_density(data=samplesBinLinkBetaPriorComplete %>%
                 tidybayes::spread_draws(p[year]),aes(p),color="black",geom='line',size=1) +
  stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
                 tidybayes::spread_draws(p[year]),aes(p),color="red",geom='line') +
  stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
                 tidybayes::spread_draws(p[site,year]),aes(p),color="blue",geom='line') +
  geom_point(data=simData,
             aes(x=fruitingPlantNumber/seedlingNumber,y=0),
             position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
  facet_wrap(~year,scales="free_y",nrow=2) +
  theme_minimal()


nohier <-apply(MCMCchains(samplesBinLinkBetaPriorComplete,params='p'),2,median)
hier1 <- apply(MCMCchains(samplesBinLinkBetaPriorPartialMean,params='phi'),2,median)
hier2 <-apply(MCMCchains(samplesBinLinkBetaPriorPartialMeanHier,params='phi'),2,median)

plot(nohier,hier1);abline(a=0,b=1)
plot(nohier,hier2);abline(a=0,b=1)
plot(hier1,hier2);abline(a=0,b=1)

samplesBinLinkBetaPriorPartialMean %>%
  tidybayes::spread_draws(theta[year,plot])

ref<-simData %>% 
  dplyr::select(year,plot) %>%
  tidyr::unite(yearPlot,c(year,plot))

theta.med<-samplesBinLinkBetaPriorPartialMean %>%
  tidybayes::spread_draws(theta[year,plot]) %>%
  tidyr::unite(yearPlot,c(year,plot)) %>%
  dplyr::filter(yearPlot %in% ref$yearPlot) %>%
  tidyr::separate(yearPlot, c("year","plot")) %>%
  dplyr::group_by(year,plot) %>%
  dplyr::summarise(theta.med = median(theta))

theta.med$year <- as.numeric(theta.med$year)
theta.med$plot <- as.numeric(theta.med$plot)

out<-simData %>%
  dplyr::left_join(theta.med,by=c("year","plot"))

par(mfrow=c(1,1))
plot(out$fruitingPlantNumber/out$seedlingNumber,
     out$theta.med)

plot(out$seedlingNumber,out$theta.med-out$fruitingPlantNumber/out$seedlingNumber)

p.summary<-samplesBinLinkBetaPriorPartialMean %>%
  tidybayes::spread_draws(p[year]) %>%
  dplyr::summarise(p.med = median(p))
p.summary$year <- as.numeric(p.summary$year)

phi.summary<-samplesBinLinkBetaPriorPartialMean %>%
  tidybayes::spread_draws(phi[year]) %>%
  dplyr::summarise(phi.med = median(phi))
phi.summary$year <- as.numeric(phi.summary$year)

ggplot() +
  geom_point(data=out,aes(x=fruitingPlantNumber/seedlingNumber,y=theta.med)) +
  geom_abline(intercept=0,slope=1,lwd=0.5) +
  geom_hline(data=p.summary,aes(yintercept=p.med),size=0.5) +
  geom_hline(data=phi.summary,aes(yintercept=phi.med),size=0.5,color='red') +
  
  facet_wrap(~year) +
  theme_minimal()


theta.med<-samplesBinLinkBetaPriorPartialMean %>%
  tidybayes::spread_draws(p[year]) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(p.med = median(p))

theta.med$year <- as.numeric(theta.med$year)
theta.med$plot <- as.numeric(theta.med$plot)

out<-simData %>%
  dplyr::left_join(theta.med,by=c("year","plot"))

par(mfrow=c(1,1))
plot(out$fruitingPlantNumber/out$seedlingNumber,
     out$theta.med)


thetaBinLinkBetaPriorPartialMode <- MCMCchains(samplesBinLinkBetaPriorPartialMode,params="theta")

# plot(density(thetaBinLinkBetaPriorPartialMode[,1]), type='n', 
#      xlim=c(-.2,1.2),ylim=c(0,30))
# for(i in 211:241){
#   lines(density(thetaBinLinkBetaPriorPartialMode[,i]))
# }

omega<-MCMCchains(samplesBinLinkBetaPriorPartialMode,params="omega")

plot(density(omega[,1]), type='n',
     xlim=c(-.2,1.2),ylim=c(0,80))
for(i in 1:10){
  lines(density(omega[,i]))
}

kappa<-MCMCchains(samplesBinLinkBetaPriorPartialMode,params="kappa")

alpha = omega*(kappa-2)+1
beta = (1-omega)*(kappa-2)+1
plot(seq(0,1,by=0.001),
     dbeta(seq(0,1,by=0.001),median(alpha),median(beta)),
     type='n',
     ylim=c(0,20))
for(i in 1:20){
  lines(seq(0,1,by=0.001),
        dbeta(seq(0,1,by=0.001),median(alpha[i]),median(beta[i])),
        type='l',lwd=.5,col='red')
}


thetaBinLinkBetaPriorPartialMean <- MCMCchains(samplesBinLinkBetaPriorPartialMean,params="theta")

plot(density(thetaBinLinkBetaPriorPartialMean[,1]), type='n', 
     xlim=c(-.2,1.2),ylim=c(0,30))
for(i in 211:241){
  lines(density(thetaBinLinkBetaPriorPartialMean[,i]))
}

phi<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params="phi")


plot(density(phi[,1]), type='n',
     xlim=c(-.2,1.2),ylim=c(0,80))
for(i in 1:10){
  lines(density(phi[,i]))
}


kappa<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params="kappa")
alpha = kappa*phi
beta = kappa*(1-phi)
plot(seq(0,1,by=0.001),dbeta(seq(0,1,by=0.001),median(alpha),median(beta)),
     type='n')
for(i in 1:20){
  lines(seq(0,1,by=0.001),
        dbeta(seq(0,1,by=0.001),median(alpha[i]),median(beta[i])),
        type='l',lwd=.5,col='red')
}

p<-MCMCchains(samplesBinLinkBetaPriorComplete,params='p')

# compare posteriors for the hyperparameters
# the parameterization via the mode seems to sample better in general
phi<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params='phi')
kappa.mean<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params='kappa')
thin <- function(x) sample(x,1000)
par(mfrow=c(2,5))
for(i in 1:10){
  plot(thin(phi[,i]),thin(kappa.mean[,i]),cex=.5,pch=16)
}


omega<-MCMCchains(samplesBinLinkBetaPriorPartialMode,params='omega')
kappa.mode<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params='kappa')
thin <- function(x) sample(x,1000)
par(mfrow=c(2,5))
for(i in 1:10){
  plot(thin(omega[,i]),thin(kappa.mode[,i]),cex=.5,pch=16)
}

plot(apply(p,2,median),apply(phi,2,median),xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)
plot(apply(p,2,median),apply(omega,2,median),xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)
plot(apply(phi,2,median),apply(omega,2,median),xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)

stat.summary <- function(x){
  
  x %<>% recover_types(simData)
  
  tmp <- x %>% 
    tidybayes::spread_draws(p[year]) %>%
    dplyr::summarise(ci.lo = quantile(theta, prob = c(.025)),
                     med = quantile(theta, prob = c(.5)),
                     ci.hi = quantile(theta, prob = c(.975)))
  return(tmp)
}

thetaBinLinkBetaPriorPartialMode <- stat.summary(samplesBinLinkBetaPriorPartialMode)
thetaBinLinkBetaPriorPartialMean <- stat.summary(samplesBinLinkBetaPriorPartialMean)

d<-data.frame(thetaBinLinkBetaPriorPartialMode,thetaBinLinkBetaPriorPartialMean) %>% 
  dplyr::mutate(abs(med.1-med)) %>%
  dplyr::left_join(simData,by=c('year','plot'))

d<-d %>% dplyr::filter(!is.na(fruitingPlantNumber))

plot(d$med,d$med.1, 
     pch=16,cex = 0.5,xlim=c(0,1),ylim=c(0,1))
# segments(x0=thetaBinLinkBetaPriorPartialMode$med,y0=thetaBinLinkBetaPriorPartialMean$ci.lo,
#          x1=thetaBinLinkBetaPriorPartialMode$med,y1=thetaBinLinkBetaPriorPartialMean$ci.hi)
# segments(x0=thetaBinLinkBetaPriorPartialMode$ci.lo,y0=thetaBinLinkBetaPriorPartialMean$med,
#          x1=thetaBinLinkBetaPriorPartialMode$ci.hi,y1=thetaBinLinkBetaPriorPartialMean$med)
abline(a=0,b=1)

mle<-simData %>%
  dplyr::group_by(year,plot) %>%
  dplyr::summarise(y = sum(fruitingPlantNumber), n = sum(seedlingNumber)) %>%
  dplyr::mutate(pHat = y/n) %>%
  dplyr::left_join(thetaBinLinkBetaPriorPartialMode)

points(x=mle$med,y=mle$pHat,pch=2)

## summary
MCMCsummary(samplesBinLinkBetaPriorPartialMean,params=c("phi","kappa"))
MCMCsummary(samplesBinLinkBetaPriorPartialMode,params=c("omega","kappa"))

samplesBinLinkBetaPriorPartialMean %<>% recover_types(simData)

d1 <- samplesBinLinkBetaPriorPartialMean %>%
  tidybayes::spread_draws(theta[year,plot]) %>%
  dplyr::filter(year==2014)

d1.summary <- d1 %>% dplyr::group_by(plot) %>%
  dplyr::summarise(mean.med=median(theta))

samplesBinLinkBetaPriorPartialMode %<>% recover_types(simData)

d2 <- samplesBinLinkBetaPriorPartialMode %>%
  tidybayes::spread_draws(theta[year,plot]) %>%
  dplyr::filter(year==2014)

d2.summary <- d2 %>% dplyr::group_by(plot) %>%
  dplyr::summarise(mode.med=median(theta))

d <- d1.summary %>%
  dplyr::left_join(d2.summary, by= "plot")
plot(d$mode.med,d$mean.med,xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)

simData %>%
  dplyr::filter(year==2014) %>%
  dplyr::left_join(d,by="plot") %>% View
pdf(file=paste0(dirFigures,"appendix-x-mle_bayes.pdf"), width=8, height=4)

