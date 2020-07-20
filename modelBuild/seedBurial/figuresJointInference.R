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

directory = "/Users/Gregor/Dropbox/dataLibrary/posteriors/"
simFiles <- paste0(directory,list.files(directory))

mcmcSamples <- readRDS(simFiles[[6]])

################################################################################
# Germination
#################################################################################

directory = "/Users/Gregor/Dropbox/dataLibrary/workflow/tidyData/"
simFiles <- paste0(directory,list.files(directory))

data <- readRDS(simFiles[[1]]) %>% dplyr::filter(siteBags=="BG"|siteBags=="BR")
siteNames = unique(data$siteBags)

################################################################################
# Posteriors
#################################################################################

theta_e <- MCMCchains(mcmcSamples,params="theta_e")
plot(density(theta_e),type='n',ylim=c(0,20))
for(i in 1:57){
lines(density(theta_e[,i]),col='gray50',add=TRUE)
}

data1<-data %>% dplyr::filter(ageBags=="1")
data2<-data %>% dplyr::filter(ageBags=="2")

df<-data.frame(cbind(data1,theta_e=apply(theta_e,2,median)))

ggplot(df) + 
  geom_point(aes(x=theta_e,intactOct/100)) +
  facet_wrap(yearBags~siteBags) +
  theme_bw()

p0_e<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p0_e[siteBags,yearBags]) %>%
  dplyr::group_by(siteBags,yearBags) %>%
  dplyr::summarise(med = median(p0_e), 
                   ci.lo = quantile(p0_e,probs=0.025), 
                   ci.hi = quantile(p0_e,probs=0.975),
                   ci.lo2 = quantile(p0_e,probs=0.25), 
                   ci.hi2 = quantile(p0_e,probs=0.75)
  ) %>% data.frame

df %>% dplyr::left_join(p0_e,by=c("siteBags","yearBags"))

ggplot() + 
  geom_point(data=df,aes(x=theta_e,intactOct/100)) +
  geom_hline(data=p0_e,aes(yintercept=med),linetype='dashed') +
  facet_wrap(yearBags~siteBags,ncol=2) +
  theme_bw()


p_e<-mcmcSamples %>%
  tidybayes::recover_types(data1) %>%
  tidybayes::spread_draws(p_e[siteBags]) %>%
  dplyr::group_by(siteBags) %>%
  dplyr::summarise(med = median(p_e), 
                   ci.lo = quantile(p_e,probs=0.025), 
                   ci.hi = quantile(p_e,probs=0.975),
                   ci.lo2 = quantile(p_e,probs=0.25), 
                   ci.hi2 = quantile(p_e,probs=0.75)
  ) %>% data.frame

ggplot() + 
  geom_point(data=df,aes(x=theta_e,intactOct/100)) +
  geom_hline(data=p0_e,aes(yintercept=med),linetype='dashed') +
  geom_hline(data=p0_e,aes(yintercept=ci.lo),linetype='dashed',alpha=0.3) +
  geom_hline(data=p0_e,aes(yintercept=ci.hi),linetype='dashed',alpha=0.3) +
  
  geom_hline(data=p_e,aes(yintercept=med),linetype='dotted') +
  geom_hline(data=p_e,aes(yintercept=ci.lo),linetype='dotted',alpha=0.3) +
  geom_hline(data=p_e,aes(yintercept=ci.hi),linetype='dotted',alpha=0.3) +
  
  facet_wrap(siteBags~yearBags,nrow=2) +
  theme_bw()

################################################################################
# Posteriors
#################################################################################

theta_4 <- MCMCchains(mcmcSamples,params="theta_4")
plot(density(theta_4),type='n',ylim=c(0,20))
for(i in 1:35){
  lines(density(theta_4[,i]),col='gray50')
}

data1<-data %>% dplyr::filter(ageBags=="1")
data2<-data %>% dplyr::filter(ageBags=="2")

df<-data.frame(cbind(data2,theta_4=apply(theta_4,2,median)))

ggplot(df) + 
  geom_point(aes(x=theta_4,totalJan/100)) +
  facet_wrap(yearBags~siteBags) +
  theme_bw()

p0_4<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p0_4[siteBags,yearBags]) %>%
  dplyr::group_by(siteBags,yearBags) %>%
  dplyr::summarise(med = median(p0_4), 
                   ci.lo = quantile(p0_4,probs=0.025), 
                   ci.hi = quantile(p0_4,probs=0.975),
                   ci.lo2 = quantile(p0_4,probs=0.25), 
                   ci.hi2 = quantile(p0_4,probs=0.75)
  ) %>% data.frame

ggplot() + 
  geom_point(data=df,aes(x=theta_4,totalJan/100)) +
  geom_hline(data=p0_4,aes(yintercept=med),linetype='dashed') +
  facet_wrap(yearBags~siteBags,ncol=2) +
  theme_bw()


p_4<-mcmcSamples %>%
  tidybayes::recover_types(data1) %>%
  tidybayes::spread_draws(p_4[siteBags]) %>%
  dplyr::group_by(siteBags) %>%
  dplyr::summarise(med = median(p_4), 
                   ci.lo = quantile(p_4,probs=0.025), 
                   ci.hi = quantile(p_4,probs=0.975),
                   ci.lo2 = quantile(p_4,probs=0.25), 
                   ci.hi2 = quantile(p_4,probs=0.75)
  ) %>% data.frame

ggplot() + 
  geom_point(data=df,aes(x=theta_4,totalJan/100)) +
  geom_hline(data=p0_4,aes(yintercept=med),linetype='dashed') +
  geom_hline(data=p0_4,aes(yintercept=ci.lo),linetype='dashed',alpha=0.3) +
  geom_hline(data=p0_4,aes(yintercept=ci.hi),linetype='dashed',alpha=0.3) +
  
  geom_hline(data=p_4,aes(yintercept=med),linetype='dotted') +
  geom_hline(data=p_4,aes(yintercept=ci.lo),linetype='dotted',alpha=0.3) +
  geom_hline(data=p_4,aes(yintercept=ci.hi),linetype='dotted',alpha=0.3) +
  
  facet_wrap(siteBags~yearBags,nrow=2) +
  theme_bw()

p0_1 <- MCMCchains(mcmcSamples,params="p0_1")
p0_2 <- MCMCchains(mcmcSamples,params="p0_2")
p0_3 <- MCMCchains(mcmcSamples,params="p0_3")
p0_e <- MCMCchains(mcmcSamples,params="p0_e")

mat <- matrix(rep(NA,dim(p0_1)[1]*dim(p0_1)[2]),nrow=dim(p0_1)[1],ncol=dim(p0_1)[2])
for(i in 1:dim(p0_1)[2]){
  mat[,i]<-p0_1[,i]*(1-p0_2[,i])*p0_3[,i]
}

par(mfrow=c(3,2))
for(i in 1:6){
hist(mat[,i],breaks=50,col='lightgray',xlim=c(0,1),ylim=c(0,2000))
hist(p0_e[,i],breaks=50,col='red',xlim=c(0,1),add=TRUE)
}

hist(mat[,1],breaks=50)

pdf("~/Downloads/testing.pdf",width=8,heigh=2)

par(mfrow=c(1,3))
i<-sample(dim(df)[1],1)
for(i in 1:35){
tmp<-MCMCchains(mcmcSamples,"p0_e") %>%
  recover_types(data) %>%
  spread_draws(p0_e[siteBags,yearBags]) %>%
  dplyr::filter(siteBags==df[i,]$siteBags&yearBags==df[i,]$yearBags) %>%
  dplyr::ungroup() %>%
  dplyr::select(p0_e) %>% 
  unlist

tmp2<-MCMCchains(mcmcSamples,"p0_4") %>%
  recover_types(data) %>%
  spread_draws(p0_4[siteBags,yearBags]) %>%
  dplyr::filter(siteBags==df[i,]$siteBags&yearBags==df[i,]$yearBags) %>%
  dplyr::ungroup() %>%
  dplyr::select(p0_4) %>% 
  unlist

hist(tmp,breaks=50,xlim=c(0,1)) 
hist(tmp2,breaks=50,col='lightgray',xlim=c(0,1))
hist(MCMCchains(mcmcSamples,"theta_4")[,i],breaks=50,xlim=c(0,1),add=TRUE)
hist(tmp*MCMCchains(mcmcSamples,"theta_4")[,i],xlim=c(0,1),
     main=paste(df[i,]$siteBags,df[i,]$yearBags));abline(v=df[i,]$totalJan[1]/100,col='red')
}
dev.off()

n=10000
mu0 <- rnorm(n,0,1)
sigma0 <- extraDistr::rhnorm(n,1)
#sigma0 <- runif(n,0,2)
#sigma <- runif(n,0,2)

mu <- boot::inv.logit(rnorm(n,mu0,sigma0))

par(mfrow=c(2,2))
hist(mu,breaks=50,col='red',freq=FALSE,ylim=c(0,4))
hist(MCMCchains(mcmcSamples,"p0_4")[,1],breaks=50,xlim=c(0,1),add=TRUE,freq=FALSE)

hist(mu,breaks=50,col='red',freq=FALSE,ylim=c(0,4))
hist(MCMCchains(mcmcSamples,"p0_4")[,2],breaks=50,xlim=c(0,1),add=TRUE,freq=FALSE)

hist(mu,breaks=50,col='red',freq=FALSE,ylim=c(0,4))
hist(MCMCchains(mcmcSamples,"p0_4")[,3],breaks=50,xlim=c(0,1),add=TRUE,freq=FALSE)

hist(mu,breaks=50,col='red',freq=FALSE,ylim=c(0,4))
hist(MCMCchains(mcmcSamples,"p0_4")[,4],breaks=50,xlim=c(0,1),add=TRUE,freq=FALSE)
################################################################################
# Raw data
#################################################################################

means<-data %>% 
  dplyr::group_by(ageBags,siteBags,yearBags) %>%
  dplyr::summarise(mu=mean(intactOct)) %>%
  dplyr::filter(ageBags!="3")

ref<-data %>% dplyr::select(siteBags,yearBags,ageBags)

out <- data %>%
  dplyr::filter(ageBags!="3")

ggplot() +
  geom_point(data=data1,aes(x=ageBags,y=intactOct)) +
  geom_point(data=data2,aes(x=ageBags,y=totalJan)) +
  #geom_point(data=means,aes(x=ageBags,y=mu),color='red') +
  facet_wrap(siteBags~yearBags,nrow=2) +
  theme_bw() 

################################################################################
# Posterior predictive
#################################################################################

iter<-dim(MCMCchains(mcmcSamples,params="y_octtot"))[1]
data1<-data %>% dplyr::filter(ageBags=="1")

y<-data1$intactOct
yrep<-as.matrix(MCMCchains(mcmcSamples,params="y_octtot")[sample(iter,100), ])

ppc_dens_overlay(y,yrep ) +
  theme_bw() +
  labs(title="Posterior predictive checks for intact seeds in October",
       caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 100 draws of the posterior")

group <-interaction( as.factor(data1$siteBags),as.factor(data1$yearBags) )

color_scheme_set("brightblue")

yrep<-as.matrix(MCMCchains(mcmcSamples,params="y_octtot")[sample(iter,1000), ])

ppc_stat_grouped(y, yrep,
                 group=group,facet_args=list(nrow=2)) +
  theme_bw() +
  labs(title="Posterior predictive checks for mean number of intact seeds in October",
       caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(y, yrep, group=group,stat="var",facet_args=list(nrow=2)) +
  theme_bw() +
  labs(title="Posterior predictive checks for variance of number of intact seeds in October",
       caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")


#
data2<-data %>% dplyr::filter(ageBags=="2")

y<-data2$totalJan
yrep<-as.matrix(MCMCchains(mcmcSamples,params="y_total2")[sample(iter,100), ])

ppc_dens_overlay(y,yrep ) +
  theme_bw() +
  labs(title="Posterior predictive checks for number of fruiting plants",
       caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")

group <-interaction( as.factor(data2$siteBags),as.factor(data2$yearBags) )

color_scheme_set("brightblue")

yrep<-as.matrix(MCMCchains(mcmcSamples,params="y_total2")[sample(iter,1000), ])

ppc_stat_grouped(y, yrep,
                 group=group,facet_args=list(nrow=2)) +
  theme_bw() +
  labs(title="Posterior predictive checks for total seeds in January",
       caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")

ppc_stat_grouped(y, yrep, group=group,stat="var",facet_args=list(nrow=2)) +
  theme_bw() +
  labs(title="Posterior predictive checks for total seeds in January",
       caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")


################################################################################
# Posteriors
#################################################################################

mcmcSummary1<-mcmcSamples1 %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(g1[siteNames]) %>%
  dplyr::summarise(med = median(g1), 
                   ci.lo = quantile(g1,probs=0.025), 
                   ci.hi = quantile(g1,probs=0.975),
                   ci.lo2 = quantile(g1,probs=0.25), 
                   ci.hi2 = quantile(g1,probs=0.75)
  )
mcmcSummary1<-cbind(mcmcSummary1,site=siteNames)

vr.site = mcmcSummary1 %>% dplyr::select(site,med)

siteIndex <- data.frame(siteIndex=unique(data$siteBags),siteBags=1:2)
yearIndex <- data.frame(yearBags=1:3,yearIndex=2006:2008)

mcmcSummary1<-mcmcSamples1 %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(g1.0[siteBags,yearBags]) %>%
  dplyr::group_by(siteBags,yearBags) %>%
  dplyr::summarise(med = median(g1.0), 
                   ci.lo = quantile(g1.0,probs=0.025), 
                   ci.hi = quantile(g1.0,probs=0.975),
                   ci.lo2 = quantile(g1.0,probs=0.25), 
                   ci.hi2 = quantile(g1.0,probs=0.75)
  )

mcmcSummary1<-mcmcSummary1 %>%
  dplyr::mutate(yearBags = as.integer(yearBags)) %>%
  #dplyr::left_join(siteIndex,by="siteBags") %>%
  dplyr::left_join(yearIndex,by="yearBags") %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(yearBags)) %>%
  #dplyr::select(-c(siteBags)) %>%
  #dplyr::rename(site = siteIndex) %>%
  dplyr::rename(site = siteBags) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)


siteIndex <- data.frame(siteIndex=unique(data$siteBags),site=1:2)
yearIndex <- data.frame(year=1:3,yearIndex=2006:2008)

interannualG1 <- mcmcSamples1 %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(g1.0[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(lambda.med = median(g1.0), 
                   ci.lo = quantile(g1.0,probs=0.27), 
                   ci.hi = quantile(g1.0,probs=0.83),
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  dplyr::select(-c(site,year)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::left_join(vr.site, by="site") %>%
  dplyr::arrange(med) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(desc(lambda.med),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  mutate(id = row_number()) %>% 
  ggplot(aes(x = id , y = lambda.med)) + 
  geom_hline(aes(yintercept=med),linetype='dotted') +
  
  geom_point(aes(color=year)) +
  
  geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +
  
  coord_flip() +
  facet_grid(site ~ ., scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ylab("Probability of germination [P(G)]") +
  ylim(c(0,1))

################################################################################
# Germination 1
#################################################################################

mcmcSummary2<-mcmcSamples2 %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(g1[siteBags]) %>%
  dplyr::group_by(siteBags) %>%
  dplyr::summarise(med = median(g1), 
                   ci.lo = quantile(g1,probs=0.025), 
                   ci.hi = quantile(g1,probs=0.975),
                   ci.lo2 = quantile(g1,probs=0.25), 
                   ci.hi2 = quantile(g1,probs=0.75)
  )
mcmcSummary2<-cbind(mcmcSummary2,site=siteNames)

vr.site = mcmcSummary2 %>% dplyr::select(site,med)

siteIndex <- data.frame(siteIndex=unique(data$siteBags),siteBags=1)
yearIndex <- data.frame(yearBags=1:3,yearIndex=2006:2008)

mcmcSummary2<-mcmcSamples2 %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(g1.0[siteBags,yearBags]) %>%
  dplyr::group_by(siteBags,yearBags) %>%
  dplyr::summarise(med = median(g1.0), 
                   ci.lo = quantile(g1.0,probs=0.025), 
                   ci.hi = quantile(g1.0,probs=0.975),
                   ci.lo2 = quantile(g1.0,probs=0.25), 
                   ci.hi2 = quantile(g1.0,probs=0.75)
  )

mcmcSummary2<-mcmcSummary2 %>%
  dplyr::mutate(yearBags = as.integer(yearBags)) %>%
  #dplyr::left_join(siteIndex,by="siteBags") %>%
  dplyr::left_join(yearIndex,by="yearBags") %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(yearBags)) %>%
  #dplyr::select(-c(siteBags)) %>%
  #dplyr::rename(site = siteIndex) %>%
  dplyr::rename(site = siteBags) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)


siteIndex <- data.frame(siteIndex=unique(data$siteBags),site=1:2)
yearIndex <- data.frame(year=1:3,yearIndex=2006:2008)

interannualG1.2 <- mcmcSamples2 %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(g1.0[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(lambda.med = median(g1.0), 
                   ci.lo = quantile(g1.0,probs=0.27), 
                   ci.hi = quantile(g1.0,probs=0.83),
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  dplyr::select(-c(site,year)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::left_join(vr.site, by="site") %>%
  dplyr::arrange(med) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(desc(lambda.med),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  mutate(id = row_number()) %>% 
  ggplot(aes(x = id , y = lambda.med)) + 
  geom_hline(aes(yintercept=med),linetype='dotted') +
  
  geom_point(aes(color=year)) +
  
  geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +
  
  coord_flip() +
  facet_grid(site ~ ., scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ylab("Probability of germination [P(G)]") +
  ylim(c(0,1))

gridExtra::grid.arrange(interannualG1,interannualG1.2)

means<-data %>% 
  dplyr::group_by(siteBags,yearBags) %>%
  dplyr::summarise(mu=mean(seedlingJan/totalJan))

data %>%
  dplyr::left_join(means,by=c("siteBags","yearBags")) %>%
  dplyr::rename(site = siteBags) %>%
  dplyr::rename(year = yearBags) %>%
  dplyr::left_join(vr.site, by="site") %>%
  dplyr::arrange(med) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  # dplyr::group_by(site) %>%
  #dplyr::arrange(desc(lambda.med),.by_group = TRUE) +
  ggplot(aes(x=year,y=seedlingJan/totalJan)) +
  geom_point() +
  geom_point(aes(x=year,y=mu),color='red') +
  coord_flip() +
  facet_grid(site ~ ., scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(
    axis.ticks.y=element_blank()) +
  ylab("Probability of germination [P(G)]") +
  ylim(c(0,1))

gridExtra::grid.arrange(interannualG1,interannualG1.2)

means<-data %>% 
  # dplyr::left_join(ref,by=c("siteBags","yearBags")) %>%
  dplyr::group_by(ageBags,siteBags,yearBags) %>%
  dplyr::summarise(mu=mean(seedlingJan)) %>%
  dplyr::filter(ageBags!="3")

ref<-data %>% dplyr::select(siteBags,yearBags,ageBags)

out <- data %>%
  dplyr::filter(ageBags!="3")

ggplot() +
  geom_point(data=out,aes(x=ageBags,y=seedlingJan)) +
  geom_point(data=means,aes(x=ageBags,y=mu),color='red') +
  facet_wrap(yearBags~siteBags) +
  theme_bw() 

v1<-MCMCsummary(mcmcSamples1,params=c("s1.0","g1.0","s2.0"))$`50%`
v2<-MCMCsummary(mcmcSamples2,params=c("s1.0","g1.0","s2.0"))$`50%`
names <- rownames(MCMCsummary(mcmcSamples1,c("s1.0","g1.0","s2.0")))
d<-data.frame(names,v1,v2)
ggplot(d,aes(x=v1,y=v2,color=names)) + geom_abline(intercept=0,slope=1)+ geom_point()

################################################################################
# Winter seed survival (s1)
#################################################################################
# 
# mcmcSummary<-mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(s1[siteBags]) %>%
#   dplyr::group_by(siteBags) %>%
#   dplyr::summarise(med = median(s1), 
#                    ci.lo = quantile(s1,probs=0.025), 
#                    ci.hi = quantile(s1,probs=0.975),
#                    ci.lo2 = quantile(s1,probs=0.25), 
#                    ci.hi2 = quantile(s1,probs=0.75)
#   )
# mcmcSummary<-cbind(mcmcSummary,site=siteNames)
# 
# 
# s1Summary <- mcmcSummary %>%
#   dplyr::rename(ci.lo95 = ci.lo,
#                 ci.hi95 = ci.hi,
#                 ci.lo50 = ci.lo2,
#                 ci.hi50 = ci.hi2)
# 
# fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/parameterSummary/")
# dir.create(file.path(fileDirectory), showWarnings = FALSE)
# saveRDS(s1Summary,file=paste0(fileDirectory,"s1Summary.RDS"))
# 
# 
# vr = position %>% 
#   dplyr::left_join(mcmcSummary,by="site")
# 
# vr.site = vr %>%
#   dplyr::select(site,med)
# 
# pdf(paste0(dirFigures,"spatial-s1.pdf"), width=8, height=6)
# 
# par(mfrow=c(1,1))
# plot(vr$easting,vr$med,ylim=c(0.0,1),
#      ylab="Seed survival probability in first winter [P(S1)]",
#      xlab="Easting (km)", 
#      cex.lab = 1.5, cex.axis = 1.5,
#      pch=16,type='n')
# 
# points(vr$easting,vr$med,ylim=c(0.0,1),
#        pch=16)
# 
# segments(x0=vr$easting, y0=vr$ci.lo, y1=vr$ci.hi)
# segments(x0=vr$easting, y0=vr$ci.lo2, y1=vr$ci.hi2,lwd=2)
# 
# dev.off()
# 
# s1<-MCMCchains(mcmcSamples,params="s1")
# 
# plot(density(s1),type='n',ylim=c(0,8))
# for(i in 1:20){
#   lines(density(s1[,i]),lwd=.5)
# }
# 
# 
# 
# siteIndex <- data.frame(siteIndex=unique(data$siteBags),siteBags=1:20)
# yearIndex <- data.frame(yearBags=1:3,yearIndex=2006:2008)
# 
# mcmcSummary<-mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(s1.0[siteBags,yearBags]) %>%
#   dplyr::group_by(siteBags,yearBags) %>%
#   dplyr::summarise(med = median(s1.0), 
#                    ci.lo = quantile(s1.0,probs=0.025), 
#                    ci.hi = quantile(s1.0,probs=0.975),
#                    ci.lo2 = quantile(s1.0,probs=0.25), 
#                    ci.hi2 = quantile(s1.0,probs=0.75)
#   )
# 
# mcmcSummary<-mcmcSummary %>%
#   dplyr::mutate(yearBags = as.integer(yearBags)) %>%
#   #dplyr::left_join(siteIndex,by="siteBags") %>%
#   dplyr::left_join(yearIndex,by="yearBags") %>%
#   dplyr::ungroup() %>%
#   dplyr::select(-c(yearBags)) %>%
#   #dplyr::select(-c(siteBags)) %>%
#   #dplyr::rename(site = siteIndex) %>%
#   dplyr::rename(site = siteBags) %>%
#   dplyr::rename(year = yearIndex) %>%
#   dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)
# 
# 
# vr = position %>% 
#   dplyr::left_join(mcmcSummary,by="site")
# 
# 
# 
# 
# tmp <- vr %>% 
#   dplyr::select(site,year,med,ci.lo,ci.hi)
# 
# s1Summary <- tmp
# write.csv(s1Summary , 
#           file = "~/Dropbox/clarkiaSeedBanks/products/dataFiles2/s1Summary.csv",
#           row.names=FALSE)
# 
# g1 <- ggplot(data=vr) +
#   geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +
#   geom_linerange(aes(x=easting,ymin=ci.lo2,ymax=ci.hi2),size=1) +
#   geom_point(aes(x=easting,y=med)) +
#   facet_wrap(~year,nrow=1) +
#   theme_bw() +
#   ylab("Seed survival probability in first winter [P(S1)]") +
#   xlab("Easting (km)") 
# 
# 
# ggsave(filename=paste0(dirFigures,"annual-s1.pdf"),
#        plot=g1,width=20,height=4)
# 
# library(pdftools)
# pdf_combine(c(paste0(dirFigures,"spatial-s1.pdf"),
#               paste0(dirFigures,"annual-s1.pdf")), 
#             output = paste0(dirFigures,"summary-s1.pdf"))
# 
# siteIndex <- data.frame(siteIndex=unique(data$siteBags),site=1:20)
# yearIndex <- data.frame(year=1:3,yearIndex=2006:2008)
# 
# interannualS1 <- mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(s1.0[site,year]) %>%
#   dplyr::group_by(site,year) %>%
#   dplyr::summarise(lambda.med = median(s1.0), 
#                    ci.lo = quantile(s1.0,probs=0.27), 
#                    ci.hi = quantile(s1.0,probs=0.83),
#   ) %>%
#   dplyr::ungroup() %>%
#   dplyr::left_join(siteIndex,by="site") %>%
#   dplyr::left_join(yearIndex,by="year") %>%
#   dplyr::select(-c(site,year)) %>%
#   dplyr::rename(site = siteIndex) %>%
#   dplyr::rename(year = yearIndex) %>%
#   dplyr::left_join(vr.site, by="site") %>%
#   dplyr::arrange(med) %>% 
#   dplyr::mutate(site=factor(site,levels=unique(site))) %>%
#   dplyr::group_by(site) %>%
#   dplyr::arrange(desc(lambda.med),.by_group = TRUE) %>%
#   dplyr::mutate(year=factor(year)) %>%
#   mutate(id = row_number()) %>% 
#   ggplot(aes(x = id , y = lambda.med)) + 
#   geom_hline(aes(yintercept=med),linetype='dotted') +
#   
#   geom_point(aes(color=year)) +
#   
#   geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +
#   
#   coord_flip() +
#   facet_grid(site ~ ., scales="free_x", space="free_x") +
#   theme_bw() +
#   theme(panel.spacing=unit(0,"pt"), 
#         panel.border=element_rect(colour="grey50", fill=NA)) +
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank()) +
#   ylab("Probability of winter seed survival [P(S1)]")
# 
# ggsave(filename=paste0(dirFigures,"interannual-S1.pdf"),
#        plot=interannualS1,width=6,height=12)
# 
# ################################################################################
# # Summer seed survival (s2)
# #################################################################################
# 
# mcmcSummary<-mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(s2[siteBags]) %>%
#   dplyr::group_by(siteBags) %>%
#   dplyr::summarise(med = median(s2), 
#                    ci.lo = quantile(s2,probs=0.025), 
#                    ci.hi = quantile(s2,probs=0.975),
#                    ci.lo2 = quantile(s2,probs=0.25), 
#                    ci.hi2 = quantile(s2,probs=0.75)
#   )
# mcmcSummary<-cbind(mcmcSummary,site=siteNames)
# 
# 
# s2Summary <- mcmcSummary %>%
#   dplyr::rename(ci.lo95 = ci.lo,
#                 ci.hi95 = ci.hi,
#                 ci.lo50 = ci.lo2,
#                 ci.hi50 = ci.hi2)
# 
# fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/parameterSummary/")
# dir.create(file.path(fileDirectory), showWarnings = FALSE)
# saveRDS(s2Summary,file=paste0(fileDirectory,"s2Summary.RDS"))
# 
# vr = position %>% 
#   dplyr::left_join(mcmcSummary,by="site")
# 
# vr.site = vr %>%
#   dplyr::select(site,med)
# 
# pdf(paste0(dirFigures,"spatial-s2.pdf"), width=8, height=6)
# 
# par(mfrow=c(1,1))
# plot(vr$easting,vr$med,ylim=c(0.0,1),
#      ylab="Seed survival probability in first summer [P(S2)]",
#      xlab="Easting (km)", 
#      cex.lab = 1.5, cex.axis = 1.5,
#      pch=16,type='n')
# 
# points(vr$easting,vr$med,ylim=c(0.0,1),
#        pch=16)
# 
# segments(x0=vr$easting, y0=vr$ci.lo, y1=vr$ci.hi)
# segments(x0=vr$easting, y0=vr$ci.lo2, y1=vr$ci.hi2,lwd=2)
# 
# dev.off()
# 
# s2<-MCMCchains(mcmcSamples,params="s2")
# 
# plot(density(s2),type='n',ylim=c(0,8))
# for(i in 1:20){
#   lines(density(s2[,i]),lwd=.5)
# }
# 
# 
# siteIndex <- data.frame(siteIndex=unique(data$siteBags),siteBags=1:20)
# yearIndex <- data.frame(yearBags=1:3,yearIndex=2006:2008)
# 
# mcmcSummary<-mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(s2.0[siteBags,yearBags]) %>%
#   dplyr::group_by(siteBags,yearBags) %>%
#   dplyr::summarise(med = median(s2.0), 
#                    ci.lo = quantile(s2.0,probs=0.025), 
#                    ci.hi = quantile(s2.0,probs=0.975),
#                    ci.lo2 = quantile(s2.0,probs=0.25), 
#                    ci.hi2 = quantile(s2.0,probs=0.75)
#   )
# 
# mcmcSummary<-mcmcSummary %>%
#   dplyr::mutate(yearBags = as.integer(yearBags)) %>%
#   #dplyr::left_join(siteIndex,by="siteBags") %>%
#   dplyr::left_join(yearIndex,by="yearBags") %>%
#   dplyr::ungroup() %>%
#   dplyr::select(-c(yearBags)) %>%
#   #dplyr::select(-c(siteBags)) %>%
#   #dplyr::rename(site = siteIndex) %>%
#   dplyr::rename(site = siteBags) %>%
#   dplyr::rename(year = yearIndex) %>%
#   dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)
# 
# 
# vr = position %>% 
#   dplyr::left_join(mcmcSummary,by="site")
# 
# 
# tmp <- vr %>% 
#   dplyr::select(site,year,med,ci.lo,ci.hi)
# 
# s2Summary <- tmp
# write.csv(s1Summary , 
#           file = "~/Dropbox/clarkiaSeedBanks/products/dataFiles2/s2Summary.csv",
#           row.names=FALSE)
# 
# g1 <- ggplot(data=vr) +
#   geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +
#   geom_linerange(aes(x=easting,ymin=ci.lo2,ymax=ci.hi2),size=1) +
#   geom_point(aes(x=easting,y=med)) +
#   facet_wrap(~year,nrow=1) +
#   theme_bw() +
#   ylab("Seed survival probability in first summer [P(S2)]") +
#   xlab("Easting (km)") 
# 
# 
# ggsave(filename=paste0(dirFigures,"annual-s2.pdf"),
#        plot=g1,width=20,height=4)
# 
# library(pdftools)
# pdf_combine(c(paste0(dirFigures,"spatial-s2.pdf"),
#               paste0(dirFigures,"annual-s2.pdf")), 
#             output = paste0(dirFigures,"summary-s2.pdf"))
# 
# siteIndex <- data.frame(siteIndex=unique(data$siteBags),site=1:20)
# yearIndex <- data.frame(year=1:3,yearIndex=2006:2008)
# 
# interannualS2 <- mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(s2.0[site,year]) %>%
#   dplyr::group_by(site,year) %>%
#   dplyr::summarise(lambda.med = median(s2.0), 
#                    ci.lo = quantile(s2.0,probs=0.27), 
#                    ci.hi = quantile(s2.0,probs=0.83),
#   ) %>%
#   dplyr::ungroup() %>%
#   dplyr::left_join(siteIndex,by="site") %>%
#   dplyr::left_join(yearIndex,by="year") %>%
#   dplyr::select(-c(site,year)) %>%
#   dplyr::rename(site = siteIndex) %>%
#   dplyr::rename(year = yearIndex) %>%
#   dplyr::left_join(vr.site, by="site") %>%
#   dplyr::arrange(med) %>% 
#   dplyr::mutate(site=factor(site,levels=unique(site))) %>%
#   dplyr::group_by(site) %>%
#   dplyr::arrange(desc(lambda.med),.by_group = TRUE) %>%
#   dplyr::mutate(year=factor(year)) %>%
#   mutate(id = row_number()) %>% 
#   ggplot(aes(x = id , y = lambda.med)) + 
#   geom_hline(aes(yintercept=med),linetype='dotted') +
#   
#   geom_point(aes(color=year)) +
#   
#   geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +
#   
#   coord_flip() +
#   facet_grid(site ~ ., scales="free_x", space="free_x") +
#   theme_bw() +
#   theme(panel.spacing=unit(0,"pt"), 
#         panel.border=element_rect(colour="grey50", fill=NA)) +
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank()) +
#   ylab("Probability of summmer seed survival [P(S2)]")
# 
# ggsave(filename=paste0(dirFigures,"interannual-S2.pdf"),
#        plot=interannualS2,width=6,height=12)
# 
# ################################################################################
# # Winter seed survival, year 2 (s3)
# #################################################################################

mcmcSummary<-mcmcSamples2 %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(s3[siteBags]) %>%
  dplyr::group_by(siteBags) %>%
  dplyr::summarise(med = median(s3),
                   ci.lo = quantile(s3,probs=0.025),
                   ci.hi = quantile(s3,probs=0.975),
                   ci.lo2 = quantile(s3,probs=0.25),
                   ci.hi2 = quantile(s3,probs=0.75)
  )
mcmcSummary<-cbind(mcmcSummary,site=siteNames)


s3Summary <- mcmcSummary %>%
  dplyr::rename(ci.lo95 = ci.lo,
                ci.hi95 = ci.hi,
                ci.lo50 = ci.lo2,
                ci.hi50 = ci.hi2)

# fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/parameterSummary/")
# dir.create(file.path(fileDirectory), showWarnings = FALSE)
# saveRDS(s3Summary,file=paste0(fileDirectory,"s3Summary.RDS"))

vr = position %>%
  dplyr::left_join(mcmcSummary,by="site")

vr.site = vr %>%
  dplyr::select(site,med)

#pdf(paste0(dirFigures,"spatial-s3.pdf"), width=8, height=6)

par(mfrow=c(1,1))
plot(vr$easting,vr$med,ylim=c(0.0,1),
     ylab="Seed survival probability in second winter [P(S3)]",
     xlab="Easting (km)",
     cex.lab = 1.5, cex.axis = 1.5,
     pch=16,type='n')

points(vr$easting,vr$med,ylim=c(0.0,1),
       pch=16)

segments(x0=vr$easting, y0=vr$ci.lo, y1=vr$ci.hi)
segments(x0=vr$easting, y0=vr$ci.lo2, y1=vr$ci.hi2,lwd=2)

#dev.off()

s3<-MCMCchains(mcmcSamples,params="s3")

plot(density(s3),type='n',ylim=c(0,8))
for(i in 1:20){
  lines(density(s3[,i]),lwd=.5)
}


siteIndex <- data.frame(siteIndex=unique(data$siteBags),siteBags=1:20)
yearIndex <- data.frame(yearBags=1:3,yearIndex=2006:2008)

mcmcSummary<-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(s3.0[siteBags,yearBags]) %>%
  dplyr::group_by(siteBags,yearBags) %>%
  dplyr::summarise(med = median(s3.0),
                   ci.lo = quantile(s3.0,probs=0.025),
                   ci.hi = quantile(s3.0,probs=0.975),
                   ci.lo2 = quantile(s3.0,probs=0.25),
                   ci.hi2 = quantile(s3.0,probs=0.75)
  )

mcmcSummary<-mcmcSummary %>%
  dplyr::mutate(yearBags = as.integer(yearBags)) %>%
  #dplyr::left_join(siteIndex,by="siteBags") %>%
  dplyr::left_join(yearIndex,by="yearBags") %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(yearBags)) %>%
  #dplyr::select(-c(siteBags)) %>%
  #dplyr::rename(site = siteIndex) %>%
  dplyr::rename(site = siteBags) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)


vr = position %>%
  dplyr::left_join(mcmcSummary,by="site")


tmp <- vr %>%
  dplyr::select(site,year,med,ci.lo,ci.hi)

# s3Summary <- tmp
# write.csv(s1Summary ,
#           file = "~/Dropbox/clarkiaSeedBanks/products/dataFiles2/s3Summary.csv",
#           row.names=FALSE)

g1 <- ggplot(data=vr) +
  geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +
  geom_linerange(aes(x=easting,ymin=ci.lo2,ymax=ci.hi2),size=1) +
  geom_point(aes(x=easting,y=med)) +
  facet_wrap(~year,nrow=1) +
  theme_bw() +
  ylab("Seed survival probability in second winter [P(S3)]") +
  xlab("Easting (km)")


# ggsave(filename=paste0(dirFigures,"annual-s3.pdf"),
#        plot=g1,width=20,height=4)

library(pdftools)
pdf_combine(c(paste0(dirFigures,"spatial-s3.pdf"),
              paste0(dirFigures,"annual-s3.pdf")),
            output = paste0(dirFigures,"summary-s3.pdf"))

siteIndex <- data.frame(siteIndex=unique(data$siteBags),site=1:20)
yearIndex <- data.frame(year=1:2,yearIndex=2006:2007)

interannualS3 <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(s3.0[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(lambda.med = median(s3.0),
                   ci.lo = quantile(s3.0,probs=0.27),
                   ci.hi = quantile(s3.0,probs=0.83),
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  dplyr::select(-c(site,year)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::left_join(vr.site, by="site") %>%
  dplyr::arrange(med) %>%
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(desc(lambda.med),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  mutate(id = row_number()) %>%
  ggplot(aes(x = id , y = lambda.med)) +
  geom_hline(aes(yintercept=med),linetype='dotted') +
  
  geom_point(aes(color=year)) +
  
  geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +
  
  coord_flip() +
  facet_grid(site ~ ., scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"),
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ylab("Probability of winter seed survival [P(S3)]")
# 
# ggsave(filename=paste0(dirFigures,"interannual-S3.pdf"),
#        plot=interannualS3,width=6,height=12)


# ################################################################################
# # ???
# #################################################################################
# ref<-seedBagExperiment %>% 
#   dplyr::select(yearBags,siteBags) %>%
#   tidyr::unite(yearSite,c(yearBags,siteBags))
# 
# nYears = length(unique(seedBagExperiment$yearBags))
# 
# sampleYear <-sample(1:nYears,1)
# par(mfrow=c(1,3))
# for(i in sampleYear){
#   # filter data to year
#   tmp = seedBagExperiment %>% 
#     dplyr::filter(yearBags==unique(seedBagExperiment$yearBags)[i]) %>%
#     dplyr::mutate(observed.p = totalJan/seedStart) %>%
#     dplyr::arrange(observed.p)
#   
#   mcmcSummary<-mcmcSamples %>%
#     tidybayes::spread_draws(theta_1[yearBags]) %>%
#     dplyr::group_by(yearBags) %>%
#     dplyr::summarise(med = median(theta_1), 
#                      ci.lo = quantile(theta_1,probs=0.025), 
#                      ci.hi = quantile(theta_1,probs=0.975))
#   
# tmp<-tmp %>%
#   dplyr::left_join(mcmcSummary,by="yearBags")
# 
# mle = sum(tmp$fruitingPlantNumber)/sum(tmp$seedlingNumber)
# 
# # plot hierarchical population-level estimate
# 
# # toggle on jitter to avoid overplotting
# jitter <- rnorm(length(tmp$observed.p),mean=0,sd=0.01)
# 
# plot(tmp$observed.p + jitter,tmp$med,
#      pch=16, cex = 1,
#      col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#      xlim=c(0,1),ylim=c(0,1),
#      xlab="Observed proportion, y[i]/n[i]",
#      ylab="theta[i] and 95% credible interval",
#      main=yearNames[i],type='n')
# 
# segments(x0=tmp$observed.p + jitter,x1=tmp$observed.p + jitter,
#          y0=tmp$ci.lo,y1=tmp$ci.hi,lwd=0.75,
#          col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5))
# points(tmp$observed.p + jitter,tmp$med,
#        pch=16, cex = 1,
#        col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5))
# 
# # one-to-one line
# abline(a=0,b=1)
# 
# # add mle estimate for the year
# abline(h = mle, lty='dashed')
# }
# 
# # compare hierarchical to non-hierarchical
# # for posterior of p
# for(i in sampleYear){
#   
#   # filter data to year
#   tmp = simData %>% 
#     dplyr::filter(year==i) %>%
#     dplyr::mutate(observed.p = fruitingPlantNumber/seedlingNumber) %>%
#     dplyr::arrange(observed.p) %>%
#     tidyr::unite(yearPlot,c(year,plot))
#   
#   mcmcSummary<-samplesBinLinkBetaPriorPartialMean %>%
#     tidybayes::spread_draws(theta[year,plot]) %>%
#     tidyr::unite(yearPlot,c(year,plot)) %>%
#     dplyr::filter(yearPlot %in% tmp$yearPlot) %>%
#     dplyr::group_by(yearPlot) %>%
#     dplyr::summarise(med = median(theta), 
#                      ci.lo = quantile(theta,probs=0.025), 
#                      ci.hi = quantile(theta,probs=0.975))
#   
#   tmp<-tmp %>%
#     dplyr::left_join(mcmcSummary,by="yearPlot")
#   
#   mle = sum(tmp$fruitingPlantNumber)/sum(tmp$seedlingNumber)
#   
#   # plot hierarchical population-level estimate
#   
#   # toggle on jitter to avoid overplotting
#   jitter <- rnorm(length(tmp$observed.p),mean=0,sd=0.01)
#   
#   plot(tmp$observed.p + jitter,tmp$med,
#        pch=16, cex = 1,
#        col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5),
#        xlim=c(0,1),ylim=c(0,1),
#        xlab="Observed proportion, y[i]/n[i]",
#        ylab="theta[i] and 95% credible interval",
#        main=yearNames[i],type='n')
#   
#   segments(x0=tmp$observed.p + jitter,x1=tmp$observed.p + jitter,
#            y0=tmp$ci.lo,y1=tmp$ci.hi,lwd=0.75,
#            col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5))
#   points(tmp$observed.p + jitter,tmp$med,
#          pch=16, cex = 1,
#          col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5))
#   
#   # one-to-one line
#   abline(a=0,b=1)
#   
#   # add mle estimate for the year
#   abline(h = mle, lty='dashed')
#   
# }
# 
# # compare hierarchical to non-hierarchical
# # for posterior of p
# for(i in sampleYear){
#   
#   # filter data to year
#   tmp = simData %>% 
#     dplyr::filter(year==i) %>%
#     dplyr::mutate(observed.p = fruitingPlantNumber/seedlingNumber) %>%
#     dplyr::arrange(observed.p) %>%
#     tidyr::unite(yearPlot,c(year,plot))
#   
#   mcmcSummary<-samplesBinLinkBetaPriorPartialMeanHier %>%
#     tidybayes::spread_draws(theta[site,year,plot]) %>%
#     tidyr::unite(yearPlot,c(year,plot)) %>%
#     dplyr::filter(yearPlot %in% tmp$yearPlot) %>%
#     dplyr::group_by(yearPlot) %>%
#     dplyr::summarise(med = median(theta), 
#                      ci.lo = quantile(theta,probs=0.025), 
#                      ci.hi = quantile(theta,probs=0.975))
#   
#   tmp<-tmp %>%
#     dplyr::left_join(mcmcSummary,by="yearPlot")
#   
#   mle = sum(tmp$fruitingPlantNumber)/sum(tmp$seedlingNumber)
#   
#   # plot hierarchical population-level estimate
#   
#   # toggle on jitter to avoid overplotting
#   jitter <- rnorm(length(tmp$observed.p),mean=0,sd=0.01)
#   
#   plot(tmp$observed.p + jitter,tmp$med,
#        pch=16, cex = 1,
#        col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        xlim=c(0,1),ylim=c(0,1),
#        xlab="Observed proportion, y[i]/n[i]",
#        ylab="theta[i] and 95% credible interval",
#        main=yearNames[i],type='n')
#   
#   segments(x0=tmp$observed.p + jitter,x1=tmp$observed.p + jitter,
#            y0=tmp$ci.lo,y1=tmp$ci.hi,lwd=0.75,
#            col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
#   points(tmp$observed.p + jitter,tmp$med,
#          pch=16, cex = 1,
#          col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
#   
#   # one-to-one line
#   abline(a=0,b=1)
#   
#   # add mle estimate for the year
#   abline(h = mle, lty='dashed')
#   
# }
# 
# dev.off()
# 
# 
# pdf(file=paste0(dirFigures,"appendix-x-figure54.pdf"), width=8, height=6)
# # compare both hierarchical models
# # purple is year-level parameters
# # gray is year- and population-level parameters
# par(mfrow=c(2,5))
# # compare hierarchical to non-hierarchical
# # for posterior of p
# for(i in 1:nYears){
#   
#   # filter data to year
#   tmp = simData %>% 
#     dplyr::filter(year==i) %>%
#     dplyr::mutate(observed.p = fruitingPlantNumber/seedlingNumber) %>%
#     dplyr::arrange(observed.p) %>%
#     tidyr::unite(yearPlot,c(year,plot))
#   
#   mcmcSummary<-samplesBinLinkBetaPriorPartialMeanHier %>%
#     tidybayes::spread_draws(theta[site,year,plot]) %>%
#     tidyr::unite(yearPlot,c(year,plot)) %>%
#     dplyr::filter(yearPlot %in% tmp$yearPlot) %>%
#     dplyr::group_by(yearPlot) %>%
#     dplyr::summarise(med.h = median(theta), 
#                      ci.lo.h = quantile(theta,probs=0.025), 
#                      ci.hi.h = quantile(theta,probs=0.975))
#   
#   mcmcSummary2<-samplesBinLinkBetaPriorPartialMean %>%
#     tidybayes::spread_draws(theta[year,plot]) %>%
#     tidyr::unite(yearPlot,c(year,plot)) %>%
#     dplyr::filter(yearPlot %in% tmp$yearPlot) %>%
#     dplyr::group_by(yearPlot) %>%
#     dplyr::summarise(med = median(theta), 
#                      ci.lo = quantile(theta,probs=0.025), 
#                      ci.hi = quantile(theta,probs=0.975))
#   
#   tmp<-tmp %>%
#     dplyr::left_join(mcmcSummary,by="yearPlot") %>%
#     dplyr::left_join(mcmcSummary2,by="yearPlot")
#   
#   mle = sum(tmp$fruitingPlantNumber)/sum(tmp$seedlingNumber)
#   
#   # plot hierarchical population-level estimate
#   
#   # toggle on jitter to avoid overplotting
#   jitter <- rnorm(length(tmp$observed.p),mean=0,sd=0.01)
#   
#   plot(tmp$observed.p + jitter,tmp$med.h,
#        pch=16, cex = 1,
#        col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
#        xlim=c(0,1),ylim=c(0,1),
#        xlab="Observed proportion, y[i]/n[i]",
#        ylab="theta[i] and 95% credible interval",
#        main=yearNames[i],type='n')
#   
#   segments(x0=tmp$observed.p + jitter,x1=tmp$observed.p + jitter,
#            y0=tmp$ci.lo.h,y1=tmp$ci.hi.h,lwd=0.75,
#            col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
#   points(tmp$observed.p + jitter,tmp$med.h,
#          pch=16, cex = 1,
#          col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
#  
#   # toggle on jitter to avoid overplotting
#   jitter <- rnorm(length(tmp$observed.p),mean=0,sd=0.01)
#   
#   # lower hierarchy
#   segments(x0=tmp$observed.p + jitter,x1=tmp$observed.p + jitter,
#            y0=tmp$ci.lo,y1=tmp$ci.hi,lwd=0.75,
#            col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5))
#   points(tmp$observed.p + jitter,tmp$med,
#          pch=16, cex = 1,
#          col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5))
#   
#   # one-to-one line
#   abline(a=0,b=1)
#   
#   # add mle estimate for the year
#   abline(h = mle, lty='dashed')
#   
# }
# dev.off()
# 
# # compare pooling for upper level parameters
# #pdf(file=paste0(dirFigures,"appendix-x-figure54-2.pdf"), width=8, height=6)
# 
# par(mfrow=c(1,1))
#   
#   tmp = simData %>% 
#     dplyr::group_by(year) %>%
#     dplyr::summarise(y = sum(fruitingPlantNumber), n = sum(seedlingNumber)) %>%
#     dplyr::mutate(mle.year = y/n) %>%
#     dplyr::arrange(mle.year)
#   
#   mcmcSummary<-samplesBinLinkBetaPriorPartialMeanHier %>%
#     tidybayes::spread_draws(phi[site,year]) %>%
#     dplyr::group_by(year) %>%
#     dplyr::summarise(med.h = median(phi), 
#                      ci.lo.h = quantile(phi,probs=0.025), 
#                      ci.hi.h = quantile(phi,probs=0.975))
#   
#   mcmcSummary2<-samplesBinLinkBetaPriorPartialMean %>%
#     tidybayes::spread_draws(phi[year]) %>%
#     dplyr::group_by(year) %>%
#     dplyr::summarise(med = median(phi), 
#                      ci.lo = quantile(phi,probs=0.025), 
#                      ci.hi = quantile(phi,probs=0.975))
#   
#   tmp<-tmp %>%
#     dplyr::left_join(mcmcSummary,by="year") %>%
#     dplyr::left_join(mcmcSummary2,by="year")
#   
#   # filter data to year
#   mle.pop = sum(simData$fruitingPlantNumber)/sum(simData$seedlingNumber)
#   
#   # plot hierarchical population-level estimate
#   
# 
#   # toggle on jitter to avoid overplotting
#   jitter <- rnorm(length(tmp$mle.year),mean=0,sd=0.01)
#   
#   plot(tmp$mle.year + jitter,tmp$med,
#        pch=16, cex = 1,
#        col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5),
#        xlim=c(0,1),ylim=c(0,1),
#        xlab="Year-level MLE, sum(y[i])/sum(n[i])",
#        ylab="phi[i] and 95% credible interval",type='n')
# 
#   # toggle on jitter to avoid overplotting
#   jitter <- rnorm(length(tmp$mle.year),mean=0,sd=0.01)
#   
#   # 2-level hierarchical model
#   segments(x0=tmp$mle.year + jitter,x1=tmp$mle.year + jitter,
#            y0=tmp$ci.lo.h,y1=tmp$ci.hi.h,lwd=0.75,
#            col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
#   points(tmp$mle.year + jitter,tmp$med.h,
#          pch=16, cex = 1,
#          col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
#   
#   # 1-level hierarchical model
#   segments(x0=tmp$mle.year - jitter,x1=tmp$mle.year - jitter,
#            y0=tmp$ci.lo,y1=tmp$ci.hi,lwd=0.75,
#            col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5))
#   points(tmp$mle.year - jitter,tmp$med,
#          pch=16, cex = 1,
#          col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5))
#   
#   # one-to-one line
#   abline(a=0,b=1)
#   
#   # add mle estimate for the population
#   abline(h = mle.pop, lty='dashed')
#   
#   
#   mcmcSummaryPop<-samplesBinLinkBetaPriorPartialMeanHier %>%
#     tidybayes::spread_draws(phi0) %>%
#     dplyr::summarise(med = median(phi0), 
#                      ci.lo = quantile(phi0,probs=0.025), 
#                      ci.hi = quantile(phi0,probs=0.975))
#   
#   # add mle estimate for the population
#   abline(h = mcmcSummaryPop$med, lty=,
#          col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
#  
#   out=seq(-.2,1.2,.1)
#   yhi = rep(mcmcSummaryPop$ci.hi,length(out))
#   ylo = rep(mcmcSummaryPop$ci.lo,length(out))
# 
#   abline(h = mcmcSummaryPop$med, 
#          col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
# 
#   polygon(x=c(out,rev(out)),y=c(ylo,rev(yhi)),
#           col=rgb(red = 0, green = 0, blue = 0, alpha = 0.1),
#           border=NA)
# 
# 
#   # -------------------------------------------------------------------
#   # -------------------------------------------------------------------
# # Compare parameterizations
#   # -------------------------------------------------------------------
#   # -------------------------------------------------------------------
#   samplesBinLinkBetaPriorPartialMean %<>% tidybayes::recover_types(simData)
#   samplesBinLinkBetaPriorPartialMeanLogit %<>% tidybayes::recover_types(simData)
#   samplesBinLinkBetaPriorPartialMeanLogitNC %<>% tidybayes::recover_types(simData)
#   
#   ref<-simData %>% 
#     dplyr::select(year,plot) %>%
#     tidyr::unite(yearPlot,c(year,plot))
#   
#   
#   tmp = simData %>% 
#    # dplyr::filter(year==i) %>%
#      dplyr::mutate(observed.p = fruitingPlantNumber/seedlingNumber) %>%
#     dplyr::arrange(observed.p) %>%
#     tidyr::unite(yearPlot,c(year,plot))
#   
#   mcmcSummary<-samplesBinLinkBetaPriorPartialMean %>%
#     tidybayes::spread_draws(theta[year,plot]) %>%
#     tidyr::unite(yearPlot,c(year,plot)) %>%
#    # dplyr::filter(yearPlot %in% tmp$yearPlot) %>%
#     dplyr::group_by(yearPlot) %>%
#     dplyr::summarise(med = median(theta), 
#                      ci.lo = quantile(theta,probs=0.025), 
#                      ci.hi = quantile(theta,probs=0.975))
#   
#   mcmcSummary2<-samplesBinLinkBetaPriorPartialMeanLogit %>%
#     tidybayes::spread_draws(alpha[year,plot]) %>%
#     tidyr::unite(yearPlot,c(year,plot)) %>%
#     dplyr::mutate(theta=boot::inv.logit(alpha)) %>%
#     #dplyr::filter(yearPlot %in% tmp$yearPlot) %>% 
#     dplyr::group_by(yearPlot) %>%
#     dplyr::summarise(med.logit = median(theta), 
#                      ci.lo.logit = quantile(theta,probs=0.025), 
#                      ci.hi.logit = quantile(theta,probs=0.975))
#   
#   mcmcSummary3<-samplesBinLinkBetaPriorPartialMeanLogitNC %>%
#     tidybayes::spread_draws(alpha[year,plot]) %>%
#     tidyr::unite(yearPlot,c(year,plot)) %>%
#     dplyr::mutate(theta=boot::inv.logit(alpha)) %>%
#     #dplyr::filter(yearPlot %in% tmp$yearPlot) %>% 
#     dplyr::group_by(yearPlot) %>%
#     dplyr::summarise(med.logit.nc = median(theta), 
#                      ci.lo.logit.nc = quantile(theta,probs=0.025), 
#                      ci.hi.logit.nc = quantile(theta,probs=0.975))
#   
#   tmp<-tmp %>%
#     dplyr::left_join(mcmcSummary,by="yearPlot") %>%
#     dplyr::left_join(mcmcSummary2,by="yearPlot") %>%
#     dplyr::left_join(mcmcSummary3,by="yearPlot")
#   
#   mle = sum(tmp$fruitingPlantNumber)/sum(tmp$seedlingNumber)
#   
#   # plot hierarchical population-level estimate
#   
#   # toggle on jitter to avoid overplotting
#   jitter <- rnorm(length(tmp$observed.p),mean=0,sd=0.01)
#   
#   par(mfrow=c(1,2))
#   plot(tmp$observed.p + jitter,tmp$med,
#        pch=16, cex = 1,
#        col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5),
#        xlim=c(0,1),ylim=c(0,1),
#        xlab="Observed proportion, y[i]/n[i]",
#        ylab="theta[i] and 95% credible interval",
#        main="",type='n')
#   
#   segments(x0=tmp$observed.p + jitter,x1=tmp$observed.p + jitter,
#            y0=tmp$ci.lo,y1=tmp$ci.hi,lwd=0.75,
#            col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5))
#   points(tmp$observed.p + jitter,tmp$med,
#          pch=16, cex = 1,
#          col=rgb(red = .5, green = 0, blue = .5, alpha = 1))
#   
#   # one-to-one line
#   abline(a=0,b=1)
#   
#   # add mle estimate for the year
#   abline(h = mle, lty='dashed')
#   
#   plot(tmp$observed.p + jitter,tmp$med.logit,
#        pch=16, cex = 1,
#        col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5),
#        xlim=c(0,1),ylim=c(0,1),
#        xlab="Observed proportion, y[i]/n[i]",
#        ylab="theta[i] and 95% credible interval",
#        main="",type='n')
#   
#   segments(x0=tmp$observed.p + jitter,x1=tmp$observed.p + jitter,
#            y0=tmp$ci.lo.logit,y1=tmp$ci.hi.logit,lwd=0.75,
#            col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5))
#   points(tmp$observed.p + jitter,tmp$med.logit,
#          pch=16, cex = 1,
#          col=rgb(red = .5, green = 0, blue = .5, alpha = 1))
#   
#   # one-to-one line
#   abline(a=0,b=1)
#   
#   # add mle estimate for the year
#   abline(h = mle, lty='dashed')
# 
#   par(mfrow=c(1,3))
#   # compare  
#   plot(tmp$med ,tmp$med.logit,
#        pch=16,cex=.75,
#        xlim=c(0,1),ylim=c(0,1),
#        xlab="theta[i] and 95% credible interval\n direct parameterization",
#        ylab="theta[i] and 95% credible interval, logit-parameterization")
#   # one-to-one line
#   abline(a=0,b=1)
#   
#   #compare
#   plot(tmp$med ,tmp$med.logit.nc,
#        pch=16,cex=.75,
#        xlim=c(0,1),ylim=c(0,1),
#        xlab="theta[i] and 95% credible interval\n direct parameterization",
#        ylab="theta[i] and 95% credible interval, non-centered logit-parameterization")
#   # one-to-one line
#   abline(a=0,b=1)
#   
#   #compare
#   plot(tmp$med.logit ,tmp$med.logit.nc,
#        pch=16,cex=.75,
#        xlim=c(0,1),ylim=c(0,1),
#        xlab="theta[i] and 95% credible interval\n logit-parameterization",
#        ylab="theta[i] and 95% credible interval, non-centered logit-parameterization")
#   # one-to-one line
#   abline(a=0,b=1)
#   
#   
#   
# dim(MCMCchains(samplesBinLinkBetaPriorPartialMean,params="theta"))
# dim(MCMCchains(samplesBinLinkBetaPriorPartialMeanLogit,params='alpha'))
# theta <- MCMCchains(samplesBinLinkBetaPriorPartialMean,params="theta") %>%
#   recover_types(simData) %>%
#   spread_draws(theta[year,plot]) %>%
#   tidyr::unite(yearPlot, c(year,plot),remove=FALSE) %>%
#   dplyr::filter(yearPlot %in% ref$yearPlot) %>% 
#   unspread_draws(theta[year,plot]) %>% 
#   dplyr::select(-c('.chain','.iteration','.draw'))
#   
# theta.logit<- MCMCchains(samplesBinLinkBetaPriorPartialMeanLogit,params="alpha") %>%
#   recover_types(simData) %>%
#   spread_draws(alpha[year,plot]) %>%
#   dplyr::mutate(theta = boot::inv.logit(alpha)) %>%
#   dplyr::select(-alpha) %>%
#   tidyr::unite(yearPlot, c(year,plot),remove=FALSE) %>%
#   dplyr::filter(yearPlot %in% ref$yearPlot) %>% 
#   unspread_draws(theta[year,plot]) %>% 
#   dplyr::select(-c('.chain','.iteration','.draw'))
#   
# delta=theta.logit-theta
# delta.summary<-apply(delta,2,quantile,probs=c(.025,.5,.975))  
# 
# test<-delta.summary %>%
#   tidybayes::spread_draws(theta[year,plot]) %>%
#   dplyr::filter(.draw==2) %>%
#   dplyr::left_join(simData,by=c("year","plot"))
# 
# par(mfrow=c(1,1))
# plot(test$seedlingNumber,test$theta)
# 
# 
# par(mfrow=c(1,1))
# 
# tmp = simData %>% 
#   dplyr::group_by(year) %>%
#   dplyr::summarise(y = sum(fruitingPlantNumber), n = sum(seedlingNumber)) %>%
#   dplyr::mutate(mle.year = y/n) %>%
#   dplyr::arrange(mle.year)
# 
# mcmcSummary<-samplesBinLinkBetaPriorPartialMean %>%
#   tidybayes::spread_draws(phi[year]) %>%
#   dplyr::group_by(year) %>%
#   dplyr::summarise(med = median(phi), 
#                    ci.lo = quantile(phi,probs=0.025), 
#                    ci.hi = quantile(phi,probs=0.975))
# 
# mcmcSummary2<-samplesBinLinkBetaPriorPartialMeanLogit %>%
#   tidybayes::spread_draws(mu[year]) %>%
#   dplyr::group_by(year) %>%
#   dplyr::summarise(med.logit = median(boot::inv.logit(mu)), 
#                    ci.lo.logit = quantile(boot::inv.logit(mu),probs=0.025), 
#                    ci.hi.logit = quantile(boot::inv.logit(mu),probs=0.975))
# 
# tmp<-tmp %>%
#   dplyr::left_join(mcmcSummary,by="year") %>%
#   dplyr::left_join(mcmcSummary2,by="year")
# 
# # plot hierarchical population-level estimate
# 
# 
# # toggle on jitter to avoid overplotting
# jitter <- rnorm(length(tmp$mle.year),mean=0,sd=0.01)
# 
# plot(tmp$mle.year + jitter,tmp$med,
#      pch=16, cex = 1,
#      col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5),
#      xlim=c(0,1),ylim=c(0,1),
#      xlab="Year-level MLE, sum(y[i])/sum(n[i])",
#      ylab="phi[i] and 95% credible interval",type='n')
# 
# # toggle on jitter to avoid overplotting
# jitter <- rnorm(length(tmp$mle.year),mean=0,sd=0.01)
# 
# # 1-level hierarchical model, direct param
# segments(x0=tmp$mle.year + jitter,x1=tmp$mle.year + jitter,
#          y0=tmp$ci.lo,y1=tmp$ci.hi,lwd=0.75,
#          col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
# points(tmp$mle.year + jitter,tmp$med,
#        pch=16, cex = 1,
#        col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
# 
# # 1-level hierarchical model, logit param
# segments(x0=tmp$mle.year - jitter,x1=tmp$mle.year - jitter,
#          y0=tmp$ci.lo.logit,y1=tmp$ci.hi.logit,lwd=0.75,
#          col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5))
# points(tmp$mle.year - jitter,tmp$med.logit,
#        pch=16, cex = 1,
#        col=rgb(red = .5, green = 0, blue = .5, alpha = 0.5))
# 
# # one-to-one line
# abline(a=0,b=1)
# 
# # add mle estimate for the population
# abline(h = mle.pop, lty='dashed')
# 
# 
# mcmcSummaryPop<-samplesBinLinkBetaPriorPartialMeanHier %>%
#   tidybayes::spread_draws(phi0) %>%
#   dplyr::summarise(med = median(phi0), 
#                    ci.lo = quantile(phi0,probs=0.025), 
#                    ci.hi = quantile(phi0,probs=0.975))
# 
# # add mle estimate for the population
# abline(h = mcmcSummaryPop$med, lty=,
#        col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
# 
# out=seq(-.2,1.2,.1)
# yhi = rep(mcmcSummaryPop$ci.hi,length(out))
# ylo = rep(mcmcSummaryPop$ci.lo,length(out))
# 
# abline(h = mcmcSummaryPop$med, 
#        col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
# 
# polygon(x=c(out,rev(out)),y=c(ylo,rev(yhi)),
#         col=rgb(red = 0, green = 0, blue = 0, alpha = 0.1),
#         border=NA)
# 
# 
# 
# tmp = simData %>% 
#   # dplyr::filter(year==i) %>%
#   dplyr::mutate(observed.p = fruitingPlantNumber/seedlingNumber) %>%
#   dplyr::arrange(observed.p) %>%
#   tidyr::unite(yearPlot,c(year,plot))
# 
#   samplesBinLinkBetaPriorComplete <- readRDS(simFiles[[5]])
#   
#   samplesBinLinkBetaPriorPartialMean <- readRDS(simFiles[[7]])
#   samplesBinLinkBetaPriorPartialMeanLogit <- readRDS(simFiles[[10]])
#   samplesBinLinkBetaPriorPartialMeanLogitNC <- readRDS(simFiles[[11]])
#   
#   samplesBinLinkBetaPriorPartialMeanHier <- readRDS(simFiles[[6]])
#   samplesBinLinkBetaPriorPartialMeanHierLogit <- readRDS(simFiles[[8]])
#   samplesBinLinkBetaPriorPartialMeanHierLogitNC <- readRDS(simFiles[[9]])
#   
# 
#   
#   # Comparison shows that the distribution of rhat for
#   # the logit parameterization is <1.05
#   # plot comparing Rhat values 
#   par(mfrow=c(1,3))  
#   hist(MCMCsummary(samplesBinLinkBetaPriorPartialMean,params=c("theta"))$Rhat,breaks=40)
#   hist(MCMCsummary(samplesBinLinkBetaPriorPartialMeanLogit,params=c("theta"))$Rhat,breaks=40)
#   hist(MCMCsummary(samplesBinLinkBetaPriorPartialMeanLogitNC,params=c("theta"))$Rhat,breaks=40)
# 
#   par(mfrow=c(1,3))  
#   MCMCsummary(samplesBinLinkBetaPriorPartialMean,params=c("phi"))
#   MCMCsummary(samplesBinLinkBetaPriorPartialMeanLogit,params=c("mu"))
#   MCMCsummary(samplesBinLinkBetaPriorPartialMeanLogitNC,params=c("mu"))
#   
#   par(mfrow=c(1,3))  
#   hist(MCMCsummary(samplesBinLinkBetaPriorPartialMeanHier,params=c("theta"))$Rhat,breaks=40)
#   hist(MCMCsummary(samplesBinLinkBetaPriorPartialMeanHierLogit,params=c("theta"))$Rhat,breaks=40)
#   hist(MCMCsummary(samplesBinLinkBetaPriorPartialMeanHierLogitNC,params=c("theta"))$Rhat,breaks=40)
#   
#   par(mfrow=c(1,3))  
#   MCMCsummary(samplesBinLinkBetaPriorPartialMeanHier,params=c("phi"))
#   MCMCsummary(samplesBinLinkBetaPriorPartialMeanHierLogit,params=c("mu"))
#   MCMCsummary(samplesBinLinkBetaPriorPartialMeanHierLogitNC,params=c("mu"))
#     
#   # look at posteriors
# phi<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params=c("phi"))
# kappa<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params=c("kappa"))
# theta<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params=c("theta"))
# 
# # see neal's funnel
# par(mfrow=c(2,5))
# for(i in 1:nYears){
# plot(boot::logit(phi[,i]),log(kappa[,i]))
# }
# 
# par(mfrow=c(1,1))
# # also the relationship of theta and kappa
# plot(boot::logit(theta[,1]),log(kappa[,1]))
# 
# 
# # look at posteriors
# mu<-MCMCchains(samplesBinLinkBetaPriorPartialMeanLogit,params=c("mu"))
# sigma<-MCMCchains(samplesBinLinkBetaPriorPartialMeanLogit,params=c("sigma"))
# alpha<-MCMCchains(samplesBinLinkBetaPriorPartialMeanLogit,params=c("alpha"))
# theta<-MCMCchains(samplesBinLinkBetaPriorPartialMeanLogit,params=c("theta"))
# 
# par(mfrow=c(1,1))
# # also the relationship of log-odds of success (alpha) and log population scale
# plot(alpha[,1],log(sigma[,1]))
# 
# 
# # look at posteriors non-centered parm
# mu<-MCMCchains(samplesBinLinkBetaPriorPartialMeanLogitNC,params=c("mu"))
# sigma<-MCMCchains(samplesBinLinkBetaPriorPartialMeanLogitNC,params=c("sigma"))
# alpha.std<-MCMCchains(samplesBinLinkBetaPriorPartialMeanLogitNC,params=c("alpha.std"))
# theta<-MCMCchains(samplesBinLinkBetaPriorPartialMeanLogitNC,params=c("theta"))
# 
# par(mfrow=c(1,1))
# plot(mu[,1],log(sigma[,1]))
# 
# # also the relationship of log-odds of success (alpha) and log population scale
# plot(alpha.std[,1],log(sigma[,1]))
# 
# # -------------------------------------------------------------------
# # -------------------------------------------------------------------
# # Posterior Predictive Checks
# 
# # notes: check group = interaction(dat$site,dat$yearStart)
# # for stat density grouped plots
# # -------------------------------------------------------------------
# # -------------------------------------------------------------------
# library(gridExtra)
# library(bayesplot)
# library(dplyr)
# 
# 
# iter<-dim(MCMCchains(mcmcSamples,params="fruitingPlantNumberSim"))[1]
# 
# #stat:skew
# skew <- psych::skew
# 
# # pdf(
# #   "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/products/seedBagsCompletePoolingPPC.pdf",
# #   onefile=TRUE,
# #   paper="USr",
# #   height = 5.5, width = 10)
# 
# # ppc_dens_overlay(simData$fruitingPlantNumber, MCMCchains(samplesBinLinkBetaPriorPartialMeanHier,params="fruitingPlantNumberSim")[sample(iter,1000), ]) +
# #   theme_bw() +
# #   labs(title="Posterior predictive checks for number of fruiting plants", 
# #   caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")
# 
# y <- simData$fruitingPlantNumber
# yrep.nh <- MCMCchains(samplesBinLinkBetaPriorComplete,params="fruitingPlantNumberSim")[sample(iter,1000), ]
# group <- as.factor(simData$year)
# 
# color_scheme_set("brightblue")
# 
# nh.mean<-ppc_stat_grouped(y, yrep.nh,
#                    group=group,facet_args=list(nrow=2)) +
#   theme_bw() +
#   labs(title="Posterior predictive checks for number of fruiting plants; non-hierarchical model", 
#        caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")
# 
# nh.var<-ppc_stat_grouped(y, yrep.nh, group=group,stat="var",facet_args=list(nrow=2)) +
#   theme_bw() +
#   labs(title="Posterior predictive checks for number of fruiting plants; non-hierarchical model", 
#        caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")
# 
# 
# yrep.hyear = MCMCchains(samplesBinLinkBetaPriorPartialMean,params="fruitingPlantNumberSim")[sample(iter,1000), ]
# 
# color_scheme_set("purple")
# 
# hieryear.mean<-ppc_stat_grouped(y, yrep.hyear,group=group,facet_args=list(nrow=2)) +
#   theme_bw() +
#   labs(title="Posterior predictive checks for number of fruiting plants; hierarchical model with year-level parameters", 
#        caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")
# 
# hieryear.var<-ppc_stat_grouped(y, yrep.hyear,group=group,stat='var',facet_args=list(nrow=2)) +
#   theme_bw()
#   labs(title="Posterior predictive checks for number of fruiting plants; hierarchical model with year-level parameters", 
#        caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")
# 
# 
#  yrep.hyearpop= MCMCchains(samplesBinLinkBetaPriorPartialMeanHier,params="fruitingPlantNumberSim")[sample(iter,1000), ]
# 
#  color_scheme_set("darkgray")
# 
#  
#  hieryearpop.mean<-ppc_stat_grouped(y, yrep.hyearpop,group=group,facet_args=list(nrow=2)) +
#   theme_bw() +
#   labs(title="Posterior predictive checks for number of fruiting plants; hierarchical model with year- and population-level parameters", 
#        caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")
# 
# hieryearpop.var<-ppc_stat_grouped(y, yrep.hyearpop,group=group,stat='var',facet_args=list(nrow=2)) +
#   theme_bw() +
#   labs(title="Posterior predictive checks for number of fruiting plants; hierarchical model with year- and population-level parameters", 
#        caption="Dark line is the density of observed data (y) and the lighter lines show the densities of Y_rep from 1000 draws of the posterior")
# 
# l = list(nh.mean,hieryear.mean,hieryearpop.mean,nh.var,hieryear.var,hieryearpop.var)
# 
# ggsave(paste0(dirFigures,"appendix-x-ppc.pdf"), marrangeGrob(grobs = l,nrow=1,ncol=1))
# 
# 
# 
# 
# par(mfrow=c(2,5))
# # compare hierarchical to non-hierarchical
# # for posterior of p
# for(i in 1:nYears){
#   # plot hierarchical population-level estimate
#   plot(density(hyearpop.p_site),
#        xlim=c(0,1),ylim=c(0,30),
#        xlab="Probability of survivorship",
#        main=yearNames[i])
#   
#   # plot hierarchical year-level estimates
#   lines(density(hyearpop.p[,i]),lty=1,lwd=1,col='red')
#   
#   lines(density(hyear.p[,i]),lty=1,lwd=1,col='blue')
#   
#   tmp = simData %>% dplyr::filter(year==i)
#   points(tmp$fruitingPlantNumber/tmp$seedlingNumber,
#          rep(0,dim(tmp)[1]),
#          pch=16,
#          col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
# }
# 
# hyear.p.sum<-apply(hyear.p,2,quantile,probs=c(.025,.5,.975))
# hyearpop.p.sum<-apply(hyearpop.p,2,quantile,probs=c(.025,.5,.975))
# 
# # show summary of posteriors with 95% credible intervals
# par(mfrow=c(1,1))
# plot(1:9-0.25,t(hyear.p.sum)[,2],pch=1,ylim=c(0,1),xlim=c(0,10))
# segments(x0=1:9-0.25,y0=t(hyear.p.sum)[,1],
#          x1=1:9-0.25,y1=t(hyear.p.sum)[,3])
# points(1:9+0.25,t(hyearpop.p.sum)[,2],pch=3)
# segments(x0=1:9+0.25,y0=t(hyearpop.p.sum)[,1],
#          x1=1:9+0.25,y1=t(hyearpop.p.sum)[,3])
# 
# yearNames = levels(simData$year)
# 
# # par(mfrow=c(2,5))
# # for(i in 1:nYears){
# #   plot(density(p[,i]),xlim=c(-.1,1.1),main=yearNames[i])
# #   #lines(density(omega[,i]),lty='dotted')
# #   lines(density(phi[,i]),lty='dashed')
# # }
# # 
# # simData$year <- as.numeric(simData$year)
# # ggplot() +
# #   stat_density(data=samplesBinLinkBetaPriorComplete %>%
# #                  tidybayes::spread_draws(p[year]),aes(p),
# #                geom='line') +
# #   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
# #                  tidybayes::spread_draws(phi[year]),aes(phi),color="red",
# #                geom='line') +
# #   facet_wrap(~year,scales="free",nrow=2) +
# #   geom_jitter(data=simData ,aes(x=fruitingPlantNumber/seedlingNumber,y=0)) +
# #   theme_minimal()
# 
# # posteriors for the complete pooling model, the mean of the beta, the mode of the beta
# # it shows that the complete pooling model (one p per year, black)
# # has the smallest credible intervals which generally
# # overlap those from the model parameterized via the mean (red)
# # the models parameterized via the mean are slightly broader credible intervals
# 
# # ggplot() +
# #   stat_density(data=samplesBinLinkBetaPriorComplete %>%
# #                  tidybayes::spread_draws(p[year]),aes(p),geom='line') +
# #   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
# #                  tidybayes::spread_draws(phi[year]),aes(phi),color="red",geom='line') +
# #   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
# #                  tidybayes::spread_draws(p[year]),aes(p),color="red",geom='line',linetype='dotted') +
# #   stat_density(data=samplesBinLinkBetaPriorPartialMode %>%
# #                  tidybayes::spread_draws(omega[year]),aes(omega),color="blue",geom='line') +
# #   stat_density(data=samplesBinLinkBetaPriorPartialMode %>%
# #                  tidybayes::spread_draws(p[year]),aes(p),color="blue",geom='line',linetype='dotted') +
# #   facet_wrap(~year,scales="free",nrow=2) +
# #   geom_jitter(data=simData ,aes(x=fruitingPlantNumber/seedlingNumber,y=0)) +
# #   theme_minimal()
# 
# # the posterior distribution for probability of success in year p (marginalized over all plots)
# # is broader than the posterior for the mean probability of success in year p
# # parameterization via the mean and mode give similar results once marginalized
# # though the mode parameterization seems to give broader posteriors
# 
# simData$year <- as.numeric(as.factor(simData$year))
# simData$plot <- as.numeric(as.factor(simData$plot))
# # ggplot() +
# #   stat_density(data=samplesBinLinkBetaPriorComplete %>%
# #                  tidybayes::spread_draws(p[year]) %>%
# #                  dplyr::filter(year==5),aes(p),lwd=1, geom='line') +
# #   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
# #                  tidybayes::spread_draws(phi[year]) %>% 
# #                  dplyr::filter(year==5),aes(phi),color="red",lwd=1, geom='line') +
# #   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
# #                  tidybayes::spread_draws(p[year]) %>% 
# #                  dplyr::filter(year==5),aes(p),
# #                color="red",linetype='dashed',lwd=0.5, geom="line") +
# #   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
# #                  tidybayes::spread_draws(theta[year,plot]) %>% 
# #                  dplyr::filter(year==5) ,aes(theta,group=plot),
# #                color="black",lwd=0.25, geom="line") +
# #   geom_point(data=simData  %>% 
# #                dplyr::filter(year==5),
# #              aes(x=fruitingPlantNumber/seedlingNumber,y=0),position=ggstance::position_dodgev(height=0.3)) +
# #   theme_minimal() +
# #   facet_wrap(~plot, scales = 'free')
# 
# 
# samplesBinLinkBetaPriorPartialMeanHier %<>% recover_types(simData)
# # plot 2
# ggplot() +
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(phi0),aes(phi0),color="red",geom='line',size=1) +
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(phi[site,year]),aes(phi,group=year),geom='line',size=0.25, position='identity') +
#   geom_point(data=simData,
#              aes(x=fruitingPlantNumber/seedlingNumber,y=0),
#              position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
#   # facet_wrap(~year,scales="free_y",nrow=2) +
#   theme_minimal()
# 
# ggplot() +
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(p_site),aes(p_site),color="red",geom='line',size=1) +
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(p[site,year]),aes(p,group=year),geom='line',size=0.25, position='identity') +
#   geom_point(data=simData,
#              aes(x=fruitingPlantNumber/seedlingNumber,y=0),
#              position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
#   facet_wrap(~year,scales="free_y",nrow=2) +
#   theme_minimal()
# 
# # compare levels of the hierarchy
# ggplot() +
#   # posterior distribution of the mean probability of success at the site
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(phi0),aes(phi0),color="red",geom='line',size=.5) +
#   # posterior distribution of the mean probability of success in each year
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(phi[site,year]),aes(phi,group=year),geom='line',size=0.25, position='identity') +
#   # posterior distribution for the probability of success at the site
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(p_site),aes(p_site),color="red",geom='line',size=1) +
#   # posterior distribution for the probability of success in each year 
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(p[site,year]),aes(p,group=year),geom='line',size=1, position='identity') +
#   geom_point(data=simData,
#              aes(x=fruitingPlantNumber/seedlingNumber,y=0),
#              position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
#   facet_wrap(~year,scales="free_y",nrow=2) +
#   theme_minimal()
# 
# # key takeaway from the plot above: 
# # the posterior for probability of success at each level (site, year-within-site) is broader than the mean probability of success at that level (compare thin to solid lines)
# 
# ggplot() +
#   # posterior distribution of the mean probability of success at the site
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(phi0),aes(phi0),color="red",geom='line',size=1) +
#   # posterior distribution of the mean probability of success in each year
#   # stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#   #                tidybayes::spread_draws(phi[site,year]),aes(phi,group=as.factor(year)),geom='line',size=0.25, position='identity') +
#   # posterior distribution for the probability of success at the site
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(p_site),aes(p_site),color="red",geom='line',size=.5) +
#   # posterior distribution for the probability of success in each year 
#   # stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#   #                tidybayes::spread_draws(p[site,year]),aes(p,group=year),geom='line',size=1, position='identity') +
#   geom_point(data=simData,
#              aes(x=fruitingPlantNumber/seedlingNumber,y=0),
#              position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
#   #facet_wrap(~year,scales="free_y",nrow=2) +
#   theme_minimal()
# 
# ggplot() +
#   # posterior distribution of the mean probability of success at the site
#   # stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#   #                tidybayes::spread_draws(phi0),aes(phi0),color="red",geom='line',size=1) +
#   # posterior distribution of the mean probability of success in each year
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(phi[site,year]),aes(phi,group=as.factor(year)),geom='line',size=0.25, position='identity') +
#   # posterior distribution for the probability of success at the site
#   # stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#   #                tidybayes::spread_draws(p_site),aes(p_site),color="red",geom='line',size=.5) +
#   # posterior distribution for the probability of success in each year 
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(p[site,year]),aes(p,group=year),geom='line',size=1, position='identity') +
#   geom_point(data=simData,
#              aes(x=fruitingPlantNumber/seedlingNumber,y=0),
#              position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
#   facet_wrap(~year,scales="free_y",nrow=2) +
#   theme_minimal()
# 
# # effect of hierarchy on estimate of means - estimates pool towards hierarchy with little data
# ggplot() +
#   # posterior distribution of the mean probability of success at the site
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(phi0),aes(phi0),color="red",geom='line',size=1) +
#   # posterior distribution of the mean probability of success in each year
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(phi[site,year]),aes(phi,group=as.factor(year)),geom='line',size=0.25, position='identity',color='red') +
#   # posterior distribution of the mean probability of success in each year
#   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
#                  tidybayes::spread_draws(phi[year]),aes(phi,group=as.factor(year)),geom='line',size=0.25, position='identity') +
#   geom_point(data=simData,
#              aes(x=fruitingPlantNumber/seedlingNumber,y=0),
#              position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
#   facet_wrap(~year,scales="free_y",nrow=2) +
#   theme_minimal()
# 
# 
# 
# 
# 
# # compare the three different fits (no hierarchy, hierarchy at 1 level, hierarchy at 2 levels)
# # difference is greatest with small datasets (year = 9 only has 1 data point)
# ggplot() +
#   stat_density(data=samplesBinLinkBetaPriorComplete %>%
#                  tidybayes::spread_draws(p[year]),aes(p),color="black",geom='line',size=1) +
#   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
#                  tidybayes::spread_draws(phi[year]),aes(phi),color="red",geom='line') +
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(phi[site,year]),aes(phi),color="blue",geom='line') +
#   geom_point(data=simData,
#              aes(x=fruitingPlantNumber/seedlingNumber,y=0),
#              position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
#   facet_wrap(~year,scales="free_y",nrow=2) +
#   theme_minimal()
# 
# # compare distribution of p
# ggplot() +
#   stat_density(data=samplesBinLinkBetaPriorComplete %>%
#                  tidybayes::spread_draws(p[year]),aes(p),color="black",geom='line',size=1) +
#   stat_density(data=samplesBinLinkBetaPriorPartialMean %>%
#                  tidybayes::spread_draws(p[year]),aes(p),color="red",geom='line') +
#   stat_density(data=samplesBinLinkBetaPriorPartialMeanHier %>%
#                  tidybayes::spread_draws(p[site,year]),aes(p),color="blue",geom='line') +
#   geom_point(data=simData,
#              aes(x=fruitingPlantNumber/seedlingNumber,y=0),
#              position=ggstance::position_dodgev(height=0.3),alpha=0.5) +
#   facet_wrap(~year,scales="free_y",nrow=2) +
#   theme_minimal()
# 
# 
# nohier <-apply(MCMCchains(samplesBinLinkBetaPriorComplete,params='p'),2,median)
# hier1 <- apply(MCMCchains(samplesBinLinkBetaPriorPartialMean,params='phi'),2,median)
# hier2 <-apply(MCMCchains(samplesBinLinkBetaPriorPartialMeanHier,params='phi'),2,median)
# 
# plot(nohier,hier1);abline(a=0,b=1)
# plot(nohier,hier2);abline(a=0,b=1)
# plot(hier1,hier2);abline(a=0,b=1)
# 
# samplesBinLinkBetaPriorPartialMean %>%
#   tidybayes::spread_draws(theta[year,plot])
# 
# ref<-simData %>% 
#   dplyr::select(year,plot) %>%
#   tidyr::unite(yearPlot,c(year,plot))
# 
# theta.med<-samplesBinLinkBetaPriorPartialMean %>%
#   tidybayes::spread_draws(theta[year,plot]) %>%
#   tidyr::unite(yearPlot,c(year,plot)) %>%
#   dplyr::filter(yearPlot %in% ref$yearPlot) %>%
#   tidyr::separate(yearPlot, c("year","plot")) %>%
#   dplyr::group_by(year,plot) %>%
#   dplyr::summarise(theta.med = median(theta))
# 
# theta.med$year <- as.numeric(theta.med$year)
# theta.med$plot <- as.numeric(theta.med$plot)
# 
# out<-simData %>%
#   dplyr::left_join(theta.med,by=c("year","plot"))
# 
# par(mfrow=c(1,1))
# plot(out$fruitingPlantNumber/out$seedlingNumber,
#      out$theta.med)
# 
# plot(out$seedlingNumber,out$theta.med-out$fruitingPlantNumber/out$seedlingNumber)
# 
# p.summary<-samplesBinLinkBetaPriorPartialMean %>%
#   tidybayes::spread_draws(p[year]) %>%
#   dplyr::summarise(p.med = median(p))
# p.summary$year <- as.numeric(p.summary$year)
# 
# phi.summary<-samplesBinLinkBetaPriorPartialMean %>%
#   tidybayes::spread_draws(phi[year]) %>%
#   dplyr::summarise(phi.med = median(phi))
# phi.summary$year <- as.numeric(phi.summary$year)
# 
# ggplot() +
#   geom_point(data=out,aes(x=fruitingPlantNumber/seedlingNumber,y=theta.med)) +
#   geom_abline(intercept=0,slope=1,lwd=0.5) +
#   geom_hline(data=p.summary,aes(yintercept=p.med),size=0.5) +
#   geom_hline(data=phi.summary,aes(yintercept=phi.med),size=0.5,color='red') +
#   
#   facet_wrap(~year) +
#   theme_minimal()
# 
# 
# theta.med<-samplesBinLinkBetaPriorPartialMean %>%
#   tidybayes::spread_draws(p[year]) %>%
#   dplyr::group_by(year) %>%
#   dplyr::summarise(p.med = median(p))
# 
# theta.med$year <- as.numeric(theta.med$year)
# theta.med$plot <- as.numeric(theta.med$plot)
# 
# out<-simData %>%
#   dplyr::left_join(theta.med,by=c("year","plot"))
# 
# par(mfrow=c(1,1))
# plot(out$fruitingPlantNumber/out$seedlingNumber,
#      out$theta.med)
# 
# 
# thetaBinLinkBetaPriorPartialMode <- MCMCchains(samplesBinLinkBetaPriorPartialMode,params="theta")
# 
# # plot(density(thetaBinLinkBetaPriorPartialMode[,1]), type='n', 
# #      xlim=c(-.2,1.2),ylim=c(0,30))
# # for(i in 211:241){
# #   lines(density(thetaBinLinkBetaPriorPartialMode[,i]))
# # }
# 
# omega<-MCMCchains(samplesBinLinkBetaPriorPartialMode,params="omega")
# 
# plot(density(omega[,1]), type='n',
#      xlim=c(-.2,1.2),ylim=c(0,80))
# for(i in 1:10){
#   lines(density(omega[,i]))
# }
# 
# kappa<-MCMCchains(samplesBinLinkBetaPriorPartialMode,params="kappa")
# 
# alpha = omega*(kappa-2)+1
# beta = (1-omega)*(kappa-2)+1
# plot(seq(0,1,by=0.001),
#      dbeta(seq(0,1,by=0.001),median(alpha),median(beta)),
#      type='n',
#      ylim=c(0,20))
# for(i in 1:20){
#   lines(seq(0,1,by=0.001),
#         dbeta(seq(0,1,by=0.001),median(alpha[i]),median(beta[i])),
#         type='l',lwd=.5,col='red')
# }
# 
# 
# thetaBinLinkBetaPriorPartialMean <- MCMCchains(samplesBinLinkBetaPriorPartialMean,params="theta")
# 
# plot(density(thetaBinLinkBetaPriorPartialMean[,1]), type='n', 
#      xlim=c(-.2,1.2),ylim=c(0,30))
# for(i in 211:241){
#   lines(density(thetaBinLinkBetaPriorPartialMean[,i]))
# }
# 
# phi<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params="phi")
# 
# 
# plot(density(phi[,1]), type='n',
#      xlim=c(-.2,1.2),ylim=c(0,80))
# for(i in 1:10){
#   lines(density(phi[,i]))
# }
# 
# 
# kappa<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params="kappa")
# alpha = kappa*phi
# beta = kappa*(1-phi)
# plot(seq(0,1,by=0.001),dbeta(seq(0,1,by=0.001),median(alpha),median(beta)),
#      type='n')
# for(i in 1:20){
#   lines(seq(0,1,by=0.001),
#         dbeta(seq(0,1,by=0.001),median(alpha[i]),median(beta[i])),
#         type='l',lwd=.5,col='red')
# }
# 
# p<-MCMCchains(samplesBinLinkBetaPriorComplete,params='p')
# 
# # compare posteriors for the hyperparameters
# # the parameterization via the mode seems to sample better in general
# phi<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params='phi')
# kappa.mean<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params='kappa')
# thin <- function(x) sample(x,1000)
# par(mfrow=c(2,5))
# for(i in 1:10){
#   plot(thin(phi[,i]),thin(kappa.mean[,i]),cex=.5,pch=16)
# }
# 
# 
# omega<-MCMCchains(samplesBinLinkBetaPriorPartialMode,params='omega')
# kappa.mode<-MCMCchains(samplesBinLinkBetaPriorPartialMean,params='kappa')
# thin <- function(x) sample(x,1000)
# par(mfrow=c(2,5))
# for(i in 1:10){
#   plot(thin(omega[,i]),thin(kappa.mode[,i]),cex=.5,pch=16)
# }
# 
# plot(apply(p,2,median),apply(phi,2,median),xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)
# plot(apply(p,2,median),apply(omega,2,median),xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)
# plot(apply(phi,2,median),apply(omega,2,median),xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)
# 
# stat.summary <- function(x){
#   
#   x %<>% recover_types(simData)
#   
#   tmp <- x %>% 
#     tidybayes::spread_draws(p[year]) %>%
#     dplyr::summarise(ci.lo = quantile(theta, prob = c(.025)),
#                      med = quantile(theta, prob = c(.5)),
#                      ci.hi = quantile(theta, prob = c(.975)))
#   return(tmp)
# }
# 
# thetaBinLinkBetaPriorPartialMode <- stat.summary(samplesBinLinkBetaPriorPartialMode)
# thetaBinLinkBetaPriorPartialMean <- stat.summary(samplesBinLinkBetaPriorPartialMean)
# 
# d<-data.frame(thetaBinLinkBetaPriorPartialMode,thetaBinLinkBetaPriorPartialMean) %>% 
#   dplyr::mutate(abs(med.1-med)) %>%
#   dplyr::left_join(simData,by=c('year','plot'))
# 
# d<-d %>% dplyr::filter(!is.na(fruitingPlantNumber))
# 
# plot(d$med,d$med.1, 
#      pch=16,cex = 0.5,xlim=c(0,1),ylim=c(0,1))
# # segments(x0=thetaBinLinkBetaPriorPartialMode$med,y0=thetaBinLinkBetaPriorPartialMean$ci.lo,
# #          x1=thetaBinLinkBetaPriorPartialMode$med,y1=thetaBinLinkBetaPriorPartialMean$ci.hi)
# # segments(x0=thetaBinLinkBetaPriorPartialMode$ci.lo,y0=thetaBinLinkBetaPriorPartialMean$med,
# #          x1=thetaBinLinkBetaPriorPartialMode$ci.hi,y1=thetaBinLinkBetaPriorPartialMean$med)
# abline(a=0,b=1)
# 
# mle<-simData %>%
#   dplyr::group_by(year,plot) %>%
#   dplyr::summarise(y = sum(fruitingPlantNumber), n = sum(seedlingNumber)) %>%
#   dplyr::mutate(pHat = y/n) %>%
#   dplyr::left_join(thetaBinLinkBetaPriorPartialMode)
# 
# points(x=mle$med,y=mle$pHat,pch=2)
# 
# ## summary
# MCMCsummary(samplesBinLinkBetaPriorPartialMean,params=c("phi","kappa"))
# MCMCsummary(samplesBinLinkBetaPriorPartialMode,params=c("omega","kappa"))
# 
# samplesBinLinkBetaPriorPartialMean %<>% recover_types(simData)
# 
# d1 <- samplesBinLinkBetaPriorPartialMean %>%
#   tidybayes::spread_draws(theta[year,plot]) %>%
#   dplyr::filter(year==2014)
# 
# d1.summary <- d1 %>% dplyr::group_by(plot) %>%
#   dplyr::summarise(mean.med=median(theta))
# 
# samplesBinLinkBetaPriorPartialMode %<>% recover_types(simData)
# 
# d2 <- samplesBinLinkBetaPriorPartialMode %>%
#   tidybayes::spread_draws(theta[year,plot]) %>%
#   dplyr::filter(year==2014)
# 
# d2.summary <- d2 %>% dplyr::group_by(plot) %>%
#   dplyr::summarise(mode.med=median(theta))
# 
# d <- d1.summary %>%
#   dplyr::left_join(d2.summary, by= "plot")
# plot(d$mode.med,d$mean.med,xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)
# 
# simData %>%
#   dplyr::filter(year==2014) %>%
#   dplyr::left_join(d,by="plot") %>% View
# pdf(file=paste0(dirFigures,"appendix-x-mle_bayes.pdf"), width=8, height=4)
# 


mcmcSummary1<-mcmcSamples1 %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p_1[siteNames]) %>%
  dplyr::summarise(med = median(p_1), 
                   ci.lo = quantile(p_1,probs=0.025), 
                   ci.hi = quantile(p_1,probs=0.975),
                   ci.lo2 = quantile(p_1,probs=0.25), 
                   ci.hi2 = quantile(p_1,probs=0.75)
  )
mcmcSummary1<-cbind(mcmcSummary1,site=siteNames)

vr.site = mcmcSummary1 %>% dplyr::select(site,med)

siteIndex <- data.frame(siteIndex=unique(data$siteBags),siteBags=1:2)
yearIndex <- data.frame(yearBags=1:3,yearIndex=2006:2008)

mcmcSummary1<-mcmcSamples1 %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p0_1[siteBags,yearBags]) %>%
  dplyr::group_by(siteBags,yearBags) %>%
  dplyr::summarise(med = median(p0_1), 
                   ci.lo = quantile(p0_1,probs=0.025), 
                   ci.hi = quantile(p0_1,probs=0.975),
                   ci.lo2 = quantile(p0_1,probs=0.25), 
                   ci.hi2 = quantile(p0_1,probs=0.75)
  )

mcmcSummary1<-mcmcSummary1 %>%
  dplyr::mutate(yearBags = as.integer(yearBags)) %>%
  #dplyr::left_join(siteIndex,by="siteBags") %>%
  dplyr::left_join(yearIndex,by="yearBags") %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(yearBags)) %>%
  #dplyr::select(-c(siteBags)) %>%
  #dplyr::rename(site = siteIndex) %>%
  dplyr::rename(site = siteBags) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)


siteIndex <- data.frame(siteIndex=unique(data$siteBags),site=1:2)
yearIndex <- data.frame(year=1:3,yearIndex=2006:2008)

interannualG1 <- mcmcSamples1 %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p0_1[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(lambda.med = median(p0_1), 
                   ci.lo = quantile(p0_1,probs=0.27), 
                   ci.hi = quantile(p0_1,probs=0.83),
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  dplyr::select(-c(site,year)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::left_join(vr.site, by="site") %>%
  dplyr::arrange(med) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(desc(lambda.med),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  mutate(id = row_number()) %>% 
  ggplot(aes(x = id , y = lambda.med)) + 
  geom_hline(aes(yintercept=med),linetype='dotted') +
  
  geom_point(aes(color=year)) +
  
  geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +
  
  coord_flip() +
  facet_grid(site ~ ., scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ylab("Probability of germination [P(G)]") +
  ylim(c(0,1))

################################################################################
# Germination 1
#################################################################################

mcmcSummary2<-mcmcSamples2 %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p_1[siteBags]) %>%
  dplyr::group_by(siteBags) %>%
  dplyr::summarise(med = median(p_1), 
                   ci.lo = quantile(p_1,probs=0.025), 
                   ci.hi = quantile(p_1,probs=0.975),
                   ci.lo2 = quantile(p_1,probs=0.25), 
                   ci.hi2 = quantile(p_1,probs=0.75)
  )
mcmcSummary2<-cbind(mcmcSummary2,site=siteNames)

vr.site = mcmcSummary2 %>% dplyr::select(site,med)

siteIndex <- data.frame(siteIndex=unique(data$siteBags),siteBags=1)
yearIndex <- data.frame(yearBags=1:3,yearIndex=2006:2008)

mcmcSummary2<-mcmcSamples2 %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p0_1[siteBags,yearBags]) %>%
  dplyr::group_by(siteBags,yearBags) %>%
  dplyr::summarise(med = median(p0_1), 
                   ci.lo = quantile(p0_1,probs=0.025), 
                   ci.hi = quantile(p0_1,probs=0.975),
                   ci.lo2 = quantile(p0_1,probs=0.25), 
                   ci.hi2 = quantile(p0_1,probs=0.75)
  )

mcmcSummary2<-mcmcSummary2 %>%
  dplyr::mutate(yearBags = as.integer(yearBags)) %>%
  #dplyr::left_join(siteIndex,by="siteBags") %>%
  dplyr::left_join(yearIndex,by="yearBags") %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(yearBags)) %>%
  #dplyr::select(-c(siteBags)) %>%
  #dplyr::rename(site = siteIndex) %>%
  dplyr::rename(site = siteBags) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)


siteIndex <- data.frame(siteIndex=unique(data$siteBags),site=1:2)
yearIndex <- data.frame(year=1:3,yearIndex=2006:2008)

interannualG1.2 <- mcmcSamples2 %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(p0_1[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(lambda.med = median(p0_1), 
                   ci.lo = quantile(p0_1,probs=0.27), 
                   ci.hi = quantile(p0_1,probs=0.83),
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(siteIndex,by="site") %>%
  dplyr::left_join(yearIndex,by="year") %>%
  dplyr::select(-c(site,year)) %>%
  dplyr::rename(site = siteIndex) %>%
  dplyr::rename(year = yearIndex) %>%
  dplyr::left_join(vr.site, by="site") %>%
  dplyr::arrange(med) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(desc(lambda.med),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  mutate(id = row_number()) %>% 
  ggplot(aes(x = id , y = lambda.med)) + 
  geom_hline(aes(yintercept=med),linetype='dotted') +
  
  geom_point(aes(color=year)) +
  
  geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +
  
  coord_flip() +
  facet_grid(site ~ ., scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ylab("Probability of germination [P(G)]") +
  ylim(c(0,1))

gridExtra::grid.arrange(interannualG1,interannualG1.2)


v1<-MCMCsummary(mcmcSamples1,params=c("s1.0"))$`50%`
v2<-MCMCsummary(mcmcSamples2,params=c("s1.0"))$`50%`
names <- rownames(MCMCsummary(mcmcSamples1,c("s1.0")))
d<-data.frame(names,v1,v2)
ggplot(d,aes(x=v1,y=v2,color=names)) + geom_abline(intercept=0,slope=1)+ geom_point() + xlim(c(0,1)) + ylim(c(0,1))


denom<-MCMCchains(mcmcSamples2,params="p0_2")[,1]*MCMCchains(mcmcSamples2,params="p0_1")[,1]*(1-MCMCchains(mcmcSamples2,params="p0_2")[,1])

num<-MCMCchains(mcmcSamples2,params="p0_3")[,1]
a<-ifelse(num/denom>1,1,num/denom)
hist(a)

p1<-data %>% dplyr::filter(ageBags==1) %>% ggplot(aes(yearBags,intactOct/seedStart)) + geom_point() + facet_wrap(ageBags~siteBags) + ylim(c(0,1))
p2<-data %>% dplyr::filter(ageBags==2) %>% ggplot(aes(yearBags,totalJan/seedStart)) + geom_point() + facet_wrap(~siteBags) + ylim(c(0,1))
gridExtra::grid.arrange(p1,p2)

ggplot(data) + geom_point(aes(ageBags,intactOct))

d1<-data.frame(data %>% dplyr::filter(ageBags==1) %>% select(siteBags,yearBags,ageBags,intactOct))
d2<-data.frame(data %>% dplyr::filter(ageBags==2) %>% select(siteBags,yearBags,ageBags,totalJan))
d1<-d1 %>% dplyr::rename(y = intactOct)
d2<-d2 %>% dplyr::rename(y = totalJan)
d<-rbind(d1,d2)

ggplot(d) + geom_point(aes(ageBags,y)) + facet_wrap(yearBags~siteBags)


par(mfrow=(c(3,1)))
MCMCchains(mcmcSamples2,params="p0_e")[,1] %>% hist(breaks=100)
MCMCchains(mcmcSamples2,params="p0_4")[,1] %>% hist(breaks=100)
MCMCchains(mcmcSamples2,params="s3.0")[,1] %>% hist(breaks=100)

MCMCchains(mcmcSamples2,params="p0_e")[,2] %>% hist(breaks=100)
MCMCchains(mcmcSamples2,params="p0_4")[,2] %>% hist(breaks=100)
MCMCchains(mcmcSamples2,params="s3.0")[,2] %>% hist(breaks=100)

MCMCchains(mcmcSamples2,params="p0_e")[,4] %>% hist(breaks=100)
MCMCchains(mcmcSamples2,params="p0_4")[,3] %>% hist(breaks=100)
MCMCchains(mcmcSamples2,params="s3.0")[,3] %>% hist(breaks=100)

MCMCchains(mcmcSamples2,params="p0_e")[,5] %>% hist(breaks=100)
MCMCchains(mcmcSamples2,params="p0_4")[,4] %>% hist(breaks=100)
MCMCchains(mcmcSamples2,params="s3.0")[,4] %>% hist(breaks=100)


