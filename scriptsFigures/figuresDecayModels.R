#################################################################################
################################################################################
################################################################################
# Code for figures to ...
#
# 
# Scripts by Gregor Siegmund
# fausto.siegmund@gmail.com
# last updated 02-18-2021
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

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/decayModel/"
modelFittingFiles <- paste0(directory,list.files(directory))

mcmcSamples <- readRDS(modelFittingFiles[[2]])
modelVariables <- readRDS(modelFittingFiles[[1]])
dirFigures = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/"

################################################################################
# Get site names
#################################################################################

directory = "/Users/Gregor/Dropbox/dataLibrary/workflow/tidyData/"
dataFiles <- paste0(directory,list.files(directory))
 
data <- readRDS(dataFiles[[1]])
siteNames = unique(data$siteBags)

################################################################################
# Spatial data
#################################################################################
position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,northing,elevation) %>%
  dplyr::mutate(easting=easting/1000,northing=northing/1000)


################################################################################
# Germination
#################################################################################

g0 <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(mu0_g[siteGermination,gIndex]) %>%
  dplyr::group_by(siteGermination,gIndex) %>%
  dplyr::summarise(med = median(mu0_g), 
                   ci.lo = HPDI(mu0_g,.89)[1], 
                   ci.hi = HPDI(mu0_g,.89)[2],
                   ci.lo2 = HPDI(mu0_g,.5)[1],
                   ci.hi2 = HPDI(mu0_g,.5)[2]
  )

siteIndex <- data.frame(siteIndex=unique(data$siteBags),siteBags=1:20) %>%
  dplyr::rename(siteName=siteIndex,siteGermination=siteBags)
yearIndex <- data.frame(yearIndex = c("2006","2007","2008"), 
                        yearGermination = 1:3 ) %>%
  dplyr::rename(yearName=yearIndex)

g0_indexed=g0 %>%
  dplyr::left_join(siteIndex,by="siteGermination") %>%
  dplyr::ungroup() %>%
  dplyr::select(siteName,gIndex,med) %>%
  dplyr::rename(site = siteName,g0.med = med)


g <-mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(mu_g[siteGermination,yearGermination,gIndex]) %>%
  dplyr::group_by(siteGermination,yearGermination,gIndex) %>%
  dplyr::summarise(med = median(mu_g), 
                   ci.lo = HPDI(mu_g,.89)[1], 
                   ci.hi = HPDI(mu_g,.89)[2],
                   ci.lo2 = HPDI(mu_g,.5)[1],
                   ci.hi2 = HPDI(mu_g,.5)[2]
  )

g_indexed <- g %>%
  dplyr::left_join(siteIndex,by="siteGermination") %>%
  dplyr::left_join(yearIndex,by="yearGermination")

g_ordered = g_indexed %>% dplyr::filter(gIndex==3) %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(siteGermination,yearGermination)) %>%
  dplyr::rename(site = siteName) %>%
  dplyr::rename(year = yearName) %>%
  dplyr::left_join(g0_indexed, by=c("site","gIndex")) %>%
  dplyr::arrange(g0.med) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site,gIndex) %>%
  dplyr::arrange(desc(med),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  mutate(id = row_number())

g_ordered = g_ordered %>%
  dplyr::mutate(med = boot::inv.logit(med), 
                ci.lo = boot::inv.logit(ci.lo),
                ci.hi = boot::inv.logit(ci.hi),
                g0.med = boot::inv.logit(g0.med))



ggplot(data=g_ordered, aes(x = id , y = med)) + 
  geom_point(aes(color=year)) +
  geom_hline(aes(yintercept=g0.med),linetype='dotted') +
  geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +
  coord_flip() +
  facet_grid(site ~ gIndex, scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ylab("Probability of germination [P(G)]") +
  ylim(c(0,1))


# # your data
g0_spatial<-g0 %>%
  dplyr::left_join(siteIndex,by="siteGermination") %>%
  dplyr::ungroup() %>%
  dplyr::rename(site = siteName,g0.med = med) %>%
  dplyr::left_join(position,by="site") %>%
  dplyr::mutate(ci.lo = boot::inv.logit(ci.lo),
                ci.hi = boot::inv.logit(ci.hi),
                g0.med = boot::inv.logit(g0.med))
  # dplyr::mutate(offset=ifelse(site=="EC",-.15,
  #                             ifelse(site=="CP3",.15,0)))

g0_spatial %>%
  ggplot(aes(x = easting , y = g0.med)) + 
  geom_point() +
  geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.5) +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  ylab("Germination probability [P(G)]") +
  xlab("Easting (km)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0),
                     breaks = c(0,.1,.2,.3,.4,.5,.6)) +
  facet_wrap(~gIndex) 
  

################################################################################
# Winter seed survival (s1)
#################################################################################

alpha <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(a[siteSurvival]) %>%
  dplyr::group_by(siteSurvival) %>%
  dplyr::summarise(med = median(a), 
                   ci.lo = HPDI(a,.89)[1], 
                   ci.hi = HPDI(a,.89)[2],
                   ci.lo2 = HPDI(a,.5)[1],
                   ci.hi2 = HPDI(a,.5)[2]
  )

siteIndex <- data.frame(siteIndex=unique(data$siteBags),siteBags=1:20) %>%
  dplyr::rename(siteName=siteIndex,siteSurvival=siteBags)

alpha_indexed=alpha %>%
  dplyr::left_join(siteIndex,by="siteSurvival") %>%
  dplyr::rename(site = siteName)


alpha_spatial<-alpha_indexed %>%
  dplyr::left_join(position,by="site") 

alpha_spatial %>%
  ggplot(aes(x = easting , y = med)) + 
  geom_abline(slope=0,intercept=1,linetype='dotted') +
  geom_point() +
  geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.5) +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  ylab("Weibull shape parameter") +
  xlab("Easting (km)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0,1.5),
                     expand = c(0, 0))# +
  #geom_smooth() 
                    # breaks = c(0,.1,.2,.3,.4,.5,.6)) 

mu0_s <- mcmcSamples %>%
  tidybayes::recover_types(data) %>%
  tidybayes::spread_draws(mu0_s[siteSurvival]) %>%
  dplyr::group_by(siteSurvival) %>%
  dplyr::summarise(med = median(mu0_s), 
                   ci.lo = HPDI(mu0_s,.89)[1], 
                   ci.hi = HPDI(mu0_s,.89)[2],
                   ci.lo2 = HPDI(mu0_s,.5)[1],
                   ci.hi2 = HPDI(mu0_s,.5)[2]
  )

siteIndex <- data.frame(siteIndex=unique(data$siteBags),siteBags=1:20) %>%
  dplyr::rename(siteName=siteIndex,siteSurvival=siteBags)

mu0s_indexed=mu0_s %>%
  dplyr::left_join(siteIndex,by="siteSurvival") %>%
  dplyr::rename(site = siteName)


mu0s_spatial<-mu0s_indexed %>%
  dplyr::left_join(position,by="site") 

mu0s_spatial %>%
  ggplot(aes(x = easting , y = med)) + 
  geom_abline(slope=0,intercept=1,linetype='dotted') +
  geom_point() +
  geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.5) +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  ylab("Germination probability [P(G)]") +
  xlab("Easting (km)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(-1.5,1.5),
                     expand = c(0, 0))# +


dev.off()
mu_s = list()
mean.longevity = list()
for(i in 1:20){
 eta <- mu0s_indexed[i,]$med
 alpha <- alpha_indexed[i,]$med
 inv.b = exp(-eta/alpha)
 mu <- exp(-(seq(0,1,by=.01)/inv.b)^alpha)
 mu_s[[i]] <- mu
# \mathrm{scale} \times \Gamma(1+(\frac{1}{\mathrm{shape}}))
 mean.longevity[[i]]=inv.b*gamma(1+1/alpha)
}
y_arrange=c()
par(mfrow=c(4,5),mar=c(2,1,1,1) + 0.1)

for(i in 1:20){
  plot(seq(0,1,by=.01),mu_s[[i]],type='n',xlim=c(0,1),ylim=c(0,1))
  lines(seq(0,1,by=.01),mu_s[[i]])
  y_arrange[i]=mu_s[[i]][101]
}

plot(position$easting,unlist(mean.longevity))
plot(position$easting[-c(9,15)],unlist(mean.longevity)[-c(9,15)],type='n')
text(position$easting[-c(9,15)],unlist(mean.longevity)[-c(9,15)],
     siteNames[-c(9,15)],cex=.75)
abline(a=coef(lm(unlist(mean.longevity)[-c(9,15)]~position$easting[-c(9,15)]))[1],
       b=coef(lm(unlist(mean.longevity)[-c(9,15)]~position$easting[-c(9,15)]))[2])
data.frame(siteNames,unlist(mean.longevity))

# index=order(y_arrange)
# offset=seq(-.05,.2,length.out=20)
# 
# for(i in 1:20){
# text(x=1.1+offset[i],y=y_arrange[index][i],siteNames[index][i],cex=.5,col='red')
# }


## DISCRETIZE
# for discretizing, may need to include it in the model 
sample.times = c(0,3,12,15,24,27,36)/36

mu_s = list()
for(i in 1:20){
  eta <- mu0s_indexed[i,]$med
  alpha <- alpha_indexed[i,]$med
  inv.b = exp(-eta/alpha)
  mu <- exp(-(sample.times/inv.b)^alpha)
  mu_s[[i]] <- mu
}

plot(sample.times,mu_s[[1]],type='n',xlim=c(0,1.35),ylim=c(0,1))
for(i in 1:20){
  points(sample.times,mu_s[[i]])
}

mu_s.mat = do.call(cbind,mu_s)
surv.prob = list()
for(i in 1:20){
  surv.prob[[i]]=mu_s.mat[c(2:7),i]/mu_s.mat[c(1:6),i]
}

plot(sample.times[1:6],surv.prob[[1]],type='n',xlim=c(0,1.35),ylim=c(0,1))
par(mfrow=c(4,5),mar=c(2,1,1,1) + 0.1)
for(i in 1:20){
  plot(sample.times[1:6],
       surv.prob[[i]],
       xlim=c(0,1),ylim=c(.625,1),
       xlab="",ylab="")
  text(.1,.975,siteNames[i])
}

# s1
par(mfrow=c(2,3))
plot(position$easting,rep(NA,20),type='n',ylim=c(0,1))
y.prob = c()
for(i in 1:20){
  points(position$easting[i],
       surv.prob[[i]][1])
  y.prob[i] = surv.prob[[i]][1]
}
abline(a=coef(lm(y.prob  ~ position$easting))[1],
       b=coef(lm(y.prob  ~ position$easting))[2])

# s2
plot(position$easting,rep(NA,20),type='n',ylim=c(0,1))
for(i in 1:20){
  points(position$easting[i],
         surv.prob[[i]][2])
  y.prob[i] = surv.prob[[i]][2]
}
abline(a=coef(lm(y.prob  ~ position$easting))[1],
       b=coef(lm(y.prob  ~ position$easting))[2])

# s3
plot(position$easting,rep(NA,20),type='n',ylim=c(0,1))
for(i in 1:20){
  points(position$easting[i],
         surv.prob[[i]][3])
  y.prob[i] = surv.prob[[i]][3]
}
abline(a=coef(lm(y.prob  ~ position$easting))[1],
       b=coef(lm(y.prob  ~ position$easting))[2])

# s4
plot(position$easting,rep(NA,20),type='n',ylim=c(0,1))
for(i in 1:20){
  points(position$easting[i],
         surv.prob[[i]][4])
  y.prob[i] = surv.prob[[i]][4]
}
abline(a=coef(lm(y.prob  ~ position$easting))[1],
       b=coef(lm(y.prob  ~ position$easting))[2])

# s5
plot(position$easting,rep(NA,20),type='n',ylim=c(0,1))
for(i in 1:20){
  points(position$easting[i],
         surv.prob[[i]][5])
  y.prob[i] = surv.prob[[i]][5]
}
abline(a=coef(lm(y.prob  ~ position$easting))[1],
       b=coef(lm(y.prob  ~ position$easting))[2])

# s6
plot(position$easting,rep(NA,20),type='n',ylim=c(0,1))
for(i in 1:20){
  points(position$easting[i],
         surv.prob[[i]][6])
  y.prob[i] = surv.prob[[i]][6]
}
abline(a=coef(lm(y.prob  ~ position$easting))[1],
       b=coef(lm(y.prob  ~ position$easting))[2])

# 
# 
# s1_p<-cbind(s1_p,site=siteNames)
# 
# s1_p.med = s1_p %>% dplyr::select(site,med)
# 
# siteIndex <- data.frame(siteIndex=unique(data$siteBags),siteBags=1:20)
# yearIndex <- data.frame(yearBags=1:3,yearIndex=2006:2008)
# 
# s1_py<-mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(s1.0[siteBags,yearBags]) %>%
#   dplyr::group_by(siteBags,yearBags) %>%
#   dplyr::summarise(med = median(s1.0), 
#                    ci.lo = HPDI(s1.0,.89)[1], 
#                    ci.hi = HPDI(s1.0,.89)[2],
#                    ci.lo2 = HPDI(s1.0,.5)[1],
#                    ci.hi2 = HPDI(s1.0,.5)[2]
#   )
# 
# s1_py<-s1_py %>%
#   dplyr::mutate(yearBags = as.integer(yearBags)) %>%
#   dplyr::left_join(yearIndex,by="yearBags") %>%
#   dplyr::ungroup() %>%
#   dplyr::select(-c(yearBags)) %>%
#   dplyr::rename(site = siteBags) %>%
#   dplyr::rename(year = yearIndex) %>%
#   dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)
# 
# 
# siteIndex <- data.frame(siteIndex=unique(data$siteBags),site=1:20)
# yearIndex <- data.frame(year=1:3,yearIndex=2006:2008)
# 
# interannualS1<-mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(s1.0[site,year]) %>%
#   dplyr::group_by(site,year) %>%
#   dplyr::summarise(lambda.med = median(s1.0), 
#                    ci.lo = HPDI(s1.0,.89)[1], 
#                    ci.hi = HPDI(s1.0,.89)[2]
#   ) %>%
#   dplyr::ungroup() %>%
#   dplyr::left_join(siteIndex,by="site") %>%
#   dplyr::left_join(yearIndex,by="year") %>%
#   dplyr::select(-c(site,year)) %>%
#   dplyr::rename(site = siteIndex) %>%
#   dplyr::rename(year = yearIndex) %>%
#   dplyr::left_join(s1_p.med, by="site") %>%
#   dplyr::arrange(med) %>% 
#   dplyr::mutate(site=factor(site,levels=unique(site))) %>%
#   dplyr::group_by(site) %>%
#   dplyr::arrange(desc(lambda.med),.by_group = TRUE) %>%
#   dplyr::mutate(year=factor(year)) %>%
#   mutate(id = row_number()) %>% 
#   ggplot(aes(x = id , y = lambda.med)) + 
#   geom_point(aes(color=year)) +
#   geom_hline(aes(yintercept=med),linetype='dotted') +
#   geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +
#   coord_flip() +
#   facet_grid(site ~ ., scales="free_x", space="free_x") +
#   theme_bw() +
#   theme(panel.spacing=unit(0,"pt"), 
#         panel.border=element_rect(colour="grey50", fill=NA)) +
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank()) +
#   ylab("Probability of seed survival in first winter [P(S1)]") +
#   ylim(c(0,1))
# 
# spatialS1 <-s1_p %>%
#   dplyr::left_join(position,by="site") %>%
#   ggplot(aes(x = easting , y = med)) + 
#   geom_point() +
#   geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +
#   theme_bw() +
#   theme(panel.spacing=unit(0,"pt"), 
#         panel.border=element_rect(colour="grey50", fill=NA)) +
#   ylab("[P(S1)]") +
#   xlab("Easting (km)")+
#   ylim(c(0,1))
# 
# # 
# # ################################################################################
# # # Summer seed survival (s2)
# # #################################################################################
# 
# s2_p <- mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(s2[siteBags]) %>%
#   dplyr::group_by(siteBags) %>%
#   dplyr::summarise(med = median(s2), 
#                    ci.lo = HPDI(s2,.89)[1], 
#                    ci.hi = HPDI(s2,.89)[2],
#                    ci.lo2 = HPDI(s2,.5)[1],
#                    ci.hi2 = HPDI(s2,.5)[2]
#   )
# 
# s2_p<-cbind(s2_p,site=siteNames)
# 
# s2_p.med = s2_p %>% dplyr::select(site,med)
# 
# siteIndex <- data.frame(siteIndex=unique(data$siteBags),siteBags=1:20)
# yearIndex <- data.frame(yearBags=1:3,yearIndex=2006:2008)
# 
# s2_py<-mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(s2.0[siteBags,yearBags]) %>%
#   dplyr::group_by(siteBags,yearBags) %>%
#   dplyr::summarise(med = median(s2.0), 
#                    ci.lo = HPDI(s2.0,.89)[1], 
#                    ci.hi = HPDI(s2.0,.89)[2],
#                    ci.lo2 = HPDI(s2.0,.5)[1],
#                    ci.hi2 = HPDI(s2.0,.5)[2]
#   )
# 
# s2_py<-s2_py %>%
#   dplyr::mutate(yearBags = as.integer(yearBags)) %>%
#   dplyr::left_join(yearIndex,by="yearBags") %>%
#   dplyr::ungroup() %>%
#   dplyr::select(-c(yearBags)) %>%
#   dplyr::rename(site = siteBags) %>%
#   dplyr::rename(year = yearIndex) %>%
#   dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)
# 
# 
# siteIndex <- data.frame(siteIndex=unique(data$siteBags),site=1:20)
# yearIndex <- data.frame(year=1:3,yearIndex=2006:2008)
# 
# interannualS2<-mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(s2.0[site,year]) %>%
#   dplyr::group_by(site,year) %>%
#   dplyr::summarise(lambda.med = median(s2.0), 
#                    ci.lo = HPDI(s2.0,.89)[1], 
#                    ci.hi = HPDI(s2.0,.89)[2]
#   ) %>%
#   dplyr::ungroup() %>%
#   dplyr::left_join(siteIndex,by="site") %>%
#   dplyr::left_join(yearIndex,by="year") %>%
#   dplyr::select(-c(site,year)) %>%
#   dplyr::rename(site = siteIndex) %>%
#   dplyr::rename(year = yearIndex) %>%
#   dplyr::left_join(s2_p.med, by="site") %>%
#   dplyr::arrange(med) %>% 
#   dplyr::mutate(site=factor(site,levels=unique(site))) %>%
#   dplyr::group_by(site) %>%
#   dplyr::arrange(desc(lambda.med),.by_group = TRUE) %>%
#   dplyr::mutate(year=factor(year)) %>%
#   mutate(id = row_number()) %>% 
#   ggplot(aes(x = id , y = lambda.med)) + 
#   geom_point(aes(color=year)) +
#   geom_hline(aes(yintercept=med),linetype='dotted') +
#   geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +
#   coord_flip() +
#   facet_grid(site ~ ., scales="free_x", space="free_x") +
#   theme_bw() +
#   theme(panel.spacing=unit(0,"pt"), 
#         panel.border=element_rect(colour="grey50", fill=NA)) +
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank()) +
#   ylab("Probability of seed survival in first summer [P(S2)]") +
#   ylim(c(0,1))
# 
# spatialS2 <-s2_p %>%
#   dplyr::left_join(position,by="site") %>%
#   ggplot(aes(x = easting , y = med)) + 
#   geom_point() +
#   geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +
#   theme_bw() +
#   theme(panel.spacing=unit(0,"pt"), 
#         panel.border=element_rect(colour="grey50", fill=NA)) +
#   ylab("[P(S2)]") +
#   xlab("Easting (km)") +
#   ylim(c(0,1))
# # ################################################################################
# # # Winter seed survival, year 2 (s3)
# # #################################################################################
# 
# s3_p <- mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(s3[siteBags]) %>%
#   dplyr::group_by(siteBags) %>%
#   dplyr::summarise(med = median(s3), 
#                    ci.lo = HPDI(s3,.89)[1], 
#                    ci.hi = HPDI(s3,.89)[2],
#                    ci.lo2 = HPDI(s3,.5)[1],
#                    ci.hi2 = HPDI(s3,.5)[2]
#   )
# 
# s3_p<-cbind(s3_p,site=siteNames)
# 
# s3_p.med = s3_p %>% dplyr::select(site,med)
# 
# siteIndex <- data.frame(siteIndex=unique(data$siteBags),siteBags=1:20)
# yearIndex <- data.frame(yearBags=1:2,yearIndex=2006:2007)
# 
# s3_py<-mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(s3.0[siteBags,yearBags]) %>%
#   dplyr::group_by(siteBags,yearBags) %>%
#   dplyr::summarise(med = median(s3.0), 
#                    ci.lo = HPDI(s3.0,.89)[1], 
#                    ci.hi = HPDI(s3.0,.89)[2],
#                    ci.lo2 = HPDI(s3.0,.5)[1],
#                    ci.hi2 = HPDI(s3.0,.5)[2]
#   )
# 
# s3_py<-s3_py %>%
#   dplyr::mutate(yearBags = as.integer(yearBags)) %>%
#   dplyr::left_join(yearIndex,by="yearBags") %>%
#   dplyr::ungroup() %>%
#   dplyr::select(-c(yearBags)) %>%
#   dplyr::rename(site = siteBags) %>%
#   dplyr::rename(year = yearIndex) %>%
#   dplyr::select(site,year,med,ci.lo,ci.hi,ci.lo2,ci.hi2)
# 
# 
# siteIndex <- data.frame(siteIndex=unique(data$siteBags),site=1:20)
# yearIndex <- data.frame(year=1:2,yearIndex=2006:2007)
# 
# interannualS3<-mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(s3.0[site,year]) %>%
#   dplyr::group_by(site,year) %>%
#   dplyr::summarise(lambda.med = median(s3.0), 
#                    ci.lo = HPDI(s3.0,.89)[1], 
#                    ci.hi = HPDI(s3.0,.89)[2]
#   ) %>%
#   dplyr::ungroup() %>%
#   dplyr::left_join(siteIndex,by="site") %>%
#   dplyr::left_join(yearIndex,by="year") %>%
#   dplyr::select(-c(site,year)) %>%
#   dplyr::rename(site = siteIndex) %>%
#   dplyr::rename(year = yearIndex) %>%
#   dplyr::left_join(s3_p.med, by="site") %>%
#   dplyr::arrange(med) %>% 
#   dplyr::mutate(site=factor(site,levels=unique(site))) %>%
#   dplyr::group_by(site) %>%
#   dplyr::arrange(desc(lambda.med),.by_group = TRUE) %>%
#   dplyr::mutate(year=factor(year)) %>%
#   mutate(id = row_number()) %>% 
#   ggplot(aes(x = id , y = lambda.med)) + 
#   geom_point(aes(color=year)) +
#   geom_hline(aes(yintercept=med),linetype='dotted') +
#   geom_linerange(aes(x=id,ymin=ci.lo,ymax=ci.hi),size=.25) +
#   coord_flip() +
#   facet_grid(site ~ ., scales="free_x", space="free_x") +
#   theme_bw() +
#   theme(panel.spacing=unit(0,"pt"), 
#         panel.border=element_rect(colour="grey50", fill=NA)) +
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank()) +
#   ylab("Probability of seed survival in first winter [P(S3)]") +
#   ylim(c(0,1))
# 
# spatialS3 <-s3_p %>%
#   dplyr::left_join(position,by="site") %>%
#   ggplot(aes(x = easting , y = med)) + 
#   geom_point() +
#   geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +
#   theme_bw() +
#   theme(panel.spacing=unit(0,"pt"), 
#         panel.border=element_rect(colour="grey50", fill=NA)) +
#   ylab("[P(S3)]") +
#   xlab("Easting (km)") +
#   ylim(c(0,1))
# 
# # combine plots
# 
# 
# 
# interannualS1.j <- interannualS1 + theme(legend.position="none") + scale_fill_manual(values=c("#1b9e77", "#d95f02", "#7570b3"))
# interannualG1.j <- interannualG1 + theme(legend.position="none") + scale_fill_manual(values=c("#1b9e77", "#d95f02", "#7570b3"))
# interannualS2.j <- interannualS2 + theme(legend.position="none") + scale_fill_manual(values=c("#1b9e77", "#d95f02", "#7570b3"))
# interannualS3.j <- interannualS3 + theme(legend.position="none") + scale_fill_manual(values=c("#1b9e77", "#d95f02", "#7570b3"))
# 
# get_legend<-function(myggplot){
#   tmp <- ggplot_gtable(ggplot_build(myggplot))
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#   legend <- tmp$grobs[[leg]]
#   return(legend)
# }
# legend= get_legend(interannualS1 + scale_fill_manual(values=c("#1b9e77", "#d95f02", "#7570b3"),name="Year"))
# 
# interannualPlot <- gridExtra::grid.arrange(interannualS1.j,interannualG1.j,interannualS2.j,interannualS3.j,nrow=1,legend)
# 
# ggsave(filename=paste0(dirFigures,"interannualBelowground.pdf"),
#        plot=interannualPlot,width=20,height=10)
# 
# 
# 
# spatialS1.j <- spatialS1 + theme(legend.position="none") +theme(axis.title.x=element_blank())
# spatialG1.j <- spatialG1 + theme(legend.position="none") +theme(axis.title.x=element_blank())
# spatialS2.j <- spatialS2 + theme(legend.position="none") +theme(axis.title.x=element_blank())
# spatialS3.j <- spatialS3 + theme(legend.position="none")
# 
# spatialPlot <- gridExtra::grid.arrange(spatialS1.j,spatialG1.j,spatialS2.j,spatialS3.j,ncol=1)
# 
# ggsave(filename=paste0(dirFigures,"spatialBelowground.pdf"),
#        plot=spatialPlot,width=5,height=10)
