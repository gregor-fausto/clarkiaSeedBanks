#################################################################################
################################################################################
################################################################################
# Code for figures to compare the following modeling approaches for the seedling survivorship data
# 
# Scripts by Gregor Siegmund
# fausto.siegmund@gmail.com
# last updated 03-07-2021
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

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedlingSurvival/"
simFiles <- paste0(directory,list.files(directory))

mcmcSamples <- readRDS(simFiles[[3]])

dirPars = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/parameterSummary/"
dir.create(file.path(dirPars), showWarnings = FALSE)

################################################################################
# Data directory
#################################################################################

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedlingSurvival/"
dataFiles <- paste0(directory,list.files(directory))

data <- readRDS(dataFiles[[2]])
censusSeedlingsFruitingPlants <- readRDS(dataFiles[[1]])

# -------------------------------------------------------------------
# Get site names and position
# -------------------------------------------------------------------

directory2 = "/Users/Gregor/Dropbox/dataLibrary/workflow/tidyData/"
dataFiles <- paste0(directory2,list.files(directory2))

data2 <- readRDS(dataFiles[[1]])
siteNames = unique(data2$siteBags)
years = as.numeric(unique(censusSeedlingsFruitingPlants$year))
yearNames = years[order(years)]

position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,dominant.surface.rock.type) %>%
  dplyr::mutate(easting=easting/1000)

################################################################################
# Summarize seedlings survival parameters in CSV
#################################################################################

siteIndex <- data.frame(siteIndex=unique(censusSeedlingsFruitingPlants$site),site=1:20)
yearIndex <- data.frame(yearIndex=as.numeric(unique(censusSeedlingsFruitingPlants$year)),
                        year=unique(unique(data$year)))

summary.fun = function(x){
  quant=quantile(x,c(.025,.975))
  mu = mean(x)
  med = median(x)
  hpdi.interval = HPDI(x, prob = .95)
  return(c(mu,quant,med,hpdi.interval))
}

parm.mu=MCMCchains(mcmcSamples,params='mu')
parm.prob=apply(parm.mu,2,boot::inv.logit)
parm.prob.sum = apply(parm.prob,2,summary.fun)

df.list = list()
for(i in 1:20){
  obj = parm.prob.sum
  index=grep(paste0("\\[",i,","),colnames(obj))
  tmp = signif(parm.prob.sum[,index],3)
  tmp.df=data.frame(site=siteNames[i],year=yearNames,(t(tmp)))
  names(tmp.df) = c("site","year","mu","ci.lo95","ci.hi95","med","hpdi.lo95","hpdi.hi95")
  rownames(tmp.df) = NULL
  df.list[[i]] = tmp.df
}

summary.df=do.call(rbind,df.list)
# 
# sigma_py<-mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(mu[site,year]) %>%
#   dplyr::mutate(p_1 = boot::inv.logit(mu)) %>%
#   dplyr::group_by(site,year) %>%
#   dplyr::summarise(mu = mean(p_1),
#                    ci.lo95 = quantile(p_1,probs=0.025), 
#                    ci.hi95 = quantile(p_1,probs=0.975),
#                    med = median(p_1), 
#                    hpdi.lo95 = HPDI(p_1,.95)[1], 
#                    hpdi.hi95 = HPDI(p_1,.95)[2]
#   )
# 
# sigma_py<-sigma_py %>%
#   dplyr::left_join(siteIndex,by="site") %>%
#   dplyr::left_join(yearIndex,by="year") %>%
#   dplyr::ungroup() %>%
#   dplyr::select(-c(site,year)) %>%
#   dplyr::rename(site = siteIndex) %>%
#   dplyr::rename(year = yearIndex) %>%
#   dplyr::select(site,year,mu,ci.lo95,ci.hi95,med,hpdi.lo95,hpdi.hi95)

# saveRDS(sigma_py,file=paste0(dirPars,"sigmaSummary.RDS"))
sigma_py = summary.df 
write.csv(sigma_py,file=paste0(dirPars,"sigmaSummary.csv"))


################################################################################
# Make summary plots
#################################################################################


summary.fun = function(x){
  med = median(x)
  hpdi.interval = HPDI(x, prob = .89)
  hpdi.interval2 = HPDI(x, prob = .5)
  return(c(med,hpdi.interval,hpdi.interval2))
}

parm.mu0=MCMCchains(mcmcSamples,params='mu0')
parm.prob0=apply(parm.mu0,2,boot::inv.logit)
parm.prob0.sum = apply(parm.prob0,2,summary.fun)

df.list = list()
for(i in 1:20){
  obj = parm.prob0.sum
  index=grep(paste0("\\[",i,"\\]"),colnames(obj))
  tmp = signif(parm.prob0.sum[,index],3)
  tmp.df=data.frame(site=siteNames[i],(t(tmp)))
  names(tmp.df) = c("site","med","ci.lo","ci.hi","ci.lo2","ci.hi2")
  rownames(tmp.df) = NULL
  df.list[[i]] = tmp.df
}

summary.pop.df=do.call(rbind,df.list)


parm.prob=apply(parm.mu,2,boot::inv.logit)
parm.prob.sum = apply(parm.prob,2,summary.fun)

df.list = list()
for(i in 1:20){
  obj = parm.prob.sum
  index=grep(paste0("\\[",i,","),colnames(obj))
  tmp = signif(parm.prob.sum[,index],3)
  tmp.df=data.frame(site=siteNames[i],year=yearNames,(t(tmp)))
  names(tmp.df) = c("site","year","med","ci.lo","ci.hi","ci.lo2","ci.hi2")
  rownames(tmp.df) = NULL
  df.list[[i]] = tmp.df
}

summary.df=do.call(rbind,df.list)


# rename
sigma_p = summary.pop.df %>% dplyr::select(site,med)
sigma_py = summary.df 

tmp = sigma_py %>%
  dplyr::rename(lambda.med = med) %>%
  dplyr::left_join(sigma_p,by = "site")


interannualSigmaDF =  tmp %>%
  dplyr::arrange(med) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(desc(lambda.med),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  mutate(id = row_number())

interannualSigma <- interannualSigmaDF %>%
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
  ylab("Probability of seedling survival to fruiting") +
  ylim(c(0,1)) +
  labs(color="Year") 

sigma_p = summary.pop.df

spatialSigma <- sigma_p %>%
  dplyr::left_join(position,by="site") %>%
  ggplot(aes(x = easting , y = med)) + 
  geom_point() +
  geom_linerange(aes(x=easting,ymin=ci.lo,ymax=ci.hi),size=.25) +
  theme_bw() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  ylab("Probability of seedling survival to fruiting") +
  xlab("Easting (km)") +
  ylim(c(0,1))

dirFigures = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/"

ggsave(filename=paste0(dirFigures,"interannualSigma.pdf"),
       plot=interannualSigma,width=6,height=12)

ggsave(filename=paste0(dirFigures,"spatialSigma.pdf"),
       plot=spatialSigma,width=6,height=12)

