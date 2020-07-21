# -------------------------------------------------------------------
# Analysis of fitness models
# Outputs reproductive success estimates
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE,max.print=100000)
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

# -------------------------------------------------------------------
# Seedling survival to fruiting
# -------------------------------------------------------------------

# load("~/Dropbox/modelsF2019/output/sigmaFit")
zc <- readRDS("~/Dropbox/dataLibrary/posteriors/seedSurvivalSamples.RDS")

# summary table
# caterpillar plots
par(mfrow=c(1,1))

MCMCplot(zc, params = c("p0_1"), horiz = FALSE,
         sz_labels = .4, sz_med = .75, sz_thin=.7, sz_thick = 2,ylim=c(0,1))
MCMCplot(zc, params = c("p_1"), horiz = FALSE,
         sz_labels = .4, sz_med = .75, sz_thin=.7, sz_thick = 2,ylim=c(0,1))

seedlingSurvivalSummary <-zc %>%
  tidybayes::spread_draws(p_1[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(med = median(p_1), 
                   ci.lo = quantile(p_1,probs=0.025), 
                   ci.hi = quantile(p_1,probs=0.975),
                   ci.lo2 = quantile(p_1,probs=0.25), 
                   ci.hi2 = quantile(p_1,probs=0.75)
  )

# COMMENT OUT BELOW
# a <- MCMCchains(zmP,params="alphaS")
# b <- MCMCchains(zmP,params="betaS")
# g <- MCMCchains(zmP,params="gammaS")
# 
# gamma_names <- dimnames(g)[[2]]
# 
# library(stringr)
# vals<-str_extract_all(gamma_names, "[0-9]+")
# siteGamma<-as.numeric(unlist(lapply(vals, `[[`, 1)))
# yearGamma<-as.numeric(unlist(lapply(vals, `[[`, 2)))
# 
# # site-by-year estimates for survival
# nsite <- dim(a)[2]
# nyear <- dim(b)[2]
# 
# sigmaEstimates<-list()
# sigma<-matrix(NA,nrow=dim(a)[1],ncol=dim(b)[2])
# 
# for(j in 1:nsite){
#   for(k in 1:nyear){
#     alphaIndex <- j
#     betaIndex <- k
#     gammaIndex <- intersect(which(siteGamma %in% j), which(yearGamma %in% k))
# 
#     sigma[,k] = boot::inv.logit(a[,alphaIndex] + b[,betaIndex] + g[,gammaIndex])
#   }
#   sigmaEstimates[[j]] <- sigma
# }

# -------------------------------------------------------------------
# Fruits per plant
# -------------------------------------------------------------------

# load("~/Dropbox/modelsF2019/output/sigmaFit")
zc <- readRDS("~/Dropbox/dataLibrary/posteriors/fruitsPerPlantSamples.RDS")

# summary table
# caterpillar plots
par(mfrow=c(1,1))

MCMCplot(zc, params = c("p0"), horiz = FALSE,
         sz_labels = .4, sz_med = .75, sz_thin=.7, sz_thick = 2,ylim=c(0,80))
MCMCplot(zc, params = c("p"), horiz = FALSE,
         sz_labels = .4, sz_med = .75, sz_thin=.7, sz_thick = 2,ylim=c(0,50))

fruitsPerPlantSummary<-zc %>%
  tidybayes::spread_draws(p[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(med = median(p), 
                   ci.lo = quantile(p,probs=0.025), 
                   ci.hi = quantile(p,probs=0.975),
                   ci.lo2 = quantile(p,probs=0.25), 
                   ci.hi2 = quantile(p,probs=0.75)
  )


# COMMENT BELOW
# a <- MCMCchains(zmP,params="alphaF")
# b <- MCMCchains(zmP,params="betaF")
# g <- MCMCchains(zmP,params="gammaF")
# 
# gamma_names <- dimnames(g)[[2]]
# 
# library(stringr)
# vals<-str_extract_all(gamma_names, "[0-9]+")
# siteGamma<-as.numeric(unlist(lapply(vals, `[[`, 1)))
# yearGamma<-as.numeric(unlist(lapply(vals, `[[`, 2)))
# 
# # site-by-year estimates for survival
# nsite <- dim(a)[2]
# nyear <- dim(b)[2]
# 
# fecEstimates<-list()
# fec<-matrix(NA,nrow=dim(a)[1],ncol=dim(b)[2])
# 
# for(j in 1:nsite){
#   for(k in 1:nyear){
#     alphaIndex <- j
#     betaIndex <- k
#     gammaIndex <- intersect(which(siteGamma %in% j), which(yearGamma %in% k))
# 
#     fec[,k] = exp(a[,alphaIndex] + b[,betaIndex] + g[,gammaIndex])
#   }
#   fecEstimates[[j]] <- fec
# }

# -------------------------------------------------------------------
# Seeds per fruit
# -------------------------------------------------------------------

# load("~/Dropbox/modelsF2019/output/sigmaFit")
zc <- readRDS("~/Dropbox/dataLibrary/posteriors/seedsPerFruitSamples.RDS")

# summary table
# caterpillar plots
par(mfrow=c(1,1))

MCMCplot(zc, params = c("p0"), horiz = FALSE,
         sz_labels = .4, sz_med = .75, sz_thin=.7, sz_thick = 2,ylim=c(0,100))
MCMCplot(zc, params = c("p"), horiz = FALSE,
         sz_labels = .4, sz_med = .75, sz_thin=.7, sz_thick = 2,ylim=c(0,100))

seedsPerFruitSummary<-zc %>%
  tidybayes::spread_draws(p[site,year]) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(med = median(p), 
                   ci.lo = quantile(p,probs=0.025), 
                   ci.hi = quantile(p,probs=0.975),
                   ci.lo2 = quantile(p,probs=0.25), 
                   ci.hi2 = quantile(p,probs=0.75)
  )

# COMMENT BELOW
# a <- MCMCchains(zmP,params="alphaP")
# b <- MCMCchains(zmP,params="betaP")
# g <- MCMCchains(zmP,params="gammaP")
# 
# gamma_names <- dimnames(g)[[2]]
# 
# library(stringr)
# vals<-str_extract_all(gamma_names, "[0-9]+")
# siteGamma<-as.numeric(unlist(lapply(vals, `[[`, 1)))
# yearGamma<-as.numeric(unlist(lapply(vals, `[[`, 2)))
# 
# # site-by-year estimates for survival
# nsite <- dim(a)[2]
# nyear <- dim(b)[2]
# 
# phiEstimates<-list()
# phi<-matrix(NA,nrow=dim(a)[1],ncol=dim(b)[2])
# 
# for(j in 1:nsite){
#   for(k in 1:nyear){
#     alphaIndex <- j
#     betaIndex <- k
#     gammaIndex <- intersect(which(siteGamma %in% j), which(yearGamma %in% k))
# 
#     phi[,k] = exp(a[,alphaIndex] + b[,betaIndex] + g[,gammaIndex])
#   }
#   phiEstimates[[j]] <- phi
# }

# -------------------------------------------------------------------
# Reproductive success
# -------------------------------------------------------------------
nsite <- 20
nyear <- 6

sigmaSummary <- seedlingSurvivalSummary %>%
  dplyr::select(site,year,med) %>%
  dplyr::rename(sigma = med)

fecSummary <- fruitsPerPlantSummary %>%
  dplyr::select(site,year,med) %>%
  dplyr::rename(fec = med)

phiSummary <- seedsPerFruitSummary %>%
  dplyr::select(site,year,med) %>%
  dplyr::rename(phi = med)

abovegroundMedians <- sigmaSummary %>% 
  dplyr::full_join(fecSummary,by=c("site","year")) %>%
  dplyr::full_join(phiSummary,by=c("site","year"))

rsEstimates <- abovegroundMedians %>% dplyr::mutate(rs = sigma*fec*phi)

# save estimates of RS with sites/years where there is missing data excluded
saveRDS(rsEstimates,file="~/Dropbox/clarkiaSeedBanks/products/dataFiles2/rsVarPosterior.RDS")


## rsEstimates

rm(list=ls(all=TRUE)) # clear R environment

# recover p0_1 using tidybayes
# write for loop filtering to site
# write for loop filtering to year
# for each site-year, multiply all draws/iterations
# save to a 


# load("~/Dropbox/modelsF2019/output/sigmaFit")
zc <- readRDS("~/Dropbox/dataLibrary/posteriors/seedSurvivalSamples.RDS")

mcmcSigma <-zc %>%
  tidybayes::spread_draws(p_1[site,year]) 

zc <- readRDS("~/Dropbox/dataLibrary/posteriors/fruitsPerPlantSamples.RDS")

mcmcFec <-zc %>%
  tidybayes::spread_draws(p[site,year])

zc <- readRDS("~/Dropbox/dataLibrary/posteriors/seedsPerFruitSamples.RDS")

mcmcPhi <-zc %>%
  tidybayes::spread_draws(p[site,year])

f<-function(chains=zc,site=i,year=j){
  tmp <- chains %>% dplyr::filter(site==i&year==j)
  return(tmp)
}

rsEstimates<-list()
iter<-max(dim(mcmcSigma)[1],dim(mcmcFec)[1],dim(mcmcPhi)[1])
rs<-matrix(NA,nrow=3000,ncol=7)

# change if number of years changes
for(i in 1:20){
  for(j in 1:7){
    rs[,j] <- f(chains=mcmcSigma,site=i,year=j)$p_1*f(chains=mcmcFec,site=i,year=j)$p*f(chains=mcmcPhi,site=i,year=j)$p
  }
  rsEstimates[[i]] <- rs
}

saveRDS(rsEstimates,file="~/Dropbox/clarkiaSeedBanks/products/dataFiles2/rsVarFullPosterior.RDS")


# rsEstimates<-list()
# rs<-matrix(NA,nrow=30000,ncol=nyear)
# 
# for(j in 1:nsite){
#   s <- sigmaEstimates[[j]]
#   f <- fecEstimates[[j]]
#   p <- phiEstimates[[j]]
# 
#   for(k in 1:nyear){
# 
#     rs[,k]<-s[,k]*f[,k]*p[,k]
#   }
#   rsEstimates[[j]] <- rs
# }

# # load file with index site*year oof missing data
# load("~/Dropbox/modelsF2019/output/missingness")
# 
# for(i in 1:(dim(missing_dat)[1])){
#   v<-as.numeric(missing_dat[i,])
#   rsEstimates[[v[1]]]<-rsEstimates[[(v[1])]][,-(v[2])]
# }
# 
# g1 <- mcmcSamples %>%
#   tidybayes::recover_types(data) %>%
#   tidybayes::spread_draws(p[site,year]) %>%
#   dplyr::group_by(site,year) %>%
#   dplyr::summarise(lambda.med = median(p), 
#                    ci.lo = quantile(p,probs=0.27), 
#                    ci.hi = quantile(p,probs=0.83),
#   ) %>%
#   
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
#         axis.ticks.y=element_blank())
# 
# ggsave(filename=paste0(dirFigures,"interannual-fruitsPerPlant.pdf"),
#        plot=g1,width=6,height=12)