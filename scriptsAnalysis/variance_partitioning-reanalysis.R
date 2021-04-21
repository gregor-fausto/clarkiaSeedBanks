# Read in library and data ------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE,max.print=100000)
# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(rjags) # jags interface
library(MCMCvis)
library(tidyverse)
library(reshape2)
library(HDInterval)
library(bayesplot)
library(rethinking)

# -------------------------------------------------------------------
# Read in samples from posterior distributions
# -------------------------------------------------------------------

sigma <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/sigma-analysis.RDS")
fec <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/tfe-analysis.RDS")
phi <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/phi-analysis.RDS")

rsPosterior = sigma*fec*phi

posterior.mode = function(x){
  x.max=max(x)
  x.min=min(x)
  dres <- density( x ,from = x.min, to = x.max)
  modeParam <- dres$x[which.max(dres$y)]
  return(modeParam)
}

# get mode from full posterior and calculate var in RS based on modes
sigma.mode = apply(sigma,2,median);
fec.mode = apply(fec,2,median);
phi.mode = apply(phi,2,median);

rs.mode=apply(sigma*fec*phi,2,median)
dat <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/dataFiles/rsMedianEstimates.RDS")

siteNames = unique(dat$site)

df.list<-list()
for(i in 1:20){
  index=grep(paste0("\\[",i,","),names(rs.mode))
  df=data.frame(site=rep(siteNames[i],15),rs=rs.mode[index],sigma=sigma.mode[index],fec=fec.mode[index],phi=phi.mode[index])
  df.list[[i]]=df
}
rsMedians = do.call(rbind,df.list)

library(tidyverse)
library(ggrepel)


# Validate geometric SD calculation ---------------------------------------


# https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) /( length(x)))
}

A <- rsMedians[rsMedians$site=="BG",]$rs
n = length(A)
mu_g = (prod(A))^(1/n)

exp(sqrt(sum((log(A/mu_g))^2)/n))

# sample implementation
gsd.sample <- function(x){
  n = length(x[!is.na(x)])
  mu = (prod(x,na.rm=TRUE))^(1/n)
  y <- exp(sqrt(sum((log(x/mu))^2,na.rm=TRUE)/(n-1)))
  return(y)
}

gsd.am <- function(x){
  n = length(x[!is.na(x)])
  mu = exp(mean(log(x),na.rm=TRUE))
  y <- exp(sqrt(sum((log(x/mu))^2,na.rm=TRUE)/(n-1)))
  return(y)
}

# lognormal implementation
gsd <- function(x){
  y <- exp(sd(log(x),na.rm=TRUE))
  return(y)
}

sd(log(A))
sqrt(sum((log(A/mu_g))^2)/n)

gsd(A)

# Test geometric standard deviation function ------------------------------------------------------

x <- c(rlnorm(n=10000,meanlog = 0 , sdlog = 1),NA)
y <- c(rlnorm())
gsd(x)
gsd.sample(x)
gsd.am(x)

x <- rnorm(n=1000,mean = 20 , sd = 1)
gsd(x)
gsd.am(x)


# Compare GSD calculations ------------------------------------------------

gsd_comparison <- rsMedians %>% 
  dplyr::group_by(site) %>%
  dplyr::summarise(gsd.log=gsd(rs),
                   gsd.am=gsd.am(rs))

plot(gsd_comparison$gsd.log,gsd_comparison$gsd.am)
abline(a=0,b=1)

ggplot(gsd_comparison,aes(x=gsd.log,y=gsd.am,label=site)) +
  geom_abline(intercept=0,slope=1,color='lightgray') +
  geom_point() +
  geom_text_repel(size=3,color="black") +
  theme_bw() +
  xlab("Geometric SD of fitness (log-normal theory)") +
  ylab("Geometric SD of fitness (sample)")# +
#  xlim(c(2,7) ) + ylim(c(2,7)) 
  

# Geometric variance and covariance ------------------------------------------------------

gvar.am <- function(x){
  n = length(!is.na(x))
  mu = exp(mean(log(x),na.rm=TRUE))
  y <- exp(sum((log(x/mu))^2,na.rm=TRUE)/(n-1))
  return(y)
}

gvar.sample <- function(x){
  n = length(!is.na(x))
  mu = (prod(x,na.rm=TRUE))^(1/m)
  y <- exp((sum((log(x/mu))^2)/(n-1)))
  return(y)
}

# note x and y *must* be same length here
# not sure of other
# which would be the equivalent of
# use = "complete.obs" in cov()
gcovar.am <- function(x,y){
  n_x = length(!is.na(x))
  n_y = length(!is.na(x))
  #n = ifelse(n_x>n_y,n_y,n_x)
  n = max(n_x,n_y)
  
  mu_x = exp(mean(log(x),na.rm=TRUE))
  mu_y = exp(mean(log(y),na.rm=TRUE))
  
  tmp <- exp(sum(((log(x)-log(mu_x))*(log(y)-log(mu_y))),na.rm=TRUE)/(n-1))
  return(tmp)
}

gcovar.sample <- function(x,y){
  n_x = length(!is.na(x))
  n_y = length(!is.na(x))
  n = ifelse(n_x>n_y,n_y,n_x)
  
  mu_x = (prod(x,na.rm=TRUE))^(1/n)
  mu_y = (prod(y,na.rm=TRUE))^(1/n)
  
  tmp <- exp(sum(((log(x)-log(mu_x))*(log(y)-log(mu_y))),na.rm=TRUE)/(n-1))
  return(tmp)
}

# Test geometric variance and covariance function ------------------------------------------------------

x <- c(rlnorm(n=10000,meanlog = 0 , sdlog = 1))
y <- c(rlnorm(n=10000,meanlog = 0 , sdlog = 2))
z <- x*y

gsd.am(x);gsd.am(y)
gvar.am(x);gvar.am(y)
gcovar.am(x,y) 

gvar.am(z)
gvar.am(x)*gvar.am(y)*(gcovar.am(x,y)^2)

# gsd = gvar = e^0 ~ 2.7 => variance=sd=1 (on log scale)
# gcovar = 1 => not covarying

# Apply to data -----------------------------------------------------------

# testing rsMedians <- rsMedians %>% dplyr::filter(site=="BR") 

# correct when excluding NAs
rsMedians<-rsMedians %>% dplyr::filter(!is.na(rs))

var_sigma = rsMedians %>% 
  dplyr::filter(!is.na(rs)) %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(var_sigma = gvar.am(sigma))

var_fec = rsMedians %>% 
  dplyr::filter(!is.na(rsMedians)) %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(var_fec = gvar.am(fec))

var_phi = rsMedians %>% 
  dplyr::filter(!is.na(rs)) %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(var_phi = gvar.am(phi))

cov_sigma.fec = rsMedians %>%
  dplyr::filter(!is.na(rs)) %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(cov_sigma.fec = gcovar.am(sigma,fec))

cov_sigma.phi = rsMedians %>%
  dplyr::filter(!is.na(rs)) %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(cov_sigma.phi = gcovar.am(sigma,phi))

cov_fec.phi = rsMedians %>%
  dplyr::filter(!is.na(rs)) %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(cov_fec.phi = gcovar.am(fec,phi))

# Sum variance components -------------------------------------------------

varianceComponents <-var_sigma %>%
  dplyr::left_join(var_fec,by="site") %>%
  dplyr::left_join(var_phi,by="site") %>%
  dplyr::left_join(cov_sigma.fec,by="site") %>%
  dplyr::left_join(cov_sigma.phi,by="site") %>%
  dplyr::left_join(cov_fec.phi,by="site")

# Sum variance components -------------------------------------------------

varianceComponents <- varianceComponents %>%
  dplyr::mutate(vars.only = var_sigma*var_fec*var_phi) %>%
  dplyr::mutate(covars.only = (cov_sigma.fec^2)*(cov_sigma.phi^2)*(cov_fec.phi^2)) %>%

  dplyr::mutate(var.total = var_sigma*var_fec*var_phi*(cov_sigma.fec^2)*(cov_sigma.phi^2)*(cov_fec.phi^2))

# Calculate variance from RS -------------------------------------------------

var_rs = rsMedians %>% 
  dplyr::group_by(site) %>%
  dplyr::summarise(var_rs = gvar.am(rs))

# Plotting -------------------------------------------------
# note that these are close but not identical. why is there a difference?

plot(var_rs$var_rs,varianceComponents$var.total);abline(a=0,b=1)

varianceComponents = varianceComponents %>%
  dplyr::left_join(var_rs,by="site") 

comparisonPlot <- ggplot(varianceComponents,aes(x=var_rs,y=var.total,label=site)) +
  geom_abline(intercept=0,slope=1,color='lightgray') +
  geom_point() +
  geom_text_repel(size=3,color="black") +
  theme_bw() +
  xlab("Geometric var of fitness\n calculated from RS") +
  ylab("Geometric var of fitness\n calculated from components") 

# dirFigures = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/appendix/varianceDecomp/"
# 
# ggsave(filename=paste0(dirFigures,"comparison.pdf"),
#        plot=comparisonPlot,width=6,height=6)

# the reason the variance of those sites is higher from the sum of components than from
# the individual variance is because of the missing data; there are years with data on sigma
# but no others; these are thrown out for covariance estimation but not for variance
# not sure how to deal with this

# variance only
varianceComponents <- varianceComponents %>% 
  dplyr::mutate(corrected=ifelse(var.total-var_rs>.1,as.character(0),as.character(1)))


variancePlot <- varianceComponents %>%
  dplyr::mutate(site=as.factor(site)) %>%
  tidyr::pivot_longer(cols=c(var_sigma,var_fec,var_phi)) %>%
  dplyr::mutate(name=ifelse(name=="var_sigma","Var(sigma)",name),
                name=ifelse(name=="var_fec","Var(fec)",name),
                name=ifelse(name=="var_phi","Var(phi)",name)) %>%
  dplyr::mutate(name = factor(name, 
                              levels = c("Var(sigma)", "Var(fec)", "Var(phi)"))) %>% 
  ggplot(aes(x=site,y=value)) +
  geom_hline(yintercept=c(exp(0)),linetype='dotted') +
  geom_point() +
  scale_shape_manual(values = c(16)) +
  theme_classic() +
  xlab("Site") + ylab("Variance") +
  facet_wrap(~name) +
 # ylim(c(0.5,7.5))  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

covariancePlot <- varianceComponents %>%
  dplyr::mutate(site=as.factor(site)) %>%
  tidyr::pivot_longer(cols=c(cov_sigma.fec,cov_sigma.phi,cov_fec.phi)) %>%
  dplyr::mutate(name=ifelse(name=="cov_sigma.fec","Cov(sigma,fec)",name),
                name=ifelse(name=="cov_sigma.phi","Cov(sigma,phi)",name),
                name=ifelse(name=="cov_fec.phi","Cov(fec,phi)",name)) %>%
  dplyr::mutate(name = factor(name, 
                              levels = c(
                                         "Cov(sigma,fec)", "Cov(sigma,phi)",
                                         "Cov(fec,phi)"))) %>%  
  ggplot(aes(x=site,y=value)) +
  geom_hline(yintercept=c(exp(0),exp(1)),linetype='dotted') +
  geom_point() +
  scale_shape_manual(values = c(16)) +
  theme_classic() +
  xlab("Site") + ylab("Covariance") +
  facet_wrap(~name) +
  ylim(c(0.5,3))

  # gsd = gvar = e^0 ~ 2.7 => variance=sd=1 (on log scale)
  # gcovar = 1 => not covarying  
  # var/cov <1 dampens

gridExtra::grid.arrange(variancePlot,covariancePlot,nrow=2)
  
dirFigures = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/analysis/"

ggsave(filename=paste0(dirFigures,"variance-decomp.pdf"),
       plot=variancePlot,width=8,height=6)

# free axis

varianceComponents %>%
  dplyr::mutate(site=as.factor(site)) %>%
  tidyr::pivot_longer(cols=c(var_sigma,var_fec,var_phi)) %>%
  dplyr::mutate(name=ifelse(name=="var_sigma","Var(sigma)",name),
                name=ifelse(name=="var_fec","Var(fec)",name),
                name=ifelse(name=="var_phi","Var(phi)",name)) %>%
  dplyr::mutate(name = factor(name, 
                              levels = c("Var(sigma)", "Var(fec)", "Var(phi)"))) %>% 
  dplyr::group_by(site) %>%
  dplyr::arrange(var.total) %>%
  ggplot(aes(x=value,y=var.total)) +
   geom_vline(xintercept=c(exp(0)),linetype='dotted') +
  geom_point() +
  theme_classic() +
  xlab("Component variance") + ylab("Total variance") +
  facet_wrap(~name,scales='free')

varianceComponents %>%
  dplyr::mutate(site=as.factor(site)) %>%
  tidyr::pivot_longer(cols=c(cov_sigma.fec,cov_sigma.phi,cov_fec.phi)) %>%
  dplyr::mutate(name=ifelse(name=="cov_sigma.fec","Cov(sigma,fec)",name),
                name=ifelse(name=="cov_sigma.phi","Cov(sigma,phi)",name),
                name=ifelse(name=="cov_fec.phi","Cov(fec,phi)",name)) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(var.total) %>%
  ggplot(aes(x=value,y=var.total)) +
  geom_vline(xintercept=c(exp(0)),linetype='dotted') +
  geom_point() +
  theme_classic() +
  xlab("Component covariance") + ylab("Total variance") +
  facet_wrap(~name,scales='free')


df <-   varianceComponents %>%
  dplyr::mutate(site=as.factor(site)) %>%
  dplyr::mutate(cov_sigma.fec=cov_sigma.fec^2,
                cov_sigma.phi=cov_sigma.phi^2,
                cov_fec.phi=cov_fec.phi^2) %>%
  tidyr::pivot_longer(cols=c(cov_sigma.fec,cov_sigma.phi,cov_fec.phi,var_sigma,var_fec,var_phi)) %>%
  dplyr::mutate(name=ifelse(name=="cov_sigma.fec","Cov(survival,fruits)",name),
                name=ifelse(name=="cov_sigma.phi","Cov(survival,seeds)",name),
                name=ifelse(name=="cov_fec.phi","Cov(fruits,seeds)",name),
                name=ifelse(name=="var_sigma","Var(survival)",name),
              name=ifelse(name=="var_fec","Var(fruits)",name),
              name=ifelse(name=="var_phi","Var(seeds)",name)) %>%
  dplyr::group_by(site) %>%
  dplyr::mutate(name = factor(name, 
                                    levels = c("Var(survival)", "Var(fruits)", "Var(seeds)",
                                               "Cov(survival,fruits)", "Cov(survival,seeds)",
                                               "Cov(fruits,seeds)"))) 

allArrayed <- df%>% 
  ggplot(aes(x=name,y=value,group=site)) +
  geom_hline(yintercept=c(exp(0)),linetype='dotted',size=.5) +
  geom_point(size=1) +
  scale_shape_manual(values = c(16)) +
  geom_line(alpha=0.5) +

  #scale_linetype_manual(values = c("dotted", "solid")) +
  theme_classic() +
  xlab("Component") + ylab("Variance or covariance") +
  theme(text = element_text(size=20))
  

ggsave(filename=paste0(dirFigures,"lineplot.pdf"),
       plot=allArrayed,width=14,height=6)

variancePlot

barPlots<-varianceComponents %>%
  dplyr::mutate(site=as.factor(site)) %>%
  tidyr::pivot_longer(cols=c(cov_sigma.fec,cov_sigma.phi,cov_fec.phi,var_sigma,var_fec,var_phi)) %>%
  dplyr::mutate(name=ifelse(name=="cov_sigma.fec","Cov(sigma,fec)",name),
                name=ifelse(name=="cov_sigma.phi","Cov(sigma,phi)",name),
                name=ifelse(name=="cov_fec.phi","Cov(fec,phi)",name),
                name=ifelse(name=="var_sigma","Var(sigma)",name),
                name=ifelse(name=="var_fec","Var(fec)",name),
                name=ifelse(name=="var_phi","Var(phi)",name)) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(var.total,.by_group=TRUE) %>%
  
  mutate(id = cur_group_id()) %>%
  dplyr::mutate(name = factor(name, 
                              levels = c("Var(sigma)", "Var(fec)", "Var(phi)",
                                         "Cov(sigma,fec)", "Cov(sigma,phi)",
                                         "Cov(fec,phi)"), ordered=TRUE)) %>%
  
  ggplot(aes(x=name,y=value,group=id)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("lightgray")) +
  geom_hline(yintercept=c(exp(0)),linetype='solid',color='red',size=1) +
  facet_wrap(~site,scales='free',nrow=2) +
  theme_classic() +
  xlab("Component") + ylab("Component variance or covariance") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(filename=paste0(dirFigures,"barPlots.pdf"),
       plot=barPlots,width=10,height=6)



varianceVsEasting<-varianceComponents %>%
  dplyr::mutate(site=as.factor(site)) %>%
  tidyr::pivot_longer(cols=c(cov_sigma.fec,cov_sigma.phi,cov_fec.phi,var_sigma,var_fec,var_phi)) %>%
  dplyr::mutate(name=ifelse(name=="cov_sigma.fec","Cov(sigma,fec)",name),
                name=ifelse(name=="cov_sigma.phi","Cov(sigma,phi)",name),
                name=ifelse(name=="cov_fec.phi","Cov(fec,phi)",name),
                name=ifelse(name=="var_sigma","Var(sigma)",name),
                name=ifelse(name=="var_fec","Var(fec)",name),
                name=ifelse(name=="var_phi","Var(phi)",name)) %>%
  dplyr::arrange(desc(var.total),.by_group=TRUE) %>%
  dplyr::mutate(name = factor(name, 
                              levels = c("Var(sigma)", "Var(fec)", "Var(phi)",
                                         "Cov(sigma,fec)", "Cov(sigma,phi)",
                                         "Cov(fec,phi)"), ordered=TRUE)) %>%
  dplyr::left_join(position,by="site") %>%
  dplyr::mutate(perc=value/var.total) %>% 
  
  ggplot(aes(x=easting,y=value,label=site)) +
  geom_hline(yintercept=c(exp(0)),linetype='dotted',color='lightgray',size=.5) +
    geom_point() +
  facet_wrap(~name,scales='free',nrow=2) +
  theme_classic() +
  xlab("Easting") + ylab("Component variance or covariance") +
  geom_text_repel(size=1.5,color="black") 
  
ggsave(filename=paste0(dirFigures,"varianceVsEasting.pdf"),
       plot=varianceVsEasting,width=10,height=6)

varianceComponents %>%
  dplyr::mutate(site=as.factor(site)) %>%
  tidyr::pivot_longer(cols=c(cov_sigma.fec,cov_sigma.phi,cov_fec.phi,var_sigma,var_fec,var_phi)) %>%
  dplyr::mutate(name=ifelse(name=="cov_sigma.fec","Cov(sigma,fec)",name),
                name=ifelse(name=="cov_sigma.phi","Cov(sigma,phi)",name),
                name=ifelse(name=="cov_fec.phi","Cov(fec,phi)",name),
                name=ifelse(name=="var_sigma","Var(sigma)",name),
                name=ifelse(name=="var_fec","Var(fec)",name),
                name=ifelse(name=="var_phi","Var(phi)",name)) %>%
  dplyr::mutate(name = factor(name, 
                              levels = rev(c("Var(sigma)", "Var(fec)", "Var(phi)",
                                         "Cov(sigma,fec)", "Cov(sigma,phi)",
                                         "Cov(fec,phi)")), ordered=TRUE)) %>%
  dplyr::mutate(perc=value/var.total) %>% 
  ggplot(aes(x=name,y=value,group=site)) +
  geom_segment( aes(x=name ,xend=name, y=0, yend=value), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme_classic() +
  facet_wrap(~site,scales='free_x') +
  xlab("Component variance or covariance") + ylab("")  +
  theme()

position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,northing,elevation) %>%
  dplyr::mutate(easting=easting/1000,northing=northing/1000)

library(MASS)

parcoord(varianceComponents[,c(2,3,4)] )

# PCA ---------------------------------------------------------------------

varianceComponentsPCA <- varianceComponents %>%
  dplyr::select(var_sigma, var_fec, var_phi,
                cov_sigma.fec,cov_sigma.phi,cov_fec.phi)

varcomp.pca <- prcomp(varianceComponentsPCA)

pdf(file=paste0(dirFigures,"pcaPlots.pdf"),width=10,height=6)

plot(varcomp.pca,type='l')

# plot(varcomp.pca$x[,1], varcomp.pca$x[,2], pch='', xlab = "PC 1", ylab = "PC 2",
#      xlim=c(-5,5),ylim=c(-1,1))
# text(varcomp.pca$x[,1], varcomp.pca$x[,2], labels=varianceComponents$site)

pcMat = data.frame(varcomp.pca$x,site=varianceComponents$site)

ggplot(data = pcMat,aes(x=PC1,y=PC2,label=site)) +
  geom_point() +
  geom_text_repel(size=3,color="black") +
  theme_bw() +
  ylim(c(-.5,1)) + xlim(c(-4,4))


fviz_pca_var(varcomp.pca, col.var = "black",
             axes = c(1,2)) + xlim(c(-2.1,2.1))

fviz_contrib(varcomp.pca, choice = "var", axes = 1, top = 10)
fviz_contrib(varcomp.pca, choice = "var", axes = 2, top = 10)


ggplot(pcClimate,aes(x=cv,y=PC1)) +
  geom_smooth(method='lm') +
  geom_point() +
  facet_wrap(variable ~ season,scales='free_x') +
  xlab("Coefficient of variation \n in climate variable from 2006-2018")

dev.off()

biplot(varcomp.pca)
text(varcomp.pca$x[,1], varcomp.pca$x[,2], labels=varianceComponents$site)

library("factoextra")


fviz_pca_var(varcomp.pca, col.var = "black",
             axes = c(1,4))

p1 <- fviz_pca_var(varcomp.pca, col.var = "black",
             axes = c(1,2))
  
p2<-ggplot() 
g1<-fviz_contrib(varcomp.pca, choice = "var", axes = 1, top = 10)
g2<-fviz_contrib(varcomp.pca, choice = "var", axes = 2, top = 10)
g3<-fviz_contrib(varcomp.pca, choice = "var", axes = 3, top = 10)
g4<-fviz_contrib(varcomp.pca, choice = "var", axes = 4, top = 10)
g5<-fviz_contrib(varcomp.pca, choice = "var", axes = 5, top = 10)
g6<-fviz_contrib(varcomp.pca, choice = "var", axes = 6, top = 10)
G<-gridExtra::grid.arrange(g1,g2,nrow=2)
gridExtra::grid.arrange(p1,G,nrow=1)


 
ggplot(pcMat,aes(x=PC1,y=PC3,label=site)) +
  geom_point() +
  geom_text_repel(size=3,color="black") +
  theme_bw() 

ggplot(pcMat,aes(x=PC1,y=PC4,label=site)) +
  geom_point() +
  geom_text_repel(size=3,color="black") +
  theme_bw() 

ggplot(pcMat,aes(x=PC1,y=PC5,label=site)) +
  geom_point() +
  geom_text_repel(size=3,color="black") +
  theme_bw() 

ggplot(pcMat,aes(x=PC1,y=PC6,label=site)) +
  geom_point() +
  geom_text_repel(size=3,color="black") +
  theme_bw() 

################################################################################
# Spatial data
#################################################################################
position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,northing,elevation) %>%
  dplyr::mutate(easting=easting/1000,northing=northing/1000)

pcMat <- pcMat %>%
  dplyr::left_join(position,by="site")

ggplot(pcMat,aes(x=easting,y=PC1,label=site)) +
  geom_point() +
  geom_text_repel(size=3,color="black") +
  theme_bw() 

ggplot(pcMat,aes(x=easting,y=PC2,label=site)) +
  geom_point() +
  geom_text_repel(size=3,color="black") +
  theme_bw() 

# Extra -------------------------------------------------------------------

data=read.csv(file = "~/Dropbox/Clarkia-LTREB/weather data/datafile_demography_site_environmental_data_2005-2019.csv",header=TRUE)

library(tidyverse)

data <- data %>%
  janitor::clean_names(case="lower_camel")

data = data %>%
  tidyr::pivot_longer(cols=contains(c("twinter","tspring","tsummer","pwinter","pspring","psummer")),
                      names_to = c("variable","year"),
                      names_pattern = "(.*)(.*)",
                      values_to = "estimate") %>% 
  tidyr::separate(variable,into = c("variable", "year"), "(?<=[a-z])(?=[0-9])") %>% 
  dplyr::mutate(year = substr(year,1,2)) %>% 
  dplyr::mutate(year = as.numeric(paste0(20,year))) %>%
  tidyr::separate(variable, c("variable","season"),1) %>%
  dplyr::select(-site) %>%
  dplyr::rename(site = site_2)

data = data %>%
  dplyr::filter(intenseDemography==1)

dataSummary<- data %>%
  dplyr::group_by(site,easting,season,variable) %>%
  dplyr::summarise(mu = mean(estimate,na.rm=TRUE),
                   sd = sd(estimate,na.rm=TRUE),
                   cv = sd/mu) %>%
  dplyr::mutate(variable = ifelse(variable=='p','precip','temp'))

dataSummary$season <- factor(dataSummary$season , levels=c("winter","spring","summer") )

pcClimate <-dataSummary %>%
  dplyr::left_join(pcMat %>% dplyr::select(-easting),by=c("site"))

ggplot(pcClimate,aes(x=mu,y=PC1)) +
  geom_smooth(method='lm') +
  geom_point() +
  facet_wrap(variable ~ season,scales='free_x')

ggplot(pcClimate,aes(x=sd,y=PC1)) +
  geom_smooth(method='lm') +
  geom_point() +
  facet_wrap(variable ~ season,scales='free_x')

ggplot(pcClimate,aes(x=cv,y=PC1)) +
  geom_smooth(method='lm') +
  geom_point() +
  facet_wrap(variable ~ season,scales='free_x')

varianceVsEastingTotal<-varianceComponents %>%
  dplyr::mutate(site=as.factor(site)) %>%
  tidyr::pivot_longer(cols=c(cov_sigma.fec,cov_sigma.phi,cov_fec.phi,var_sigma,var_fec,var_phi)) %>%
  dplyr::mutate(name=ifelse(name=="cov_sigma.fec","Cov(sigma,fec)",name),
                name=ifelse(name=="cov_sigma.phi","Cov(sigma,phi)",name),
                name=ifelse(name=="cov_fec.phi","Cov(fec,phi)",name),
                name=ifelse(name=="var_sigma","Var(sigma)",name),
                name=ifelse(name=="var_fec","Var(fec)",name),
                name=ifelse(name=="var_phi","Var(phi)",name)) %>%
  dplyr::arrange(desc(var.total),.by_group=TRUE) %>%
  dplyr::mutate(name = factor(name, 
                              levels = c("Var(sigma)", "Var(fec)", "Var(phi)",
                                         "Cov(sigma,fec)", "Cov(sigma,phi)",
                                         "Cov(fec,phi)"), ordered=TRUE)) %>%
  dplyr::left_join(position,by="site") %>%
  dplyr::mutate(perc=value/var.total) %>% 
  dplyr::select(site,var.total,easting) %>% unique %>%
  
  ggplot(aes(x=easting,y=var.total,label=site)) +
  geom_hline(yintercept=c(exp(0)),linetype='dotted',color='lightgray',size=.5) +
  geom_point() +
  #  facet_wrap(~name,scales='free',nrow=2) +
  theme_classic() +
  xlab("Easting") + ylab("Total variance in per-capita RS") +
  geom_text_repel(size=1.5,color="black") 

ggsave(filename=paste0(dirFigures,"varianceVsEastingTotal.pdf"),
       plot=varianceVsEastingTotal,width=10,height=6)

pdf(file=paste0(dirFigures,"climatevarianceVsFitnessvarianceTotal.pdf"),width=10,height=6)

varianceComponents %>%
  dplyr::left_join(dataSummary,by="site") %>%
  ggplot(aes(x=cv,y=var.total,label=site)) +
  geom_hline(yintercept=c(exp(0)),linetype='dotted',color='lightgray',size=.5) +
  geom_point() +
  facet_wrap(variable ~ season,scales='free_x') +
  theme_bw() +
  xlab("Coefficient of variation \n for climate variables") + ylab("Total variance in per-capita RS") +
  geom_text_repel(size=1.5,color="black")  +
  geom_smooth(method='lm')

varianceComponents %>%
  dplyr::left_join(dataSummary,by="site") %>%
  ggplot(aes(x=cv,y=var_sigma,label=site)) +
  geom_hline(yintercept=c(exp(0)),linetype='dotted',color='lightgray',size=.5) +
  geom_point() +
  facet_wrap(variable ~ season,scales='free_x') +
  theme_bw() +
  xlab("Coefficient of variation \n for climate variables") + ylab("Variance in seedling survival to fruiting") +
  geom_text_repel(size=1.5,color="black")  +
  geom_smooth(formula="y~x+x^2",method='lm')

varianceComponents %>%
  dplyr::left_join(dataSummary,by="site") %>%
  ggplot(aes(x=cv,y=var_fec,label=site)) +
  geom_hline(yintercept=c(exp(0)),linetype='dotted',color='lightgray',size=.5) +
  geom_point() +
  facet_wrap(variable ~ season,scales='free_x') +
  theme_bw() +
  xlab("Coefficient of variation \n for climate variables") + ylab("Variance in fruits per plant") +
  geom_text_repel(size=1.5,color="black")  +
  geom_smooth(method='lm')

varianceComponents %>%
  dplyr::left_join(dataSummary,by="site") %>%
  ggplot(aes(x=cv,y=var_phi,label=site)) +
  geom_hline(yintercept=c(exp(0)),linetype='dotted',color='lightgray',size=.5) +
  geom_point() +
  facet_wrap(variable ~ season,scales='free_x') +
  theme_bw() +
  xlab("Coefficient of variation \n for climate variables") + ylab("Variance in seeds per fruit") +
  geom_text_repel(size=1.5,color="black")  +
  geom_smooth(method='lm')


varianceComponents %>%
  dplyr::left_join(dataSummary,by="site") %>%
  ggplot(aes(x=cv,y=cov_sigma.fec,label=site)) +
  geom_hline(yintercept=c(exp(0)),linetype='dotted',color='lightgray',size=.5) +
  geom_point() +
  facet_wrap(variable ~ season,scales='free_x') +
  theme_bw() +
  xlab("Coefficient of variation \n for climate variables") + ylab("Covariance of survival and fruits per plant") +
  geom_text_repel(size=1.5,color="black")  +
  geom_smooth(method='lm')

varianceComponents %>%
  dplyr::left_join(dataSummary,by="site") %>%
  ggplot(aes(x=cv,y=cov_sigma.phi,label=site)) +
  geom_hline(yintercept=c(exp(0)),linetype='dotted',color='lightgray',size=.5) +
  geom_point() +
  facet_wrap(variable ~ season,scales='free_x') +
  theme_bw() +
  xlab("Coefficient of variation \n for climate variables") + ylab("Covariance of survival and seeds per fruit") +
  geom_text_repel(size=1.5,color="black")  +
  geom_smooth(method='lm')

varianceComponents %>%
  dplyr::left_join(dataSummary,by="site") %>%
  ggplot(aes(x=cv,y=cov_fec.phi,label=site)) +
  geom_hline(yintercept=c(exp(0)),linetype='dotted',color='lightgray',size=.5) +
  geom_point() +
  facet_wrap(variable ~ season,scales='free_x') +
  theme_bw() +
  xlab("Coefficient of variation \n for climate variables") + ylab("Covariance of fruits per plant and seeds per fruit") +
  geom_text_repel(size=1.5,color="black") +
  geom_smooth(method='lm')

dev.off()

# 
# 
# # https://quant.stackexchange.com/questions/17236/geometric-variance
# # Geometric variance is the interest rate per period over a n period time frame you need to compound to get some growth. 
# 
# ###
# rsMedians$log.rs = log(rsMedians$rs)
# 
# 
# 
# var.rs=rsMedians %>%
#   dplyr::group_by(site) %>%
#   dplyr::summarise(var_rs=var(rs,na.rm=TRUE))
#   dplyr::summarise(var.sigma = var(med_sigma,na.rm=TRUE),
#                    var.fec = var(med_fec,na.rm=TRUE),
#                    var.phi = var(med_phi,na.rm=TRUE))
# 
#   plot(as.factor(var.rs$site),var.rs$var_rs)
#   
# 
#   
# cov_mat <-  rsMedians %>%
#     tidyr::pivot_longer(cols = c("med_sigma", "med_fec", "med_phi")) %>%
#     group_by(site) %>%
#     summarise(cov_sigmafec = cov(value[name=="med_sigma"],value[name=="med_fec"],use="pairwise.complete.obs"),
#               cov_sigmaphi = cov(value[name=="med_sigma"],value[name=="med_phi"],use="pairwise.complete.obs"),
#               cov_fecphi = cov(value[name=="med_fec"],value[name=="med_phi"],use="pairwise.complete.obs"))  
#     
# plot(cov_mat$cov_sigmafec,cov_mat$cov_sigmaphi,xlim=c(-1,1),ylim=c(-1,1))
# plot(cov_mat$cov_sigmafec,cov_mat$cov_fecphi,xlim=c(-1,1),ylim=c(-25,25))
# plot(cov_mat$cov_fecphi,cov_mat$cov_sigmaphi,xlim=c(-25,25),ylim=c(-1,1))


A = matrix(c(1,0,4,
              0,2,3,
              3,0,),nrow=3)
n = c(1,1,1)

A %*% n

eigen(A)
eigen(t(A))

w = eigen(A)$vectors[,1] # right
v = eigen(t(A))$vectors[,1] # left

sens = matrix(NA,nrow=3,ncol=3)
for(i in 1:3){
  for(j in 1:3){
 sens[i,j] = w[i] %*% v[j]
 }
}
