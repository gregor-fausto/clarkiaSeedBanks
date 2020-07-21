library(tidyverse)
library(ggrepel)

rsMedians <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/products/dataFiles/rsMedianEstimates.RDS")

# https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

A <- rsMedians[rsMedians$site=="BG",]$rs
n = length(A)
mu_g = (prod(A))^(1/n)

exp(sqrt(sum((log(A/mu_g))^2)/n))

# sample implementation
gsd.sample <- function(x){
  n = length(x)
  mu = (prod(x))^(1/n)
  y <- exp(sqrt(sum((log(x/mu))^2)/n))
  return(y)
}

gsd.am <- function(x){
  n = length(x)
  mu = exp(mean(log(x),na.rm=TRUE))
  y <- exp(sqrt(sum((log(x/mu))^2,na.rm=TRUE)/n))
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

# test gsd vs. sample gsd 
x <- c(rlnorm(n=10000,meanlog = 0 , sdlog = 1),NA)
gsd(x)
gsd.sample(x)
gsd.am(x)

x <- rnorm(n=1000,mean = 20 , sd = 1)
gsd(x)
gsd.am(x)


### Compare two methods of calculating gs

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
  ylab("Geometric SD of fitness (sample)") +
  xlim(c(2,7) ) + ylim(c(2,7)) 
  
###
rsMedians$log.rs = log(rsMedians$rs)



var.rs=rsMedians %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(var_rs=var(rs,na.rm=TRUE))
  dplyr::summarise(var.sigma = var(med_sigma,na.rm=TRUE),
                   var.fec = var(med_fec,na.rm=TRUE),
                   var.phi = var(med_phi,na.rm=TRUE))

  plot(as.factor(var.rs$site),var.rs$var_rs)
  

  
cov_mat <-  rsMedians %>%
    tidyr::pivot_longer(cols = c("med_sigma", "med_fec", "med_phi")) %>%
    group_by(site) %>%
    summarise(cov_sigmafec = cov(value[name=="med_sigma"],value[name=="med_fec"],use="pairwise.complete.obs"),
              cov_sigmaphi = cov(value[name=="med_sigma"],value[name=="med_phi"],use="pairwise.complete.obs"),
              cov_fecphi = cov(value[name=="med_fec"],value[name=="med_phi"],use="pairwise.complete.obs"))  
    
plot(cov_mat$cov_sigmafec,cov_mat$cov_sigmaphi,xlim=c(-1,1),ylim=c(-1,1))
plot(cov_mat$cov_sigmafec,cov_mat$cov_fecphi,xlim=c(-1,1),ylim=c(-25,25))
plot(cov_mat$cov_fecphi,cov_mat$cov_sigmaphi,xlim=c(-25,25),ylim=c(-1,1))
