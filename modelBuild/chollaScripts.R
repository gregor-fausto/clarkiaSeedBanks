library(tidyverse)
library(scales)
library(R2jags)
library(popbio)
invlogit<-function(x){exp(x)/(1+exp(x))}
volume <- function(h, w, p){
  (1/3)*pi*h*(((w + p)/2)/2)^2
}

cholla <- read.csv(url("https://raw.githubusercontent.com/texmiller/cholla_climate_IPM/master/cholla_demography_20042018.csv"))
cholla.clim <- cholla %>% 
  filter(Year_t >= min(cholla$Year_t,na.rm=T),
         Year_t <= max(cholla$Year_t,na.rm=T)) %>% 
  filter(str_sub(Plot,1,1)!="H") %>% 
  dplyr::mutate(year_int = Year_t - (min(Year_t)-1),
         plot_int = as.integer(as.factor(Plot)),
         vol_t = log(volume(h = Height_t, w = Width_t, p = Perp_t)),
         vol_t1 = log(volume(h = Height_t1, w = Width_t1, p = Perp_t1)),
         standvol_t = (vol_t - mean(vol_t,na.rm=T))/sd(vol_t,na.rm=T),
         standvol_t1 = (vol_t1 - mean(vol_t1,na.rm=T))/sd(vol_t1,na.rm=T))
## prep vital rate data sets

surv_dat <- cholla.clim %>% 
  filter(!is.na(standvol_t),
         !is.na(Survival_t1))

surv_dat<-surv_dat %>% dplyr::filter(plot_int>8)

cholla.dat<-list(N.plots = max(surv_dat$plot_int),
                 N.years = max(surv_dat$year_int),
                 
                 surv.N.obs = nrow(surv_dat),
                 surv.plot = surv_dat$plot_int,
                 surv.year = surv_dat$year_int,
                 surv.size = surv_dat$standvol_t,
                 surv.y = surv_dat$Survival_t1
)


inits<-function(){list(
                       surv.mu=rnorm(1),

                       surv.sigma.plot=rlnorm(1),

                       surv.sigma.year=rlnorm(1),
                       
                       surv.zsize=rbinom(1,1,0.5)
                       
)}

parameters<-c("surv.mu",
              "surv.bsize",
              "surv.zsize",
              "surv.sigma.plot",
              "surv.sigma.year",
              "surv.fit",
              "surv.fit.new")

## MCMC settings
ni<-40000
nb<-5000
nt<-10
nc<-3

allrates.out<-jags(data=cholla.dat,inits=inits,parameters.to.save=parameters,model.file="/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/jagsScripts/chollaJags.R",
                   n.thin=nt,n.chains=nc,n.burnin=nb,n.iter=ni,DIC=T)


# tuning (n.adapt)
n.adapt=3000
jm = jags.model("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/jagsScripts/chollaJags.R", data = cholla.dat, inits = inits,
                n.chains = 3, n.adapt = n.adapt)

# burn-in (n.update)
update(jm, n.iter = nb)
# chain (n.iter)
zc = coda.samples(jm, variable.names = c(parameters), n.iter = ni, thin = nt)

MCMCvis::MCMCsummary(zc,params=parameters)
MCMCvis::MCMCtrace(zc,params="surv.bsize")
plot(MCMCvis::MCMCchains(zc,params="surv.fit"),MCMCvis::MCMCchains(zc,params="surv.fit.new"))
abline(0,1, col='darkgray',lwd=3)

plot(allrates.out$BUGSoutput$sims.list$surv.fit,
     allrates.out$BUGSoutput$sims.list$surv.fit.new,
     xlab="SSQ for actual data",
     ylab="SSQ for perfect (new) data")
abline(0,1, col='darkgray',lwd=3)
