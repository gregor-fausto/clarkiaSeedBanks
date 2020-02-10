
library(bayesplot)
color_scheme_set("brightblue")


rm(list=ls(all=TRUE)) # clear R environment

load(file="/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/output/seedbagfit.rds")
load(file="/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/output/seedbagdata.rds")

getParameters <- function(codaObject,x) { codaObject[,stringr::str_detect(colnames(codaObject),x)] }

## Seeds in seed bags in January
parameter = "alphaS1"

  merged.data.frame = Reduce(function(...) merge(..., all=T), zc)
  
  aS1 <- getParameters(merged.data.frame,parameter)
  ps<-apply(aS1,2,boot::inv.logit)


tmp <- matrix(nrow=3000,ncol=535)
for(i in 1:length(data$n)){
  tmp[,i]<-rbinom(n=length(ps[,data$site[i]]),size=data$n[i],p=ps[,data$site[i]])
}

samps<-sample(3000,1000)

bayesplot::ppc_dens_overlay(data$yt, tmp[samps,]) +
  theme_bw() + xlim(c(0,100)) + labs(title="Posterior predictive checks for seeds counted in seed bags in January", 
                                     caption="Dark line is the density of observed data (y) and the lighterlines show the densities of Y_rep from 1000 draws of the posterior")


bayesplot::ppc_stat_grouped(data$yt, tmp[samps,], group=interaction(data$site)) +
  theme_bw() + labs(title="Posterior predictive checks for the mean of seeds counted in seed bags in January", 
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")  

bayesplot::ppc_stat_grouped(data$yt, tmp[samps,],group=data$site, stat="var") +
  theme_bw() + labs(title="Posterior predictive checks for the variance of seeds counted in seed bags in January", 
                    caption="the bar is the observed value of test statistic T(y) and the histograms show T(Y_rep) from 1000 draws of the posterior")

# Seed survival (g1)
# yg.sim[i] ~ dbinom(pg[i]*(p[i])^(1/3), yt[i]) 


# Seed survival (s2)

# Seed survival (s3)

ps[i] = ilogit(alphaS1[site[i]])

# yt.sim[i] ~ dbinom(ps[i], n[i]) 
# yg.sim[i] ~ dbinom(pg[i]*(p[i])^(1/3), yt[i]) 
# yo.sim[i] ~ dbinom(pr[i], yt[i]-yg[i]) 


# yt2.sim[i] ~ dbinom(ps2[i]*(1-pg2[i])*pr2[i]*ps3[i], n2[i])

load(file="/Users/Gregor/Dropbox/modelsF2019/output/seedbagfit")
zc[[1]][,stringr::str_detect(colnames(zc[[1]]),pattern="yt.sim")]
