## model checking for the viability model

load("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/output/viabilityModelFit.rds")

library(MCMCvis)

MCMCsummary(zc, params = c("mu.b","sigma.b"))

load("/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/output/viabilityModelData.rds")

library(ggplot2)

getParameters <- function(codaObject,x) { codaObject[,stringr::str_detect(colnames(codaObject),x)] }

## Seeds in seed bags in January
parameter = "p"

merged.data.frame = Reduce(function(...) merge(..., all=T), zc)

aS1 <- getParameters(merged.data.frame,parameter)
ps<-apply(aS1,2,boot::inv.logit)


tmp <- matrix(nrow=3000,ncol=544)
for(i in 1:length(data$nv)){
  tmp[,i]<-rbinom(n=length(ps[,data$site[i]]),size=data$nv[i],p=ps[,data$site[i]])
}

samps<-sample(3000,1000)

yv.sim =     # yv.sim[i] ~ dbinom(p[i], nv[i]) 


ggplot() +
  geom_histogram(data=data %>% as.data.frame,aes(yv)) +
  facet_wrap(~site) +
  theme_bw()

par(mfrow=c(4,5))
for(i in 1:data$nsites){
  d <- data %>% as.data.frame %>% dplyr::filter(site==i)
  hist(d$yv, breaks = 20, freq = FALSE) 
  lines(density(tmp[,data$site==i]), col = "red")
}