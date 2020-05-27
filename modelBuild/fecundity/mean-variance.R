# -------------------------------------------------------------------
# Mean-variance relationship for total fruit equivalents
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)

library(dplyr)

countFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countFruitsPerPlantAllPlots.RDS")

fruits<-countFruitsPerPlantAllPlots %>% 
  dplyr::filter(year %in% 2007:2012) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(mu=mean(countFruitNumberPerPlant,na.rm=TRUE),
                   v=var(countFruitNumberPerPlant,na.rm=TRUE),
                   n=n())

## mean variance relationship
plot(fruits$mu,fruits$v,pch=16,cex=.5,col=fruits$year,
     xlab="Mean (site*year)", ylab="Variance (site*year)",
     main="Fruits per plant",
     xlim=c(0,16))

rhs <- function(x, b0, b1) {
  b0 + x^b1
}

m.1 <-nls(v ~ rhs(mu, intercept, power), 
        data = fruits, 
        start = list(intercept = 0, power = 1), trace = T)

summary(m.1)

s <- seq(0, 16, length = 100)

plot(fruits$v ~ fruits$mu, pch=1 , 
     xlab="Mean (site*year)", ylab="Variance (site*year)",
     main="Fruits per plant: fitted power model, with intercept",
     xlim=c(0,16))
abline(h = 0, lty = 1, lwd = 0.5)
lines(s, predict(m.1, list(mu = s)), lty = 1, col = 'blue')

# -------------------------------------------------------------------
# Mean-variance relationship for undamaged fruits per plant
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)

library(dplyr)

countFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countUndamagedDamagedFruitsPerPlantAllPlots.RDS")

fruits<-countFruitsPerPlantAllPlots %>% 
  dplyr::filter(year %in% 2013:2018) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(mu=mean(countUndamagedFruitNumberPerPlant,na.rm=TRUE),
                   v=var(countUndamagedFruitNumberPerPlant,na.rm=TRUE),
                   n=n())

## mean variance relationship
plot(fruits$mu,fruits$v,pch=16,cex=.5,col=fruits$year,
     xlab="Mean (site*year)", ylab="Variance (site*year)",
     main="Undamaged fruits per plant",
     xlim=c(0,40))

rhs <- function(x, b0, b1) {
  b0 + x^b1
}

m.1 <-nls(v ~ rhs(mu, intercept, power), 
          data = fruits, 
          start = list(intercept = 0, power = 1), trace = T)

summary(m.1)

s <- seq(0, 40, length = 100)

plot(fruits$v ~ fruits$mu, pch=1 , 
     xlab="Mean (site*year)", ylab="Variance (site*year)",
     main="Fruits per plant: fitted power model, with intercept",
     xlim=c(0,40))
abline(h = 0, lty = 1, lwd = 0.5)
lines(s, predict(m.1, list(mu = s)), lty = 1, col = 'blue')



## seeds
load(file="phiIndDF.RData")

seeds<-phiIndDF %>% 
  dplyr::filter(year %in% 2007:2011) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(mu=mean(response,na.rm=TRUE),
                   v=var(response,na.rm=TRUE),
                   n=n())

plot(seeds$mu,seeds$v,pch=16,cex=.5,col=seeds$year,
     xlab="Mean (site*year)", ylab="Variance (site*year)",
     main="Seeds per fruit",
     xlim=c(0,60))

m.2 <-nls(v ~ rhs(mu, intercept, power), 
          data = seeds, 
          start = list(intercept = 0, power = 1), trace = T)

summary(m.2)

s <- seq(0, 60, length = 100)

plot(seeds$v ~ seeds$mu, pch=1 , 
     xlab="Mean (site*year)", ylab="Variance (site*year)",
     main="Seeds per fruit: fitted power model, with intercept",
     xlim=c(0,60))
abline(h = 0, lty = 1, lwd = 0.5)
lines(s, predict(m.2, list(mu = s)), lty = 1, col = 'blue')
