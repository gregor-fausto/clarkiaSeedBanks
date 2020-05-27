# compare estimates 
# this shows that estimates with small sample sizes are pooled towards the mean in Bayesian estimates relative to frequentist estimates

data= read.csv(file="~/Dropbox/clarkiaSeedBanks/products/dataFiles/seedlingSurvivalSummary-comparison.csv",header=TRUE)
head(data)

plot(data$X,data$med,xlim=c(0,1),ylim=c(0,1),
     xlab="Frequentist Estimates", ylab="Bayesian Estimates")
abline(a=0,b=1)

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedlingSurvival/"
simFiles <- paste0(directory,list.files(directory))
data1 = readRDS(simFiles[1])
data.sum = data1 %>%
  dplyr::mutate(year=as.numeric(year)) %>%
  dplyr::filter(seedlingNumber>0) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n=n())

data = data %>% dplyr::left_join(data.sum,by=c("site","year"))

plot(data$X,data$med,xlim=c(0,1),ylim=c(0,1),
     xlab="Frequentist Estimates", ylab="Bayesian Estimates",
     cex=data$n/30)
abline(a=0,b=1)


# compare F 
data= read.csv(file="~/Dropbox/clarkiaSeedBanks/products/dataFiles/fruitsPerPlantSummary-comparison.csv",header=TRUE)
head(data)

plot(data$X,data$med,xlim=c(0,25),ylim=c(0,25),
     xlab="Frequentist Estimates", ylab="Bayesian Estimates")
abline(a=0,b=1)

directory = "/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/fruitsPerPlantAllPlots/"
simFiles <- paste0(directory,list.files(directory))
data1 = readRDS(simFiles[1])
data.sum = data1 %>%
  dplyr::mutate(year=as.numeric(year)) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n=n())

data = data %>% dplyr::left_join(data.sum,by=c("site","year"))

par(mfrow=c(2,1))
plot(data$X,data$med,xlim=c(0,25),ylim=c(0,25),
     xlab="Frequentist Estimates", ylab="Bayesian Estimates",
     cex=1)
abline(a=0,b=1)

plot(data$X,data$med,xlim=c(0,25),ylim=c(0,25),
     xlab="Frequentist Estimates", ylab="Bayesian Estimates",
     cex=data$n/mean(data$n))
abline(a=0,b=1)
