df<-read.csv(file="~/Downloads/Datafile_Survivorship_Fecundity.csv",header=TRUE)
dfSite <- read.csv(file="~/Downloads/Datafile_Site_Environment_corr.csv",header=TRUE)

library(janitor)
df<-df %>%
  janitor::clean_names(case="lower_camel") %>%
  dplyr::filter(year<2010)
dfSite<-dfSite %>%
  janitor::clean_names(case="lower_camel") %>%
  dplyr::select(site,easting) %>% unique %>%
  dplyr::filter(site %in% unique(df$site))

dim(df)

dfFiltered<-df %>%
  dplyr::filter(fruitingSeedling==FALSE)

## MLE

dfMLE<-dfFiltered %>%
  dplyr::mutate(prop = fruitingPlantNumber/seedlingNumber) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(prop_mu = mean(prop,na.rm=TRUE))

d1<-dfFiltered %>% dplyr::filter(site=="LCW") #%>%
  #dplyr::filter(!is.na(fruitingSeedling))
#d1[d1$fruitingSeedling==TRUE,]$seedlingNumber=3
negLogLik <- function(x,y,n) sum(-log(dbinom(y, n, x)))
opt<-optim(0.5, negLogLik, lower = 0, upper = 1, y = d1$fruitingPlantNumber, n = d1$seedlingNumber, method = "Brent")
p.seq <- seq(0, 0.99, 0.01)
plot(p.seq, sapply(p.seq, negLogLik, y = d1$fruitingPlantNumber, n = d1$seedlingNumber ), type="l")
points(opt$par,opt$value)
opt

d1.list <-split(d1, d1$year)

optima<-c()
for(i in 1:4){
  tmp<-d1.list[[i]]
  #opt<-optim(0.5, negLogLik, lower = 0, upper = 1, y = tmp$fruitingPlantNumber, n = tmp$seedlingNumber, method = "Brent")
  opt<-optimize(negLogLik,lower=0, upper=1,y = tmp$fruitingPlantNumber, n = tmp$seedlingNumber)
  optima[i] <- opt$minimum
}
optima %>% mean

library(ggplot2)

ggplot() +
  geom_point(data=dfFiltered,aes(x=site,y=fruitingPlantNumber/seedlingNumber),cex=.2) +
  geom_point(data=dfMLE,aes(x=site,y=prop_mu),color="green") +
  theme_bw()

## glm

model <- glm(cbind(fruitingPlantNumber,seedlingNumber) ~ 0 + site*as.factor(year) -site -as.factor(year) , data = dfFiltered, family="binomial")
summary(model)

x =unlist(stringr::str_split(names(coef(model)),"site"))
df2<-data.frame(site=x[x!=""],p=boot::inv.logit(coef(model)))

plot(dfMLE$prop_mu,df2$p);abline(a=0,b=1)

ggplot() +
  geom_point(data=dfFiltered,aes(x=site,y=fruitingPlantNumber/seedlingNumber),cex=.2) +
  geom_point(data=dfMLE,aes(x=site,y=prop_mu),color="green") +
  geom_point(data=df2,aes(x=site,y=p),color="red") +
  theme_bw()
  
## glmer

library(lme4)
dfFiltered$siteYear <- paste0(dfFiltered$year,dfFiltered$site)
model <- glmer(cbind(fruitingPlantNumber,seedlingNumber) ~ 0 + (1|siteYear), data = dfFiltered, family="binomial")

x =((ranef(model))[[1]])
df3<-data.frame(site=x[x!=""],p=boot::inv.logit(fixef(model)))

plot(dfMLE$prop_mu,boot::inv.logit(x$'(Intercept)'));abline(a=0,b=1)
plot(df2$p,boot::inv.logit(x$'(Intercept)'));abline(a=0,b=1)

View(data.frame(df2$p,boot::inv.logit(x$'(Intercept)')))

ggplot() +
  geom_point(data=dfFiltered,aes(x=site,y=fruitingPlantNumber/seedlingNumber),cex=.2) +
  geom_point(data=dfMLE,aes(x=site,y=prop_mu),color="green") +
  geom_point(data=df2,aes(x=site,y=p),color="red") +
  geom_point(data=df3,aes(x=site,y=p),color="blue") +
  theme_bw()

ggplot() +
  geom_point(data=dfFiltered %>% dplyr::left_join(dfSite,by="site"),aes(x=easting,y=fruitingPlantNumber/seedlingNumber),cex=.2) +
  geom_point(data=dfMLE %>% dplyr::left_join(dfSite,by="site"),aes(x=easting,y=prop_mu),color="green") +
  geom_point(data=df2 %>% dplyr::left_join(dfSite,by="site"),aes(x=easting,y=p),color="red") +
  geom_point(data=df3 %>% dplyr::left_join(dfSite,by="site"),aes(x=easting,y=p),color="blue") +
  geom_smooth(data=dfMLE %>% dplyr::left_join(dfSite,by="site"),aes(x=easting,y=prop_mu),color='green' ,method='lm') +
  geom_smooth(data=df2 %>% dplyr::left_join(dfSite,by="site"),aes(x=easting,y=p),color="red",method="lm") +
  geom_smooth(data=df3 %>% dplyr::left_join(dfSite,by="site"),aes(x=easting,y=p),color="blue",method='lm') +
  theme_bw()
