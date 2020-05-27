## FIGURE 1

years = 100
n = 50
mu = .5
sigma = 0.5

p = rnorm(n = years, mean = mu, sd = sigma)
u = boot::inv.logit(p)

S = c()
for(i in 1:years){
S[i] = rbinom(n=1,size=n,prob=u[i])
}


plot(S/n,u,
     xlim=c(0,1),
     ylim=c(0,1),
     pch=16,col='red',
     xlab="Estimate p(t)",
     ylab="True p(t)")

p.hat = (S/n)
lm1<-lm(u ~ p.hat)
abline(a=coef(lm1)[1],b=coef(lm1)[2],col='red')

abline(a=0,b=1)

se.p.hat = sqrt((p.hat*(1-p.hat))/n)

# gould
mix.var = var(p.hat)-mean(se.p.hat^2)
gould.p.hat = mean(p.hat) + sqrt(mix.var/(mix.var+se.p.hat^2))*(p.hat-mean(p.hat))

points(gould.p.hat,u,col='blue')

lm1<-lm(u ~ gould.p.hat)
abline(a=coef(lm1)[1],b=coef(lm1)[2],col='blue')

# burnham
w = 1/(sigma^2+se.p.hat^2)
mean.p.hat = sum(w*p.hat)/sum(w)
burnham.p.hat = mean.p.hat + sqrt(mix.var/(mix.var+se.p.hat^2))*(p.hat-mean.p.hat)

points(burnham.p.hat,u,col='purple',cex=0.5)

lm1<-lm(u ~ burnham.p.hat)
abline(a=coef(lm1)[1],b=coef(lm1)[2],col='purple')

# shrinkage
gould.p.hat.shrinkage = mean(p.hat) + mix.var/(mix.var+se.p.hat^2)*(p.hat-mean(p.hat))

points(gould.p.hat.shrinkage,u,col='green')

lm1<-lm(u ~ gould.p.hat.shrinkage)
abline(a=coef(lm1)[1],b=coef(lm1)[2],col='green')

## FIGURE 2
## mean var

experiment <- function(years=25, N = 100){

mu = .5
sigma = 0.5

p = rnorm(n = years, mean = mu, sd = sqrt(sigma))
u = boot::inv.logit(p)

return(list(u,p))
}

sample.var<-apply(replicate(n=1000,experiment(years = 25, N = 100)[[2]]),2,var)
plot(sample.var)

u = experiment(years = 25, N = 100)[[1]]
p = replicate(n=100,experiment(years = 25, N = 100)[[2]])


S = c()
for(i in 1:years){
  S[i] = rbinom(n=1,size=n,prob=u[i])
}

var(u)
se.u = 1/(sqrt(n*p.hat*(1-p.hat)))

### data

data<-readRDS("~/Dropbox/dataLibrary/postProcessingData/censusSeedlingsFruitingPlants.RDS")

library(dplyr)

temp <- data %>%
  dplyr::filter(site=="BG")

temp = temp[temp$fruitplNumber<=temp$seedlingNumber,]

plot(as.factor(temp$year),temp$fruitplNumber/temp$seedlingNumber)


glm1 <- glm(cbind(fruitplNumber,seedlingNumber)~0 + as.factor(year),data=temp,family="binomial")

summary(glm1)
cis<-apply(confint(glm1),2,boot::inv.logit)



plot(temp$year,temp$fruitplNumber/temp$seedlingNumber)
abline(h=mean(temp$fruitplNumber/temp$seedlingNumber,na.rm=TRUE),lty='dotted')

segments(x0=2006:2017,y0=cis[,1],y1=cis[,2],col="red",lwd=1.5)
points(2006:2017,boot::inv.logit(coef(glm1)),pch=16,col='red')

year.averages = temp %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(avg = mean(fruitplNumber/seedlingNumber,na.rm=TRUE))

points(2006:2017,year.averages$avg[1:12],pch=16,col='blue')


### glmer

library(lme4)

glm1 <- glm(cbind(fruitplNumber,seedlingNumber)~0 + as.factor(year),data=temp,family="binomial")
glmer1 <- glmer(cbind(fruitplNumber,seedlingNumber)~1 + (1|year),data=temp,family="binomial")

summary(glmer1)

plot(temp$year,temp$fruitplNumber/temp$seedlingNumber)

abline(h=mean(temp$fruitplNumber/temp$seedlingNumber,na.rm=TRUE),lty='dotted')
abline(h=boot::inv.logit(fixef(glmer1)),lty='dotted',col='red')

points(2006:2017,boot::inv.logit(fixef(glmer1) + ranef(glmer1)$year[,1]),pch=16,col='red')
points(2006:2017,boot::inv.logit(coef(glm1)),pch=16,col='orange')


year.averages = temp %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(avg = mean(fruitplNumber/seedlingNumber,na.rm=TRUE))

points(2006:2017,year.averages$avg[1:12],pch=16,col='blue')

temp %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(n())
