weibull <- function(x,shape,scale){
  return(exp(-(x/scale)^shape))
}

A=matrix(c(0,0,0,
         1,0,0,
         1,0,0,
         1,1,0,
         1,1,0,
         1,1,1),nrow=6,byrow=TRUE)
g = c(.5,.5,.5)
((1-g[1])^A[,1])*((1-g[2])^(A[,2]))*((1-g[3])^(A[,3]))

exponential <- function(x,lambda){
  return(exp(-x*lambda))
}

x=seq(0,10,by=.01)
plot(x,y,type='n',ylim=c(0,1))
y=exponential(x,lambda=1)
lines(x,y,col='red')
y=exponential(x,lambda=.5)
lines(x,y,col='red')
y=exponential(x,lambda=10)
lines(x,y,col='red')

points(-log(.5),.5,pch=16)

y=weibull(x,shape=1,scale=1)
lines(x,y,col='black',lty='dotted')
y=weibull(x,shape=1,scale=.1)
lines(x,y,col='black',lty='dotted')
y=weibull(x,shape=1,scale=2)
lines(x,y,col='black',lty='dotted')

points((-log(.5))^(1/.1),.5,pch=16)
points((-log(.5))^(1/10),.5,pch=16)

#abline(v=1)

## PLOT 2

par(mfrow=c(1,2))
x=seq(0,10,by=.01)
plot(x,y,type='n',ylim=c(0,1))
y=exponential(x,lambda=1)
lines(x,y,col='red')

y=weibull(x,shape=1,scale=1)
lines(x,y,col='black',lty='dotted')

y=weibull(x,shape=1,scale=2)
lines(x,y,col='black',lty='dotted')

y=weibull(x,shape=1,scale=.5)
lines(x,y,col='black',lty='dotted')

plot(x,y,type='n',ylim=c(0,1))
y=exponential(x,lambda=1)
lines(x,y,col='red')

y=weibull(x,shape=.5,scale=1)
lines(x,y,col='black',lty='dashed')

y=weibull(x,shape=.5,scale=2)
lines(x,y,col='black',lty='dashed')

y=weibull(x,shape=.5,scale=.5)
lines(x,y,col='black',lty='dashed')


### GERMINATION
par(mfrow=c(1,3))
A=matrix(c(0,0,0,
           1,0,0,
           1,0,0,
           1,1,0,
           1,1,0,
           1,1,1),nrow=6,byrow=TRUE)

g = c(.2,.2,.2)

weibull <- function(x,shape,scale){
  return(exp(-(x/scale)^shape))
}

t<-c(3,12,15,24,27,36)
y<-((1-g[1])^A[,1])*((1-g[2])^(A[,2]))*((1-g[3])^(A[,3]))

plot(t,y,pch=16,ylim=c(0,1))

y<-((1-g[1])^A[,1])*((1-g[2])^(A[,2]))*((1-g[3])^(A[,3]))*weibull(t,1,36)
points(t,y,pch=1,ylim=c(0,1))
#lines(t,weibull(t,1,36))

par(mfrow=c(1,1))
g = c(.1,.1,.1)
y<-((1-g[1])^A[,1])*((1-g[2])^(A[,2]))*((1-g[3])^(A[,3]))*weibull(t,.25,36)

plot(t,y,pch=16,ylim=c(0,1))

g = c(.1,.2,.3)
y<-((1-g[1])^A[,1])*((1-g[2])^(A[,2]))*((1-g[3])^(A[,3]))*weibull(t,.25,36)
points(t,y,pch=16,col='red')

g = c(.1,.05,.025)
y<-((1-g[1])^A[,1])*((1-g[2])^(A[,2]))*((1-g[3])^(A[,3]))*weibull(t,.25,36)
points(t,y,pch=16,col='green')
