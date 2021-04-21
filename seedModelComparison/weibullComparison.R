
exponential <- function(x,lambda){
  return(exp(-x*lambda))
}

weibull <- function(x,shape,scale){
  return(exp(-(x/scale)^shape))
}

c.exp <- function(x,a,b){
  return(1/(1+b*x)^a)
}

all.colors = c("#d73027","#fc8d59","#fee090","gray75","#e0f3f8","#91bfdb","#4575b4")

# -------------------------------------------------------------------
# Exponential; vary annual survival
# -------------------------------------------------------------------
t.max = 5
x=seq(0,t.max,by=.01)

par(mfcol=c(2,3),mar=c(2,2,2,2))
p.annual.vec = c(3/4,1/2,1/4)
colors=all.colors[c(3,4,5)]
offset = c(-.2,0,.2)

plot(NA,NA,type='n',ylim=c(0,1),xlim=c(0,t.max),xlab='',ylab='')
for(i in 1:3){
  p.annual = p.annual.vec[i]
  p.surv = ((p.annual)^(1:t.max))/((p.annual)^(0:(t.max-1)))
  y=exponential(x,lambda=-log(p.annual))
  lines(x,y,lwd=2,col=colors[i])
  #segments(x0=0,x1=1:t.max,y0=(p.annual)^(1:t.max),lty='dotted')
}

plot(NA,NA,type='n',xlim=c(0,t.max+.5),ylim=c(0,1),xlab='',ylab='')
for(i in 1:3){
  p.annual = p.annual.vec[i]
  p.surv = ((p.annual)^(1:t.max))/((p.annual)^(0:(t.max-1)))
  segments(x0=1:t.max+offset[i],y0=0,y1=p.surv,lty='dotted',col=colors[i])
  points(x=1:t.max+offset[i],y=p.surv,pch=19,col=colors[i])
}

# -------------------------------------------------------------------
# Weibull; p = 1/2
# vary shape and scale together
# -------------------------------------------------------------------
sol = function(x){
  return(1/(-log(.5))^(1/x))
}
wb.df = data.frame(scale = c(sol(.5),sol(1),sol(1.5)),shape=c(.5,1,1.5))
#par(mfcol=c(2,1),mar=c(2,2,2,2))

colors=all.colors[c(2,4,6)]
offset = c(-.2,0,.2)

tmp=wb.df

plot(NA,NA,type='n',ylim=c(0,1),xlim=c(0,t.max),xlab='',ylab='')
for(j in 1:3){
  y=weibull(x,shape=tmp$shape[j],scale=tmp$scale[j])
  lines(x,y,lwd=2,col=colors[j])
}

plot(NA,NA,type='n',xlim=c(0,t.max+.5),ylim=c(0,1),xlab='',ylab='')
for(j in 1:3){
  p.annual=weibull(1:t.max,shape=tmp$shape[j],scale=tmp$scale[j])
  p.surv = p.annual[1:t.max]/c(1,p.annual[1:(t.max-1)])
  segments(x0=1:t.max+offset[j],y0=0,y1=p.surv,lty='dotted',col=colors[j])
  points(x=1:t.max+offset[j],y=p.surv,pch=19,col=colors[j])
}

# -------------------------------------------------------------------
# Continuous exponential; p = 1/2
# vary shape and scale together
# -------------------------------------------------------------------
sol = function(x,p=1/2){
  return((1/p)^(1/x)-1)
}
ce.df = data.frame(a=c(.25,100,1),b = c(sol(.25),sol(100),sol(1)))
#par(mfcol=c(2,1),mar=c(2,2,2,2))

colors=all.colors[c(1,4,7)]
offset = c(-.2,0,.2)

tmp=ce.df

plot(NA,NA,type='n',ylim=c(0,1),xlim=c(0,t.max),xlab='',ylab='')
for(j in 1:3){
  y=c.exp(x,a=tmp$a[j],b=tmp$b[j])
  lines(x,y,lwd=2,col=colors[j])
}

plot(NA,NA,type='n',xlim=c(0,t.max+.5),ylim=c(0,1),xlab='',ylab='')
for(j in 1:3){
  p.annual=c.exp(1:t.max,a=tmp$a[j],b=tmp$b[j])
  p.surv = p.annual[1:t.max]/c(1,p.annual[1:(t.max-1)])
  segments(x0=1:t.max+offset[j],y0=0,y1=p.surv,lty='dotted',col=colors[j])
  points(x=1:t.max+offset[j],y=p.surv,pch=19,col=colors[j])
}


# -------------------------------------------------------------------
# Multipanel exponential
# -------------------------------------------------------------------
t.max = 5
x=seq(0,t.max,by=.01)

par(mfcol=c(2,3),mar=c(2,2,2,2))
scale.vec = c(3,2,1)
for(i in 1:3){
  
  plot(NA,NA,type='n',ylim=c(0,1),xlim=c(0,t.max),xlab='',ylab='')
  y=weibull(x,shape=1,scale=scale.vec[i])
  lines(x,y,col='black',lwd=2)
  p.annual=weibull(1:t.max,shape=1,scale=scale.vec[i])
  p.surv = p.annual[1:t.max]/c(1,p.annual[1:(t.max-1)])
  segments(x0=0,x1=1:t.max,y0=p.annual,lty='dotted')
  
  plot(NA,NA,type='n',xlim=c(0,t.max+.5),ylim=c(0,1),xlab='',ylab='')
  segments(x0=1:t.max,y0=0,y1=p.surv,lty='dotted')
  points(x=1:t.max,y=p.surv,pch=19)
}

t.max = 5
x=seq(0,t.max,by=.01)

sol = function(x){
  return(1/(-log(.5))^(1/x))
}
sol(1)

# -------------------------------------------------------------------
# Multipanel Weibull
# -------------------------------------------------------------------
wb.df=expand.grid(scale=c(-1/(log(1/2))),shape=c(1,1,1))
wb.df = data.frame(scale = c(sol(.5),sol(1),sol(1.5)),shape=c(.5,1,1.5))
par(mfcol=c(2,3),mar=c(2,2,2,2))
colors=c("#1b9e77","#d95f02","#7570b3")
offset = c(-.2,0,.2)
for(i in 1:3){
  
  tmp=wb.df[wb.df$scale==unique(wb.df$scale)[i],]
  
  plot(NA,NA,type='n',ylim=c(0,1),xlim=c(0,t.max),xlab='',ylab='')
  for(j in 1:3){
    y=weibull(x,shape=tmp$shape[j],scale=tmp$scale[j])
    lines(x,y,lwd=2,col=colors[j])
  }
  
  plot(NA,NA,type='n',xlim=c(0,t.max+.5),ylim=c(0,1),xlab='',ylab='')
  for(j in 1:3){
    p.annual=weibull(1:t.max,shape=tmp$shape[j],scale=tmp$scale[j])
    p.surv = p.annual[1:t.max]/c(1,p.annual[1:(t.max-1)])
    segments(x0=1:t.max+offset[j],y0=0,y1=p.surv,lty='dotted',col=colors[j])
    points(x=1:t.max+offset[j],y=p.surv,pch=19,col=colors[j])
  }
}







par(mfcol=c(2,3),mar=c(2,2,2,2))
shape.vec = c(.5,1,1.5)
for(i in 1:3){

  plot(NA,NA,type='n',ylim=c(0,1),xlim=c(0,t.max),xlab='',ylab='')
  y=weibull(x,shape=shape.vec[i],scale=1)
  lines(x,y,col='black',lwd=2)
  p.annual=weibull(1:t.max,shape=shape.vec[i],scale=1)
  p.surv = p.annual[1:t.max]/c(1,p.annual[1:(t.max-1)])
  segments(x0=0,x1=1:t.max,y0=p.annual,lty='dotted')
  
  plot(NA,NA,type='n',xlim=c(0,t.max),ylim=c(0,1),xlab='',ylab='')
  segments(x0=1:t.max,y0=0,y1=p.surv,lty='dotted')
  points(x=1:t.max,y=p.surv,pch=19)
}



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


A=matrix(c(0,0,0,
           1,0,0,
           1,0,0,
           1,1,0,
           1,1,0,
           1,1,1),nrow=6,byrow=TRUE)
g = c(.5,.5,.5)
((1-g[1])^A[,1])*((1-g[2])^(A[,2]))*((1-g[3])^(A[,3]))


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

