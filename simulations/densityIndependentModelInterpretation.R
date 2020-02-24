#################################################################################
################################################################################
################################################################################
# Code for *explaining* results from
# Optimizing Reproduction in a Randomly Varying Environment
# Dan Cohen
# Journal of Theoretical Biology 
# 1966
# Vol 12, p. 119-129
#
# Scripts by Gregor Siegmund
# fausto.siegmund@gmail.com
# last updated 02-24-2020
#################################################################################
#################################################################################
#################################################################################

# equation 4 from Cohen 1966
f<-function(G=g,D=d,Y=y,P=py){
  Fitness=P*log10((1-G)*(1-D)+G*Y)
  return(Fitness)
}

# script to calculate long term population growth rate
# for different combinations of parameters
v<-function(G=g,D=d,Y=y,P=py){
  vec <- c()
  for(i in 1:length(g)){
    vec[i]<-sum(f(G=g[i],D=D,Y=Y,P=P))
  }
  return(vec)
}

# script to calculate optimal g
max_g <- function(y,py,d=d[2]){
  out<-unlist(v(G=g,D=d,Y=y,P=py))
  g[which(out %in% max(out))];max(out)
  x<-seq(-7,7,by=.001)

  emp=1-ecdf(log(y,base=10))(x)
  gmax=g[which(out %in% max(out))]
  return(gmax)
}


## Figure 7 from Cohen
## updated with more even distribution of fitness
par(mfrow=c(3,1))
g=seq(0,1,by=.01)
d=c(0,.5,.8)

seqLen=100
fromVal=0.01; toVal=100
y=10^seq(log(fromVal,base=10), log(toVal,base=10), length.out = seqLen)
ly = log(y,base=10)
py = dnorm(ly,mean=0,sd=1)

plot(g,v(G=g,D=0,Y=y,P=py),type='n',ylim=c(-40,10),xlab="G",ylab="Long term expectation of population growth rate")
abline(h=0,col='gray')
for(i in 1:length(d)){
  lines(g,unlist(v(G=g,D=d[i],Y=y,P=py)),type='l',pch=16)
  out<-unlist(v(G=g,D=d[i],Y=y,P=py))
  points(g[which(out %in% max(out))],max(out),col='red',pch=16)
}
text(.5,7,"s=1")
text(.75,1.5,"s=.5")
text(.95,-1.5,"s=.2")

plot(ly,py,type='h',lwd=2,xlab=expression('Logarithm of fitness [log'[10]*'(Y)]'),ylab="Probability of fitness Y [P(Y)]")
plot(y,py,type='h',lwd=2,xlab="Fitness [Y]",ylab="Probability of fitness Y [P(Y)]")

## Effect of seed survivorship on optimal germination fraction
par(mfrow=c(2,1))
dvals = seq(0,1,by=.1)
g=seq(0,1,by=.01)
seqLen=100
fromVal=0.01; toVal=100
y=10^seq(log(fromVal,base=10), log(toVal,base=10), length.out = seqLen)
ly = log(y,base=10)

py = dnorm(ly,mean=0,sd=1)
plot(1-dvals,sapply(dvals,max_g,y=y,py=py),xlim=c(0,1),ylim=c(0,1),
     xlab="Seed survivorship [S]", ylab="Optimal germination fraction [G*]",type='l')
#text(.875,min(sapply(dvals,max_g,y=y,py=py))-.01,expression(paste(mu, "=0")))

py = dnorm(ly,mean=-.5,sd=1)
lines(1-dvals,sapply(dvals,max_g,y=y,py=py),lty='dashed')
#text(.875,min(sapply(dvals,max_g,y=y,py=py))-.01,expression(paste(mu, "=-0.5")))

py = dnorm(ly,mean=.5,sd=1)
lines(1-dvals,sapply(dvals,max_g,y=y,py=py),lty='dotted')
#text(.875,min(sapply(dvals,max_g,y=y,py=py))-.01,expression(paste(mu, "=0.5")))


# fitness distributions for the panel above
py = dnorm(ly,mean=0,sd=1)
plot(ly,py,type='n', xlab=expression('Logarithm of fitness [log'[10]*'(Y)]'),ylab= "Probability of fitness, Y [P(Y)]")
abline(v=0,col='gray')
lines(ly,py,type='l',lwd=2)

py = dnorm(ly,mean=-.5,sd=1)
lines(ly,py,type='l',lwd=2,lty='dashed')

py = dnorm(ly,mean=.5,sd=1)
lines(ly,py,type='l',lwd=2,lty='dotted')

## Effect of variance in fitness on optimal germination fraction
seqLen=100
fromVal=0.01; toVal=100
y=10^seq(log(fromVal,base=10), log(toVal,base=10), length.out = seqLen)
ly = log(y,base=10)
py=list()
vals=seq(.1,5,by=.25)
for(i in 1:length(vals)){
  py[[i]] = dnorm(ly,mean=0,sd=vals[i])
}


# plot variance in fitness fraction
# par(mfrow=c(1,2))
# plot(y,py[[i]],type='n',lwd=2,ylim=c(0,1))
# for(i in 1:length(vals)){
#   lines(y,py[[i]],type='l',lwd=1,col=i)
# }
# 
# plot(ly,py[[i]],type='n',lwd=2,ylim=c(0,1))
# for(i in 1:length(vals)){
#   lines(ly,py[[i]],type='l',lwd=1,col=i)
# }
#
# plot(g,v(G=g,D=0,Y=y,P=py[[1]]),type='n',ylim=c(-45,12))
# abline(h=0,col='gray')
# for(i in 1:length(vals)){
#   lines(g,unlist(v(G=g,D=d[2],Y=y,P=py[[i]])),type='l',pch=16)
#   out<-unlist(v(G=g,D=d[2],Y=y,P=py[[i]]))
#   points(g[which(out %in% max(out))],max(out),col='red',pch=16)
# }


par(mfrow=c(1,1))
plot(vals^2,lapply(py,max_g,y=y,d=.5),pch=16,ylim=c(0,1),col='gray',
     xlab="Fitness variance [var(Y)]",ylab="Optimal germination fraction [G*]")

points(vals^2,lapply(py,max_g,y=y,d=.2),pch=16,col='red')
points(vals^2,lapply(py,max_g,y=y,d=.8),pch=16)
text(18,.85,"Seed survivorship=0.2")
text(18,.66,"Seed survivorship=0.5")
text(18,.57,"Seed survivorship=0.8")

