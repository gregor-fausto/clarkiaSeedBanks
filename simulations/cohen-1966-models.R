#################################################################################
################################################################################
################################################################################
# Code for reproducing figures from 
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

################################################################################
################################################################################
# Figure 1
################################################################################
################################################################################

par(mfrow=c(1,1))

# Equation 5 
# long term population growth rate
# when Yt = c(0, Y)
# note use of log10
e<-function(G=g,D=d,Y=y,P=py){
  E=(1-P)*log10((1-G)*(1-D)+G*Y[1])+(P)*log10((1-G)*(1-D)+G*Y[2])
  return(E)
}

# script to compute long-term population growth rate
# for combinations of P_y, Y, D, and G
# varying G
v<-function(G=g,D=d,Y=y,P=py){
  l<-list()
  for(i in 1:length(g)){
    l[[i]]<-e(G=g[i],D=D,Y=Y,P=P)
  }
  return(l)
}

# data frame of parameter combinations to create Figure 1
df<-data.frame(py=c(.8,.8,.1,.1,.8,.8,.1,.1),
               y=c(500,500,500,500,5,5,5,5),
               d=c(.1,.8,.1,.8,.1,.8,.1,.8))

# create sequence of germination fractions
g=seq(0.01,1,by=.01)

setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
pdf(file="cohenFigure1.pdf", width=8, height=8)
par(mar=c(5, 6, 4, 2) + 0.1)
plot(g,unlist(v(G=g,D=df$d[1],Y=c(0,df$y[1]),P=df$py[1])),
     type='n',ylim=c(-2.5,2.5),
     xlab="Germination fraction (g)",
     ylab="Long-term expectation of growth rate \n f(g) for varying P, Y, and D")
abline(h=0,col='gray')
for(i in 1:(dim(df)[1])){
  lines(g,unlist(v(G=g,D=df$d[i],Y=c(0,df$y[i]),P=df$py[i])),type='l',pch=16)
}
par(op)
dev.off()

################################################################################
################################################################################
# Figure 2
# one of the parameter combinations appears mislabeled in Cohen
################################################################################
################################################################################

# Equation 5 
# long term population growth rate
# when Yt = c(0, Y)
# note use of log10
e<-function(G=g,D=d,Y=y,P=py){
  E=(1-P)*log10((1-G)*(1-D)+G*Y[1])+(P)*log10((1-G)*(1-D)+G*Y[2])
  return(E)
}

# script to compute long-term population growth rate
# for combinations of P_y, Y, D, and G
# varying D
v<-function(G=g,D=d,Y=y,P=py){
  l<-list()
  for(i in 1:length(d)){
    l[[i]]<-e(G=G,D=d[i],Y=Y,P=P)
  }
  return(l)
}

# data frame of parameter combinations to create Figure 2
df<-data.frame(py=c(.8,.8,.1,.1,.8,.8,.1,.1),
               y=c(500,500,500,500,5,5,5,5),
               g=c(.1,.8,.1,.8,.1,.8,.1,.8))

# create sequence of decay fractions
d=seq(0,1,by=.05)

setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
pdf(file="cohenFigure2.pdf", width=8, height=8)
par(mar=c(5, 6, 4, 2) + 0.1)
plot(-log(1-d),unlist(v(G=df$g[1],D=d,Y=c(0,df$y[1]),P=df$py[1])),
     type='n',ylim=c(-2,2),
     xlab="-ln(1-D)",
     ylab="Long-term expectation of growth rate \n f[-ln(1-d)] for varying P, Y, and G")
abline(h=0,col='gray')
for(i in 1:(dim(df)[1])){
  lines(-log(1-d),unlist(v(G=df$g[i],D=d,Y=c(0,df$y[i]),P=df$py[i])),type='l',pch=16)
}
dev.off()

################################################################################
################################################################################
# Figure 3
# cant quite recreate figure 3 
################################################################################
################################################################################

# Equation 5 
# long term population growth rate
# when Yt = c(0, Y)
# note use of log10
e<-function(G=g,D=d,Y=y,P=py){
  E=(1-P)*log10((1-G)*(1-D))+(P)*log10((1-G)*(1-D)+G*Y[2])
  return(E)
}


# script to compute long-term population growth rate
# for combinations of P_y, Y, D, and G
# varying ln(Y)
v<-function(G=g,D=d,Y=y,P=py){
  l<-list()
  for(i in 1:length(Y)){
    l[[i]]<-e(G=G,D=D,Y=c(0,Y[i]),P=P)
  }
  return(l)
}

# data frame to create Figure 3
df<-data.frame(g=c(.5,.1,.5,.1,.5,.95,.95,.1,.95),
               d=c(.1,.1,.3,.3,.8,.1,.3,.8,.8))

# create combination of parameters 
y = c(2:10 %o% 10^(0:3))
py = .5

# setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
# pdf(file="cohenFigure3.pdf", width=8, height=8)
par(mar=c(5, 6, 4, 2) + 0.1)
plot(log(y),unlist(v(G=df$g[1],D=df$d[1],Y=y,P=py)),type='n',ylim=c(-0.8,1.7),
     xlab="ln(Y)",
     ylab="Long-term expectation of growth rate \n f[ln(Y)] for P = 0.5")
abline(h=0,col='gray')
for(i in 1:(dim(df)[1])){
  lines(log(y),unlist(v(G=df$g[i],D=df$d[i],Y=y,P=py)),type='l',pch=16)
}
# dev.off()

################################################################################
################################################################################
# Figure 4
################################################################################
################################################################################

# Equation 5 
# long term population growth rate
# when Yt = c(0, Y)
# note use of log10
e<-function(G=g,D=d,Y=y,P=py){
  E=(1-P)*log10((1-G)*(1-D))+(P)*log10((1-G)*(1-D)+G*Y[2])
  return(E)
}

# script to compute long-term population growth rate
# for combinations of Py, Y, D, and G
# varying Py
v<-function(G=g,D=d,Y=y,P=py){
  l<-list()
  for(i in 1:length(P)){
    l[[i]]<-e(G=G,D=D,Y=c(0,Y),P=py[i])
  }
  return(l)
}


# data frame for figure 4
df<-data.frame(y=c(100,5,100,5,100,5,100,5),
               d=c(.1,.1,.8,.8,.1,.1,.8,.8),
               g=c(.3,.3,.3,.3,.95,.95,.95,.95))

# generate grid of py 
py = seq(0,1,by=.1)

setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
pdf(file="cohenFigure4.pdf", width=8, height=8)
par(mar=c(5, 6, 4, 2) + 0.1)
plot(py,unlist(v(G=df$g[1],D=df$d[1],Y=df$y[1],P=py)),
     type='n',ylim=c(-2,2), xlab="P(Y)",
     ylab="Long-term expectation of growth rate \n f[P(Y)] for varying Y, D, and G")
abline(h=0,col='gray')
for(i in 1:(dim(df)[1])){
  lines(py,unlist(v(G=df$g[i],D=df$d[i],Y=df$y[i],P=py)),type='l',pch=16)
}
dev.off()

################################################################################
################################################################################
# Figure 5
################################################################################
################################################################################

# equation 8 in Cohen
gmax<-function(D=d,Y=y,P=py){
  g = P - (1-P)*((1-D)/(Y+D-1)) 
  return(g)
}

df<-data.frame(y=c(500,10,2,2,2),
               d=c(.8,.3,.8,.3,.1))

py = seq(0,1,by=.1)

setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
pdf(file="cohenFigure5.pdf", width=8, height=8)
par(mar=c(5, 6, 4, 2) + 0.1)
plot(py,unlist(gmax(D=df$d[1],Y=df$y[1],P=py)),type='n',ylim=c(-1,1),
     xlab='P(good year)',ylab="G*")
abline(h=0,col='gray')
for(i in 1:(dim(df)[1])){
  lines(py,unlist(gmax(D=df$d[i],Y=df$y[i],P=py)),type='l',pch=16)
}
dev.off()


################################################################################
################################################################################
# Figure 7
################################################################################
################################################################################

# equation 4
f<-function(G=g,D=d,Y=y,P=py){
  Fitness=P*log10((1-G)*(1-D)+G*Y)
  return(Fitness)
}

# function development
# g=.5;d=0;y=c(0,10);py=c(0.2,.8)
# f(G=g,D=d,Y=y,P=py)
# c(f(G=g,D=d,Y=y[1],P=py[1]),f(G=g,D=d,Y=y[2],P=py[2]))
# 
# g=0;d=0;y=c(.1,.2,.5,1,2,5,10);py=dnorm(log(y,base=10),mean=0,sd=1)
# sum(f(G=g,D=d,Y=y,P=py))

v<-function(G=g,D=d,Y=y,P=py){
  vec <- c()
  for(i in 1:length(g)){
    vec[i]<-sum(f(G=g[i],D=D,Y=Y,P=P))
  }
  return(vec)
}


# sequence of g and d for plots
g=seq(0,1,by=.01)
d=c(0,.5,.8)

# 7A
# y = c(.1,.2,.5,1,2,5,10)
# ly = log(y,base=10)
# py = c(0.05,.1,.2,.3,.2,.1,.05)
# plot(dnorm(ly,mean=0,sd=1),type='h',lwd=2)
# 
# par(mfrow=c(1,1))
# plot(g,v(G=g,D=0,Y=y,P=py),type='n',ylim=c(-1,.25))
# abline(h=0,col='gray')
# for(i in 1:length(d)){
#   lines(g,unlist(v(G=g,D=d[i],Y=y,P=py)),type='l',pch=16)
#   out<-unlist(v(G=g,D=d[i],Y=y,P=py))
#   points(g[which(out %in% max(out))],max(out),col='red',pch=16)
# }

# Panel A
setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
pdf(file="cohenFigure7A.pdf", width=4, height=8)
par(mfrow=c(3,1))
y = c(.1,.2,.5,1,2,5,10)
ly = log(y,base=10)
py = c(0.05,.1,.2,.3,.2,.1,.05)
#plot(dnorm(ly,mean=0,sd=1),type='h',lwd=2)

# otpimal germination fraction plot
plot(g,v(G=g,D=0,Y=y,P=py),type='n',ylim=c(-1,.25),
     xlab="Germination fraction (G)",
     ylab="Long-term expectation of growth rate f[G]")
abline(h=0,col='gray')
for(i in 1:length(d)){
  lines(g,unlist(v(G=g,D=d[i],Y=y,P=py)),type='l',pch=16)
  out<-unlist(v(G=g,D=d[i],Y=y,P=py))
  points(g[which(out %in% max(out))],max(out),col='red',pch=16)
}

# log base 10 distribution of fitness
plot(log(y,base=10),py,type="h",
     lwd = 4,ylim=c(0,1),
     xlab="Log10(Y)",
     ylab="P(Y)")

# untransformed distribution of fitness
plot(y,py,type="h", lwd = 4,
xlab="Y",
ylab="P(Y)")
dev.off()

# Panel B
setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
pdf(file="cohenFigure7B.pdf", width=4, height=8)
par(mfrow=c(3,1))
y=c(.01,.02,.2,2,20,200)
ly = log(y,base=10)
py = c(0.15,.2,.3,.2,.1,.05)
#py=dnorm(log(y,base=10),mean=log(1),sd=log(2.25))

plot(g,v(G=g,D=0,Y=y,P=py),type='n',ylim=c(-1,.25),
     xlab="Germination fraction (G)",
     ylab="Long-term expectation of growth rate f[G]")
abline(h=0,col='gray')
for(i in 1:length(d)){
  lines(g,unlist(v(G=g,D=d[i],Y=y,P=py)),type='l',pch=16)
  out<-unlist(v(G=g,D=d[i],Y=y,P=py))
  points(g[which(out %in% max(out))],max(out),col='red',pch=16)
}

plot(log(y,base=10),py,type="h", lwd = 4,ylim=c(0,1),
     xlab="Log10(Y)",
     ylab="P(Y)")

plot(y,py,type="h", lwd = 4,
     xlab="Y",
     ylab="P(Y)")

dev.off()

# Panel C
setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
pdf(file="cohenFigure7C.pdf", width=4, height=8)
par(mfrow=c(3,1))
y=c(.01,.1,1,10,100,1000)
ly = log(y,base=10)
py = c(0.15,.2,.3,.2,.1,.05)
#py=dnorm(log(y,base=10),mean=log(1),sd=log(2.25))

plot(g,v(G=g,D=0,Y=y,P=py),type='n',ylim=c(-1,.35),
     xlab="Germination fraction (G)",
     ylab="Long-term expectation of growth rate f[G]")
abline(h=0,col='gray')
for(i in 1:length(d)){
  lines(g,unlist(v(G=g,D=d[i],Y=y,P=py)),type='l',pch=16)
  out<-unlist(v(G=g,D=d[i],Y=y,P=py))
  points(g[which(out %in% max(out))],max(out),col='red',pch=16)
}

plot(log(y,base=10),py,type="h", lwd = 4,ylim=c(0,1),
     xlab="Log10(Y)",
     ylab="P(Y)")

plot(y,py,type="h", lwd = 4,
     xlab="Y",
     ylab="P(Y)")

dev.off()

# Panel D
setwd("~/Dropbox/clarkiaSeedBanks/products/figures")
pdf(file="cohenFigure7D.pdf", width=4, height=8)
par(mfrow=c(3,1))
y=c(.01,.2,10,500,1000)
ly = log(y,base=10)
py = c(0.15,.2,.3,.2,.15)
#py=dnorm(log(y,base=10),mean=log(1),sd=log(2.25))

plot(g,v(G=g,D=0,Y=y,P=py),type='n',ylim=c(-1,1.1),
     xlab="Germination fraction (G)",
     ylab="Long-term expectation of growth rate f[G]")
abline(h=0,col='gray')
for(i in 1:length(d)){
  lines(g,unlist(v(G=g,D=d[i],Y=y,P=py)),type='l',pch=16)
  out<-unlist(v(G=g,D=d[i],Y=y,P=py))
  points(g[which(out %in% max(out))],max(out),col='red',pch=16)
}

plot(log(y,base=10),py,type="h", lwd = 4,ylim=c(0,1),
     xlab="Log10(Y)",
     ylab="P(Y)")

plot(y,py,type="h", lwd = 4,
     xlab="Y",
     ylab="P(Y)")

dev.off()


# par(mfcol=c(2,3))
# 
# seqLen = 50
# y = c(.1,.2,.5,1,2,5,10)
# ly = log(y,base=10)
# py = c(0.05,.1,.2,.3,.2,.1,.05)
# plot(dnorm(ly,mean=0,sd=1),type='h',lwd=2)
# 
# par(mfrow=c(1,1))
# plot(g,v(G=g,D=0,Y=y,P=py),type='n',ylim=c(-2,1))
# abline(h=0,col='gray')
# for(i in 1:length(d)){
#   lines(g,unlist(v(G=g,D=d[i],Y=y,P=py)),type='b',pch=16)
#   out<-unlist(v(G=g,D=d[i],Y=y,P=py))
#   points(g[which(out %in% max(out))],max(out),col='red',pch=16)
# }
# 
# max_g <- function(y,py){
#   out<-unlist(v(G=g,D=d[2],Y=y,P=py))
#   g[which(out %in% max(out))];max(out)
#   x<-seq(-7,7,by=.001)
#   #plot(x,1 - ecdf(log(y,base=10))(x),ylim=c(0,1),type='l')
#   
#   emp=1-ecdf(log(y,base=10))(x)
#   gmax=g[which(out %in% max(out))]
#   return(gmax)
# }
# 
# 
# 
# vals<-exp(seq(log(1,base), log(1000), length.out = seqLen))
# 
# y = list()
# py = list()
# for(i in 1:seqLen){
#   y[[i]] = rnorm(n=100,mean=log(1),sd=log(vals[i]))
#   y[[i]]=10^y[[i]]
#   py[[i]]=dnorm(log(y[[i]],base=10),mean=log(1),sd=log(vals[i]))
# }
# 
# par(mfrow=c(1,1))
# plot(g,v(G=g,D=0,Y=y[[seqLen]],P=py[[seqLen]]),type='n')
# abline(h=0,col='gray')
# for(i in 1:length(py)){
#   lines(g,unlist(v(G=g,D=d[3],Y=y[[i]],P=py[[i]])),type='l',pch=16)
#   out<-unlist(v(G=g,D=d[3],Y=y[[i]],P=py[[i]]))
#   points(g[which(out %in% max(out))],max(out),col='red',pch=16)
# }
# 
# max_g <- function(y,py){
#   out<-unlist(v(G=g,D=d[2],Y=y,P=py))
#   g[which(out %in% max(out))];max(out)
#   x<-seq(-7,7,by=.001)
#   #plot(x,1 - ecdf(log(y,base=10))(x),ylim=c(0,1),type='l')
#   
#   emp=1-ecdf(log(y,base=10))(x)
#   gmax=g[which(out %in% max(out))]
#   return(gmax)
# }
# 
# g_max = mapply(max_g,y,py)
# plot(log(vals),g_max,ylim=c(0,1))
# abline(h=0.5)


### CODE COMMENTED OUT
# par(mfrow=c(3,1))
# y=rnorm(n=100,mean=log(1),sd=log(10));
# y=10^y
# py=dnorm(log(y,base=10),mean=log(1),sd=log(10))
# 
# plot(y,py,type="h", lwd = 1)
# plot(log(y,base=10),py,type="h", lwd = 1,ylim=c(0,1))
# 
# plot(g,v(G=g,D=0,Y=y,P=py),type='n',ylim=c(-40,10))
# abline(h=0,col='gray')
# for(i in 1:length(d)){
#   lines(g,unlist(v(G=g,D=d[i],Y=y,P=py)),type='l',pch=16)
#   out<-unlist(v(G=g,D=d[i],Y=y,P=py))
#   points(g[which(out %in% max(out))],max(out),col='red',pch=16)
# }
# 
# #
# par(mfrow=c(3,1))
# y=rnorm(n=100,mean=log(1),sd=log(10));
# y=10^y
# py=dnorm(log(y,base=10),mean=log(1),sd=log(10))
# 
# plot(g,v(G=g,D=0,Y=y,P=py),type='n',ylim=c(-20,15))
# abline(h=0,col='gray')
# for(i in 1:length(d)){
#   lines(g,unlist(v(G=g,D=d[i],Y=y,P=py)),type='b',pch=16)
#   out<-unlist(v(G=g,D=d[i],Y=y,P=py))
#   points(g[which(out %in% max(out))],max(out),col='red',pch=16)
# }
# 
# out<-unlist(v(G=g,D=d[3],Y=y,P=py))
# g[which(out %in% max(out))];max(out)
# x<-seq(-7,7,by=.1)
# #plot(x,1 - ecdf(log(y,base=10))(x),ylim=c(0,1),type='l')
# 
# emp=1-ecdf(log(y,base=10))(x)
# gmax=g[which(out %in% max(out))]
# min(which(emp<gmax))
# 
# plot(y,py,type="h", lwd = 1)
# abline(v=y[min(which(emp<gmax))],col='red')
# plot(log(y,base=10),py,type="h", lwd = 1,ylim=c(0,1))
# abline(v=log(y,base=10)[min(which(emp<gmax))],col='red')
# 
# # own figure
# seqLen=7
# fromVal=0.01; toVal=200
# y=10^seq(log(fromVal,base=10), log(toVal,base=10), length.out = seqLen)
# y=c(.01,.02,.2,2,20,200)
# #y = c(.1,.2,.5,1,2,5,10)
# ly = log(y,base=10)
# py = c(0.15,.2,.3,.2,.1,.05)
# #py = dnorm(ly,mean=ly[3],sd=2)
# 
# par(mfrow=c(2,1))
# 
# plot(ly,py,type='h',lwd=2)
# 
# plot(g,v(G=g,D=0,Y=y,P=py),type='n',ylim=c(-1,.5))
# abline(h=0,col='gray')
# for(i in 1:length(d)){
#   lines(g,unlist(v(G=g,D=d[i],Y=y,P=py)),type='b',pch=16)
#   out<-unlist(v(G=g,D=d[i],Y=y,P=py))
#   points(g[which(out %in% max(out))],max(out),col='red',pch=16)
# }
# 
# 
# seqLen=7
# fromVal=0.01; toVal=200
# y=10^seq(log(fromVal,base=10), log(toVal,base=10), length.out = seqLen)
# y=c(.01,.02,.2,2,20,200)
# #y = c(.1,.2,.5,1,2,5,10)
# ly = log(y,base=10)
# py = c(0.15,.2,.3,.2,.1,.05)
# #py = dnorm(ly,mean=ly[3],sd=2)
# 
# par(mfrow=c(2,1))
# 
# plot(ly,py,type='h',lwd=2)
# 
# plot(g,v(G=g,D=0,Y=y,P=py),type='n',ylim=c(-1,.5))
# abline(h=0,col='gray')
# for(i in 1:length(d)){
#   lines(g,unlist(v(G=g,D=d[i],Y=y,P=py)),type='b',pch=16)
#   out<-unlist(v(G=g,D=d[i],Y=y,P=py))
#   points(g[which(out %in% max(out))],max(out),col='red',pch=16)
# }
# 
# 
# seqLen=100
# fromVal=0.01; toVal=100
# y=10^seq(log(fromVal,base=10), log(toVal,base=10), length.out = seqLen)
# #y=c(.01,.02,.2,2,20,200)
# #y = c(.1,.2,.5,1,2,5,10)
# ly = log(y,base=10)
# #py = c(0.15,.2,.3,.2,.1,.05)
# py = dnorm(ly,mean=ly[3],sd=2)
# 
# par(mfrow=c(3,1))
# 
# plot(y,py,type='h',lwd=2)
# plot(ly,py,type='h',lwd=2)
# 
# plot(g,v(G=g,D=0,Y=y,P=py),type='n',ylim=c(-20,2))
# abline(h=0,col='gray')
# for(i in 1:length(d)){
#   lines(g,unlist(v(G=g,D=d[i],Y=y,P=py)),type='b',pch=16)
#   out<-unlist(v(G=g,D=d[i],Y=y,P=py))
#   points(g[which(out %in% max(out))],max(out),col='red',pch=16)
# }
