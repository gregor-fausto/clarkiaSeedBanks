decay = function(t,rho=1){
  1*exp(-rho*t)
}

plot(seq(0,5,length.out=10),seq(0,1,length.out=10),type='n')
x=seq(0,5,by=.1)
rho = seq(0,100,length.out=100)
for(i in 1:100){
lines(x,decay(x,rho[i]))  
}

rho=1
plot(seq(0,5,by=.1),decay(seq(0,5,by=.1),rho=1),ylim=c(0,1),type='l')
s1=decay(1,rho=1)
g1 = .2
s2=decay(1,rho=1)*(1-g1)*decay(2,rho=1)

points(1,s1)
points(1,g1)
points(2,s2)
lines(seq(0,5,by=.1),decay(seq(0,5,by=.1),rho=1),lty=2)


continuousExp = function(t,alpha=1,beta=1){
  1/(1+beta*t)^alpha
}

plot(seq(0,5,length.out=10),seq(0,1,length.out=10),type='n')
x=seq(0,5,by=.1)
alpha = seq(0,10,length.out=10)
beta = seq(0,10,length.out=10)
for(i in 1:10){
  lines(x,continuousExp(x,alpha[i],beta[i]))  
}


weibullRes = function(t,alpha=1,beta=1){
  exp(-(t/beta)^alpha)
}

plot(seq(0,5,length.out=10),seq(0,1,length.out=10),type='n')
x=seq(0,5,by=.1)
alpha = 1
beta = seq(.1,10,length.out=10)
for(i in 1:10){
  lines(x,weibullRes(x,alpha,beta[i]))  
}



plot(decay(seq(0,5,by=.1),rho=1),ylim=c(0,1),type='l')
lines(continuousExp(seq(0,5,by=.1),alpha=1,beta=1),ylim=c(0,1),type='l',lty='dotted')
abline(h=1/.1)



compoundExp = function(t, kappa = 1, rho = 1){
  (kappa*(kappa/rho)^kappa)/((t + (kappa/rho))^(kappa+1))
}

lines(compoundExp(seq(0,10,by=.1),kappa=100,rho=1),col="red")

weibull = function(t, kappa = 1, rho = 1){
  (kappa*rho*(rho*t)^(kappa-1))*exp(-(rho*t)^kappa)
}

lines(weibull(seq(0,10,by=.1),kappa=1,rho=1),col="blue")


loglog = function(t, kappa = 1, rho = 1){
  (kappa*(t^(kappa-1))*(rho)^(kappa))/(1+(t*rho)^kappa)^2
}

lines(loglog(seq(0,10,by=.1),kappa=1,rho=1),col="purple")
