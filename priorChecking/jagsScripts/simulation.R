par(mfrow=c(1,2))
b=exp(c(-3,-2,-1,0,1,2,3,4,5))
x=seq(0,1,by=.1)
plot(x,exp(-b[1]*x),type='n',ylim=c(0,1))
for(i in 1:length(b)){
  lines(x,exp(-b[i]*x))
}
#min(exp(-b*x))
exp(0*x)
exp(-1*x)
exp(-2*x)
exp(-3*x)

## MEAN LIFETIME
# tau = 1/lambda
# the lambda parameter can be used to
# calculate the mean lifetime of seeds
latent = seq(-4,4,by=.5)
b = exp(latent)
latent
b
# so for time scaled to 36 months, the average lifetime is
1/b
# in months
1/b*36
# in years
1/b*36/12

(log(2)/b)*36/12

f<-function(x) 1/(x^2)
f(c(1/8,.25,1/3,.5,1,2))
# so a b of ~.05 is an average lifetime of 60 years
# and a b of 150 is an average lifetime of .24 of a month, or 1 week


b2=b/36
x=seq(0,36,by=.1)
plot(x,exp(-b2[1]*x),type='n',ylim=c(0,1))
for(i in 1:length(b2)){
  lines(x,exp(-b2[i]*x))
}
#min(exp(-b2*x))

# x=seq(-1,1)
# plot(x,exp(x),type='l')

par(mfrow=c(1,2))
hist(exp(rnorm(1000,0,1)))
hist(exp(rnorm(1000,0,1/sqrt(10))))
hist(exp(rnorm(1000,0,.1)))

     
#par(mfrow=c(2,2))
eta_s=rnorm(10000,0,1)
alpha = rlnorm(10000,0,.5)
#hist(exp(-(eta_s/alpha)  ),breaks=100)
plot(alpha,eta_s);text(x=4,y=-4,signif(cor(eta_s,alpha),4))

#hist(exp(-(1/exp(-eta_s/alpha))^alpha),breaks=100)

par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,1),ylim=c(0,1))
t=seq(0,1,by=.01)

for(i in 1:100){
  index<-sample(1:1000,1)
  lines(t,exp(-(t/exp(-eta_s[i]/alpha[i]))^alpha[i]),
        col='gray',lwd=.5)
}

eta_s=rnorm(10000,0,1)
alpha = rgamma(10000,2,2)
#hist(exp(-(eta_s/alpha)  ),breaks=100)
#plot(alpha,eta_s);text(x=4,y=-4,signif(cor(eta_s,alpha),4))

#hist(exp(-(1/exp(-eta_s/alpha))^alpha),breaks=100)


par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,1),ylim=c(0,1))
t=seq(0,1,by=.01)

for(i in 1:100){
  index<-sample(1:1000,1)
  lines(t,exp(-(t/exp(-eta_s[i]/alpha[i]))^alpha[i]),
        col='gray',lwd=.5)
}

eta=rnorm(10000,0,1)
alpha=rlnorm(10000,0,1)

eta = 1*eta
alpha=exp()

eta_surv[i] ~ dnorm(mu_s[siteSurvival[i],yearSurvival[i]],tau_s[siteSurvival[i],yearSurvival[i]])
inv.b[i] <- exp(-(eta_surv[i])/(a[siteSurvival[i]]))
# model without germination pieces
mu_survival[i] <- exp(-(months[i]/inv.b[i])^a[siteSurvival[i]])
mu[i] <- theta_c[siteSurvival[i],compIndex[i]]*exp(-(months[i]/inv.b[i])^a[siteSurvival[i]])

transformed parameters {
  real y;
  vector[9] x;
  
  y = 3.0 * y_raw;
  x = exp(y/2) * x_raw;
}
model {
  y_raw ~ std_normal(); // implies y ~ normal(0, 3)
  x_raw ~ std_normal(); // implies x ~ normal(0, exp(y/2))
}

parameters{
  real mu;
  real<lower=0> sig;
  real log_y_raw;
}

transformed parameters {
  real y = exp(log_y_raw * sig + mu);
}

model {
  log_y_raw ~ normal(0, 1);
  //do something more with y
}