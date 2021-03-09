# Validate geometric SD calculation ---------------------------------------


# https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) /( length(x)))
}

# sample implementation
gsd.sample <- function(x){
  n = length(x[!is.na(x)])
  mu = (prod(x,na.rm=TRUE))^(1/n)
  y <- exp(sqrt(sum((log(x/mu))^2,na.rm=TRUE)/(n-1)))
  return(y)
}

gsd.am <- function(x){
  n = length(x[!is.na(x)])
  mu = exp(mean(log(x),na.rm=TRUE))
  y <- exp(sqrt(sum((log(x/mu))^2,na.rm=TRUE)/(n-1)))
  return(y)
}

# lognormal implementation
gsd <- function(x){
  y <- exp(sd(log(x),na.rm=TRUE))
  return(y)
}

# Test geometric standard deviation function ------------------------------------------------------

x <- c(rlnorm(n=10000,meanlog = 4 , sdlog = 1),NA)
y <- c(rpois(n=10000,lambda=x),NA)
gsd(x)
gsd.sample(x)
gsd.am(x)

# Estimate of the GSD seems particularly sensitive around .5
# becoming more sensitive as the mean on the log scale increases
# such that the estimate of the GSD seems to depend on the number of NAs?
delta=seq(0,2,by=.01)
tmp=c()
for(i in 1:length(delta)){tmp[i]=gsd(y+delta[i])}
plot(delta,tmp,type='l')
abline(h=gsd(y),v=.5)

tmp=c()
for(i in 1:length(delta)){tmp[i]=gsd.sample(y+delta[i])}
lines(delta,tmp,col='red')

tmp=c()
for(i in 1:length(delta)){tmp[i]=gsd.am(y+delta[i])}
lines(delta,tmp,col='red')


x <- rnorm(n=1000,mean = 20 , sd = 1)
gsd(x)
gsd.am(x)

n=10000
b = sample(10:20,n/2,replace=TRUE)
a = c(rep(15,n/2))
var(c(a,b))
var(b)
