#################################################################################
################################################################################
################################################################################
# ESS Germination Strategies in Randomly Varying Environments. II. Reciprocal Yield-Law Models
# Stephen Ellner
# Theoretical Population Biology
# 1985
# Vol 28, 80-116
#
# Scripts by Gregor Siegmund
# fausto.siegmund@gmail.com
# last updated 02-24-2020
#################################################################################
#################################################################################
#################################################################################


# function for total yield as a function
# of initial seedling density
totalYield <- function(K=k, b=b, X=x){
  f = (X*K)/(K/b+X)
  return(f)
}

plot(x=c(),y=c(),type='n',xlim=c(0,1),ylim=c(0,1),
     ylab="Total seed yield [xY(x)]",
     xlab="Number of seedlings [x]")
abline(h=0.8)
abline(a=0,b=20)

x=seq(0,1,by=0.01)
lines(x,totalYield(K=.8,b=20,X=x))

## model

# # equation year-by year
# e<-function(G=g,D=d,Y=y,P=py){
#   E=(1-P)*log10((1-G)*(1-D)+G*Y[1])+(P)*log10((1-G)*(1-D)+G*Y[2])
#   return(E)
# }

t = 1000
alpha.t = rnorm(n=t,mean=0,sd=1)

layout(matrix(c(1,2), 2, 1, byrow = TRUE))

# Panel 1
k.mean = .5
k.sd = 1

k.t = k.mean+k.sd*alpha.t
plot(density(k.t), xlim=c(-8,8))
# Panel 2
k.mean = .5
k.sd = 2

k.t = k.mean+k.sd*alpha.t
lines(density(k.t),lty='dotted')

plot(density((alpha.t)))


###
t = 1000
alpha.t = rnorm(n=t,mean=0,sd=1)

layout(matrix(c(1,2), 2, 1, byrow = TRUE))

# Panel 1
k.mean = .5
k.sd = 1

k.t = k.mean+k.sd*alpha.t
plot(density(k.t), xlim=c(-8,8))
# Panel 2
k.mean = .5
k.sd = 2

k.t = k.mean+k.sd*alpha.t
lines(density(k.t),lty='dotted')

plot(density((alpha.t)))


t=100
alpha.t = rnorm(n=t,mean=0,sd=1)

yFun<-function(Y.bar = y,y.sd=y.sd,alpha.t=alpha.t){
  y.t = Y.bar+y.sd*alpha.t
  return(y.t)
}

yFun(Y.bar=.1,y.sd=1,alpha.t=alpha.t)

# equation 4 from Cohen 1966
f<-function(G=g,S=s,Y=y){
  Fitness=log10((1-G)*(S)+G*Y)
  return(Fitness)
}

# script to calculate long term population growth rate
# for different combinations of parameters
v<-function(G=g,S=s,Y=y){
  vec <- c()
  for(i in 1:length(g)){
    vec[i]<-sum(f(G=g[i],S=S,Y=Y))
  }
  return(vec)
}

# script to calculate optimal g
max_g <- function(y,s=s){
  out<-unlist(v(G=g,S=s,Y=y))
  g[which(out %in% max(out))];max(out)
  x<-seq(-7,7,by=.001)
  
  emp=1-ecdf(log(y,base=10))(x)
  gmax=g[which(out %in% max(out))]
  return(gmax)
}
g=.1
max_g(y=yFun(Y.bar=.1,y.sd=1,alpha.t=alpha.t),
.5)
