
df=read.csv(file="~/Downloads/estimates-2011.csv",header=TRUE)
df=df[1:20,]
df

position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)
siteNames=unique(position$site)

library(MASS)

nsites = 20
nyears = 15

g1.min=min(df$g1)
g1.max=max(df$g1)
s2s3.min=min(df$s2*df$s3)
s2s3.max=max(df$s2*df$s2)
#mu = mu = rnorm(n=nsites,mean=25,sd=0)

sim <- function(rho=0,nsites=20){
  # Defines the  sequence for stochastic trials.
  reps=1000
  
  # a and b are hyperparameters of the gamma distribution 
  # that define both the expected value and variance.   
  # a = 6
  # b = 1.75
  
  # alpha and beta are hyperparameters of the beta distribution that define both the expected value and
  # variance.  
  alpha =  1
  beta =  1
  
  # Defines the temporal correlation between the two parameters.
  rho = rho
  
  # Generates standard multivariate normal data with correlation structure defined by rho.
  Z <- mvrnorm(n=reps,mu=c(0,0), 
               matrix(data=c(1,rho,rho,1),
                      nrow=2,ncol=2))
  
  # Apply the Normal CDF function to Z to obtain data that is uniform on the interval [0,1], but still correlated.
  U = pnorm(Z)
  
  # x is gamma distributed
  # X = qgamma(U[,1],shape=a,rate=b) 
  # range of s2s3

   X = qunif(U[,1],s2s3.min,s2s3.max) 
  
  # y is beta distributed
  #X <- cbind(X,qbeta(U[,2],shape1=alpha,shape2=beta) )
   #range of g1

   
  X <- cbind(X,qunif(U[,2],g1.min,g1.max) )
  
  # gamma marginal of multivariate X.
  # hist(X[,1])
  # beta marginal of multivariate X
  # hist(X[,2])
  
  # plot(X[,1],X[,2])
  
  nsites = nsites
  nyears = 10
  X[1:nsites,1:2]
  
  # # generate samples for each site with a different variance
  # d<-lapply(X[1:nsites,1],rnorm,n=nyears,mean=rnorm(n=nsites,mean=mu,sd=0))
  # # check that the variances are appropriately sampled
  # cbind(unlist(lapply(d, sd)),X[1:nsites,1])
  # 
  # d<-data.frame(d)
  # names(d) <- 1:20
  # 
  # # check variances again
  # cbind(apply(d,2, sd),X[1:nsites,1])
  # 
  # dim(d)
  # 
  # dt<-cbind(X[1:nsites,1:2],t(d))
  dt= X[1:nsites,1:2]
  return(dt)
}

pdf(
  "~/Dropbox/clarkiaSeedBanks/products/figures/analysis/simulation-1.pdf",
  height = 4, width = 8)
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

c0<-replicate(8,sim(rho=-.5,nsites=20))


layout.matrix <- matrix(c(1, 2, 3, 4, 5,
                          6,7,8, 9,9), nrow = 2, ncol = 5)

layout(mat = layout.matrix,
       heights = c(1,1), # Heights of the two rows
       widths = c(1,1,1,1, 3)) # Widths of the two columns

#layout.show(9)

#par(fig=c(0,10,5.5,10)/10)
# simulate 8 sample datasets


#par(mfrow=c(2,4))
#par(mar=c(4,4,2,1))
#par(fig=c(0,10,5.5,10)/10)
par(
    mai = c(.5, .5, 0.1, 0.1),
    mar=c(0.5, 0.5, 0.2, 0.2),
    oma = c(4, 4, 2, 0.2))
for(i in 1:dim(c0)[3]){
  plot(NA,NA,xlim=c(s2s3.min,s2s3.max),ylim=c(g1.min,g1.max),
       axes=FALSE,frame=TRUE,
       xlab='',ylab='')
  tmp=c0[,,i]
  points(tmp[,1],tmp[,2],cex=1,col=Okabe_Ito[i],pch=16)
  text(x=.45,y=.09,round(cor(tmp[,1],tmp[,2]),2),pos=4,cex=.75)
  ifelse(i%in%c(1,2),axis(2L),NA)
  ifelse(i%in%c(2,4,6,8),axis(1L),NA)
}
mtext("Seed survival", side = 1, outer = TRUE, line = 2.2,adj=.3)
mtext("Germination", side = 2, outer = TRUE, line = 2.2)
mtext("Simulated data", side = 3, outer = TRUE, line = 0.5, adj =0)
#dev.off()



# pdf(
#   "~/Dropbox/clarkiaSeedBanks/products/figures/analysis/simulation-2.pdf",
#   height = 4, width = 4)
# show correlation of 10,000 sample datasets
par(
  mai = c(.5, .5, 0.1, 0.1),
  mar=c(0.5, 0.5, 0.2, 0.2),
  oma = c(4, 4, 0.2, 0.2),new=TRUE)
c0<-replicate(10000,sim(rho=-.5,nsites=20))
cor.vec=c()
for(i in 1:dim(c0)[3]){
  tmp=c0[,,i]
  cor.vec[i]=cor(tmp[,1],tmp[,2])
}

hist(cor.vec,main='',col='gray80',border='white',breaks=20,xlim=c(-1,1),
     xlab="Distribution of simulated correlations",ylab='',yaxt='n')
abline(v=quantile(cor.vec,c(.025,.5,.975)),lty='dashed',col='orange',lwd=2)
mtext("Distribution of simulated correlations",side=1,line=2.5)
dev.off()




##########################################3
# Simulation for g1 and RS
##########################################3

library(MASS)

nsites = 20
nyears = 15

y1=sample(c(.01,2),15,replace=TRUE,prob=c(.1,.9))
y2=sample(c(.01,2),15,replace=TRUE,prob=c(.9,1.))
y1
y2

gsd.am <- function(x){
  x=x+.5
  n = length(x[!is.na(x)])
  mu = exp(mean(log(x),na.rm=TRUE))
  y <- exp(sqrt(sum((log(x/mu))^2,na.rm=TRUE)/(n-1)))
  return(y)
}

lowFitnessYears <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/output/lowFitnessYears.RDS")

lowFitness=lowFitnessYears %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(n=n())

df$rs=df$sigma*df$fec*df$phi

df=df %>%
  dplyr::left_join(lowFitness,by="site") %>%
  dplyr::mutate(n=ifelse(is.na(n),0,n))

gsd_vec = c()
for(i in 1:20){
  rs_seq=rep(df$rs[i],15-df$n[i])
  rs_seq=c(rs_seq,rep(0,df$n[i]))
  gsd_vec[i] = gsd.am(rs_seq)
}

gsd_vec

mu.rs=(log((df$rs)))
sd.rs=sd(log(df$rs))
rs.min=min(df$rs)
rs.max=max(df$rs)
gsd.min=min(gsd_vec)
gsd.max=max(gsd_vec)
g1.min=min(df$g1)
g1.max=max(df$g1)

rlnorm(15,log(df$rs[1]),sdlog=1)


mu.log = mean(c(1.7,.06))
sd.log = sd(c(1.7,.06))
#mu = mu = rnorm(n=nsites,mean=25,sd=0)

sim <- function(rho=0,nsites=20){
  # Defines the  sequence for stochastic trials.
  reps=1000
  
  # # a and b are hyperparameters of the gamma distribution
  # # that define both the expected value and variance.
  # a = 10
  # b = 1.75
  
  # Eckhart and Geber 2005 
  # average of 1.7 + .06
  # log(lifetimefitness) over 2 years
  a = mu.log
  b = sd.log
  # 
  # # alpha and beta are hyperparameters of the beta distribution that define both the expected value and
  # # variance.
  # alpha =  1
  # beta =  1
  # 
  # Defines the temporal correlation between the two parameters.
  rho = rho
  
  # Generates standard multivariate normal data with correlation structure defined by rho.
  Z <- mvrnorm(n=reps,mu=c(0,0), 
               matrix(data=c(1,rho,rho,1),
                      nrow=2,ncol=2))
  
  # Apply the Normal CDF function to Z to obtain data that is uniform on the interval [0,1], but still correlated.
  U = pnorm(Z)
  
  # x is gamma distributed
 # X = qgamma(U[,1],shape=a,rate=b) 
  #X = qlnorm(U[,1],meanlog=mu.log,sdlog=sd.log)
  #X = qlnorm(U[,1],meanlog=a,sdlog=b)
  X <- qunif(U[,1],gsd.min,gsd.max) 
  
  # y is beta distributed
  #X <- cbind(X,qbeta(U[,2],shape1=alpha,shape2=beta) )
  #range of g1
  X <- cbind(X,qunif(U[,2],g1.min,g1.max) )
  
  # gamma marginal of multivariate X.
  # hist(X[,1])
  # beta marginal of multivariate X
  # hist(X[,2])
  
  # plot(X[,1],X[,2])
  
  # nsites = nsites
  # nyears = 15
  # X[1:nsites,1:2]
  # 
  # # generate samples for each site with a different variance
  # d<-lapply(X[1:nsites,1],rnorm,n=nyears,mean=rnorm(n=nsites,mean=mu,sd=0))
  # # check that the variances are appropriately sampled
  # cbind(unlist(lapply(d, sd)),X[1:nsites,1])
  # 
  # d<-data.frame(d)
  # names(d) <- 1:20
  # 
  # # check variances again
  # cbind(apply(d,2, sd),X[1:nsites,1])
  # 
  # dim(d)

 # dt<-cbind(X[1:nsites,1:2],t(d))
  dt= X[1:nsites,1:2]
  return(dt)
}


Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

pdf(
  "~/Dropbox/clarkiaSeedBanks/products/figures/analysis/simulation-2.pdf",
  height = 4, width = 8)

layout.matrix <- matrix(c(1, 2, 3, 4, 5,
                          6,7,8, 9,9), nrow = 2, ncol = 5)

layout(mat = layout.matrix,
       heights = c(1,1), # Heights of the two rows
       widths = c(1,1,1,1, 3)) # Widths of the two columns

#layout.show(9)

#par(fig=c(0,10,5.5,10)/10)
# simulate 8 sample datasets

c0<-replicate(8,sim(rho=-.5,nsites=20))

#par(mfrow=c(2,4))
#par(mar=c(4,4,2,1))
#par(fig=c(0,10,5.5,10)/10)
par(
  mai = c(.5, .5, 0.1, 0.1),
  mar=c(0.5, 0.5, 0.2, 0.2),
  oma = c(4, 4, 2, 0.2))
for(i in 1:dim(c0)[3]){
  plot(NA,NA,xlim=c(0,c0[,1,i] %>% max+1),ylim=c(g1.min,g1.max),
       axes=FALSE,frame=TRUE,
       xlab='',ylab='')
  tmp=c0[,,i]
  points(tmp[,1],tmp[,2],cex=1,col=Okabe_Ito[i],pch=16)
  text(x=-.5,y=.09,round(cor(tmp[,1],tmp[,2]),2),pos=4,cex=.75)
  ifelse(i%in%c(1,2),axis(2L),NA)
 # ifelse(i%in%c(1,3,5,7),axis(3L),NA)
  ifelse(i%in%c(2,4,6,8),axis(1L),NA)
}
mtext("Geometric SD of per-capita RS", side = 1, outer = TRUE, line = 2.2,adj=.3)
mtext("Germination", side = 2, outer = TRUE, line = 2.2)
mtext("Simulated data", side = 3, outer = TRUE, line = 0.5, adj =0)

par(
  mai = c(.5, .5, 0.1, 0.1),
  mar=c(0.5, 0.5, 0.2, 0.2),
  oma = c(4, 4, 0.2, 0.2),new=TRUE)
c0<-replicate(10000,sim(rho=-.5,nsites=20))
cor.vec=c()
for(i in 1:dim(c0)[3]){
  tmp=c0[,,i]
  cor.vec[i]=cor(tmp[,1],tmp[,2])
}

hist(cor.vec,main='',col='gray80',border='white',breaks=20,xlim=c(-1,1),
     xlab="Distribution of simulated correlations",ylab='',yaxt='n')
abline(v=quantile(cor.vec,c(.025,.5,.975)),lty='dashed',col='orange',lwd=2)
mtext("Distribution of simulated correlations",side=1,line=2.5)
dev.off()