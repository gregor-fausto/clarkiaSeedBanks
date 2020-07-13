v1 = .6
v2 = .8
#v1 = v2^(1/2)
vx = v2/v1


plot(c(0,1,2),c(1,v1,v2),ylim=c(0,1),pch=16)

lines(x=c(0,1),y=c(1,v1))
lines(x=c(0,2),y=c(1,v2),lty=2)
lines(x=c(1,2),y=c(v1,v2),lty=3)

lines(x=c(0,1),y=c(1,(v1+v2^(1/2))/2),lty=3)
lines(x=c(0,2),y=c(1,(v1+v2^(1/2))/2)^2,lty=3)

monthly<-function(x,p){
  tmp<-x/12
  p^(tmp)
}

ms <- 1:12
#points(ms/12,monthly(ms,p=v1))
points(ms/12,monthly(ms,p=((v1+v2^(1/2))/2)))

ms <- 1:24
monthly<-function(x,p){
  tmp<-x/24
  p^(tmp)
}
#points(ms/12,monthly(ms,p=v2))
points(ms/12,monthly(ms,p=((v1+v2^(1/2))/2)^(2)),pch=16,cex=.5)

(v1+v2^(1/2))/2
0.6928203
0.8944272