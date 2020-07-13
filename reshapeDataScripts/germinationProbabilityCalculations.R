# number of seeds at start of experiment
total = 100
# number intact but not germinated
intact = 60
# number germinated
germinant = 20
# number intact in october
october.intact = 30
# number intact in january 2
january.2.intact = 10
# number germinated in january 2
january.2.germ = 10


# method 1
prob.germination = function( n = intact, y = germinant, p.v = v ){
  p = y/(n*(p.v^(1/3))+y)
  return(p)
}

prob.s1 = function( n = intact, y =germinant, y.o = october.intact, p.v = v, t = total ){
  p = ((p.v^(1/3))*n+y)/(t)
  return(p)
}

prob.s2 = function( n = intact, y =germinant, y.o = october.intact, p.v = v ){
  p = (y.o*p.v)/(n*(p.v^(1/3)))
  return(p)
}

prob.s3 = function( n = intact, y =germinant, y.o = october.intact, y.j = january.2.intact, y.g2 = january.2.germ, p.v = p , t = total){
  
  p.vj = ifelse(p.v2 < p.v, p.v^(2/3)*p.v2^(1/3), p.v2^(1/3))
  
  p = ((p.vj)*y.j+y.g2)/(y.o*p.v)
  return(p)
}
 


# method 2
prob.germination2 = function( n = intact, y =germinant, p.v =v){
  g = y/(y+n)
  p = g/(1-((1-(p.v^(1/3)))*(1-g)))
  return(p)
}

# probability of being viable and intact
prob.s1.2 = function( n = intact, y =germinant, y.o = october.intact, p.v = v , t = total){
  p.i = (n+y)/t
  g = y/(y+n)
  p = ((1-g)*(p.v^(1/3))+g)*p.i
  return(p)
}

prob.s2.2 = function( n = intact, y =germinant, y.o = october.intact, p.v = v ){
  p.o = y.o/n
  p = p.o*p.v^(2/3)
  return(p)
}


prob.s3.2 = function( n = intact, y =germinant, y.o = october.intact, y.j = january.2.intact, y.g2 = january.2.germ, p.v = p , t = total){
  
  p.i = (y.j+y.g2)/t
  
  g = y.g2/(y.g2+y.j)

  ps1=prob.s1.2( n = n, y =y, y.o = y.o, p.v = p.v , t = t )
  pg1=prob.germination2( n = n, y = y, p.v = p.v)
  ps2=prob.s2.2(n = n, y = y, y.o = y.o, p.v = p.v)
  
  p.vj = ifelse(p.v2 < p.v, p.v^(2/3)*p.v2^(1/3), p.v2^(1/3))
  #p.vj = p.v^(2/3)*p.v2^(1/3)
  
  ps3 = prob.s1.2( n = y.j, y = y.g2, p.v = p.vj, t = t)
  pg2 = prob.germination2( n = y.j, y = y.g2, p.v = p.v2)
  
  #p = ((1-g)*(p.vj^(1/3))+g)*p.i*(1/(ps1*(1-pg1)*ps2))
  p = ((1-pg2)*(p.vj)+pg2)*ps3*(1/(ps1*(1-pg1)*ps2))
  
  return(p)
}

prob.s1.2()

# if PV2(p,r)<PV1(p,r)
# PVjan=(PV1(p,r)^.66667)*(PV2(p,r)^.33333);
# else
#   PVjan=PV2(p,r)^.33333;



prob.s3(n = intact, y =germinant, y.o = october.intact, y.j = january.2.intact, y.g2 = january.2.germ, p.v = p , t = total)
prob.s3.2(n = intact, y =germinant, y.o = october.intact, y.j = january.2.intact, y.g2 = january.2.germ, p.v = p , t = total)


# pdf(file="/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/appendix/comparison.pdf",
#     width=9,height=3)
par(mfrow=c(1,4))
p = seq(0,1,by=0.001)
plot(p, prob.s1(n = intact, y = germinant, y.o = october.intact, p.v = p ),
     xlim=c(0,1), ylim=c(0,1),
     cex = .5, pch = 16,
     xlab="Probability of viability",
     ylab="Probability of survival",type='l',
     main = "s1")
abline(h=(intact+germinant)/(100),lty='dashed')

p = seq(0,1,by=0.025)
points(p, prob.s1.2(n = intact, y = germinant, y.o = october.intact, p.v = p),
     xlim=c(0,1), ylim=c(0,1),
     cex = .5, pch = 4,col='red')

p = seq(0,1,by=0.001)
plot(p, prob.germination( n = intact, y = germinant, p.v = p ),
     xlim=c(0,1), ylim=c(0,1),
     cex = .5, pch = 16,
     xlab="Probability of viability",
     ylab="Probability of germination",type='l',
     main = "g1")
abline(h=germinant/(intact+germinant),lty='dashed')

p = seq(0,1,by=0.025)
points(p, prob.germination2( n = intact, y = germinant, p.v = p ),
     xlim=c(0,1), ylim=c(0,1),
     cex = .5, pch = 4,col='red')


p = seq(0,1,by=0.001)
plot(p, prob.s2(n = intact, y = germinant, y.o = october.intact, p.v = p ),
     xlim=c(0,1), ylim=c(0,1),
     cex = .5, pch = 16,
     xlab="Probability of viability",
     ylab="Probability of survival",type='l',
     main = "s2")
abline(h=(october.intact)/(intact),lty='dashed')

p = seq(0,1,by=0.025)
points(p, prob.s2.2(n = intact, y = germinant, y.o = october.intact, p.v = p),
     xlim=c(0,1), ylim=c(0,1),
     cex = .5, pch = 4,col='red')

# year 2 viability
p.v2 = 1

#par(mfrow=c(1,1))
p = seq(0,1,by=0.001)
plot(p, prob.s3(n = intact, y =germinant, y.o = october.intact, y.j = january.2.intact, y.g2 = january.2.germ, p.v = p , t = total ),
     xlim=c(0,1), ylim=c(0,1),
     cex = .5, pch = 16,
     xlab="Probability of viability",
     ylab="Probability of survival",type='l',
     main='s3')
abline(h=(january.2.intact+january.2.germ)/(october.intact),lty='dashed')
abline(h=1,lty='solid')
# abline(v=2/3)

p = seq(0,1,by=0.025)
points(p, prob.s3.2(n = intact, y =germinant, y.o = october.intact, y.j = january.2.intact, y.g2 = january.2.germ, p.v = p , t = total),
       xlim=c(0,10), ylim=c(0,1),
       cex = .5, pch = 4,col='red')
#dev.off()



