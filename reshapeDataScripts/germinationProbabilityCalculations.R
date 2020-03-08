# number of seeds at start of experiment
total = 100
# number intact but not germinated
intact = 60
# number germinated
germinant = 10
# number intact in october
october.intact = 20
# number intact in january 2
january.intact = 5

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

prob.s3 = function( n = intact, y =germinant, y.o = october.intact, y.j = january.intact, p.v = p , t = total){
  s1 = prob.s1(n = n, y =y, y.o = y.o, p.v = p.v, t = t )
  g1c = (1-prob.germination(n = n, y = y, p.v = p.v))
  s2 = prob.s2(n=n, y=y, y.o=y.o, p.v=p.v)
  s3 = january.intact/october.intact
  p = s1*g1c*s2*s3
  return(p)
}
# 
# I=find(strcmp(Pop,PopNow) & Round==r & Age==a & ~isnan(JanIntact)); % Note: JanTot,Germ,Intact either all estimated or all missing
# NBs3(p,r)=length(I);
# % est. total no. of ungerm. seeds still alive after 1 yr. - in Oct.
# SSs3(p,r)=100*s1(p,r)*(1-g1(p,r))*s2(p,r)*NBs3(p,r); 
# if PV2(p,r)<PV1(p,r)
# PVjan=(PV1(p,r)^.66667)*(PV2(p,r)^.33333);
# else
#   PVjan=PV2(p,r)^.33333;
# end
# SSg2(p,r)=sum( JanGerm(I)+ PVjan*JanIntact(I) ); % SS for estimating g2 is no. of VIABLE seeds in Jan.
# s3(p,r)=min(SSg2(p,r)/SSs3(p,r),1);
# g2(p,r)=sum( JanGerm(I) )/SSg2(p,r);     


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


par(mfrow=c(4,1))
p = seq(0,1,by=0.01)
plot(p, prob.s1(n = intact, y = germinant, y.o = october.intact, p.v = p ),
     xlim=c(0,1), ylim=c(0,1),
     cex = .5, pch = 16,
     xlab="Probability of viability",
     ylab="Probability of survival",type='l')
abline(h=(intact+germinant)/(100),lty='dashed')

p = seq(0,1,by=0.05)
points(p, prob.s1.2(n = intact, y = germinant, y.o = october.intact, p.v = p),
     xlim=c(0,1), ylim=c(0,1),
     cex = .5, pch = 4,col='red')

p = seq(0,1,by=0.01)
plot(p, prob.germination( n = intact, y = germinant, p.v = p ),
     xlim=c(0,1), ylim=c(0,1),
     cex = .5, pch = 16,
     xlab="Probability of viability",
     ylab="Probability of germination",type='l')
abline(h=germinant/(intact+germinant),lty='dashed')

p = seq(0,1,by=0.05)
points(p, prob.germination2( n = intact, y = germinant, p.v = p ),
     xlim=c(0,1), ylim=c(0,1),
     cex = .5, pch = 4,col='red')


p = seq(0,1,by=0.01)
plot(p, prob.s2(n = intact, y = germinant, y.o = october.intact, p.v = p ),
     xlim=c(0,1), ylim=c(0,1),
     cex = .5, pch = 16,
     xlab="Probability of viability",
     ylab="Probability of survival",type='l')
abline(h=(october.intact)/(intact),lty='dashed')

p = seq(0,1,by=0.05)
points(p, prob.s2.2(n = intact, y = germinant, y.o = october.intact, p.v = p),
     xlim=c(0,1), ylim=c(0,1),
     cex = .5, pch = 4,col='red')

p = seq(0,1,by=0.01)
plot(p, prob.s3(n = intact, y =germinant, y.o = october.intact, y.j = january.intact, p.v = p, t = total ),
     xlim=c(0,1), ylim=c(0,1),
     cex = .5, pch = 16,
     xlab="Probability of viability",
     ylab="Probability of survival",type='l')
abline(h=january.intact/(total),lty='dashed')

# p = seq(0,1,by=0.05)
# points(p, prob.s2.2(n = intact, y = germinant, y.o = october.intact, p.v = p),
#        xlim=c(0,10), ylim=c(0,1),
#        cex = .5, pch = 4,col='red')



