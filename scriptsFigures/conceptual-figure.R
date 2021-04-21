######################### 
## Code for conceptual figure describing seed bank data
## Gregor Siegmund
#########################


######################### 
## Parameter values
#########################
# number of seeds starting the trials
n = 100
# scale parameter
lambda = 1
# translate scale parameter to terms for Weibull
beta=1/lambda
# shape parameter
alpha=1

# conditional probability of germination (all *intact* seeds have x% germination probability)
g = c(.25,.25,.25)

# viability at the end of the year for age 1, 2, and 3 seeds
v = c(.9,.8,.7)

######################### 
## Weibull survival function
#########################

f.exp = function(t,lam=lambda){
  return(exp(-(t/beta)^alpha))
}

######################### 
## Calculate age-specific probabilities and survivorship schedule
#########################

## sample times
t.sample = c(0,4,4,12,16,16,24,28,28,36)/36

# Calculate age-specific survival probability
Oct_0 = 1

Janpre_1 = f.exp(t=t.sample[2])
Jangerm_1 = g[1]
Janpost_1 = (1-g[1])
Oct_1 = f.exp(t=t.sample[4])/f.exp(t=t.sample[3])

Janpre_2 = f.exp(t=t.sample[5])/f.exp(t=t.sample[4])
Jangerm_2 = g[2]
Janpost_2 = (1-g[2])
Oct_2 = f.exp(t=t.sample[7])/f.exp(t=t.sample[6])

Janpre_3 = f.exp(t=t.sample[8])/f.exp(t=t.sample[7])
Jangerm_3 = g[3]
Janpost_3 = (1-g[3])
Oct_3 = f.exp(t=t.sample[10])/f.exp(t=t.sample[9])

# Calculate survivorship schedule
l.x=c()
l.x[1]=1
l.x[2]=Janpre_1
l.x[3]=Janpost_1*Janpre_1
l.x[4]=Oct_1*Janpost_1*Janpre_1
l.x[5]=Janpre_2*Oct_1*Janpost_1*Janpre_1
l.x[6]=Janpost_2*Janpre_2*Oct_1*Janpost_1*Janpre_1
l.x[7]=Oct_2*Janpost_2*Janpre_2*Oct_1*Janpost_1*Janpre_1
l.x[8]=Janpre_3*Oct_2*Janpost_2*Janpre_2*Oct_1*Janpost_1*Janpre_1
l.x[9]=Janpost_3*Janpre_3*Oct_2*Janpost_2*Janpre_2*Oct_1*Janpost_1*Janpre_1
l.x[10]=Oct_3*Janpost_3*Janpre_3*Oct_2*Janpost_2*Janpre_2*Oct_1*Janpost_1*Janpre_1
l.x[11]=NA

# unconditional germination probabilities in black
g.uncond = c(
  Jangerm_1*Janpre_1,
  Jangerm_2*Janpre_2*Oct_1*Janpost_1*Janpre_1,
  Jangerm_3*Janpre_3*Oct_2*Janpost_2*Janpre_2*Oct_1*Janpost_1*Janpre_1)

######################### 
## Calculate survivorship, pmf, and hazard
#########################
x_i = 0:9
survival = l.x[1:10]
pmf = c(1,l.x[1:9]-l.x[2:10])
hazard = c(NA,pmf[2:10]/survival[1:9])
# plot(x_i,hazard,ylim=c(0,1),
#      pch=16,col='gray50')

persistence.df = data.frame(x = x_i,
                            survival = survival,
                            pmf = pmf,
                            hazard = hazard)

s = 1:6
hazard = pmf[c(2,4,5,7,8,10)]/survival[c(1,3,4,6,7,9)]
#plot(s,hazard)

persistence.model.df = data.frame(s = s, p.s = 1 - hazard)


x_i = 0:9
# l.x.viability = c(l.x[1],
#                   l.x[2]*(g[1] + (1-g[1])*v[1]^(1/3)),
#                   l.x[3]*v[1]^(1/3),
#                   l.x[4]*v[1],
#                   #   l.x[5]*(g[2]+(1-g[2])*v[1]*(ifelse(v[2]/v[1]>1,1,v[2]/v[1]))^(1/3)),
#                   l.x[5]*(g[2]+(1-g[2])*v[1]*(v[2]/v[1])^(1/3)),
#                   #   l.x[6]*v[1]*(ifelse(v[2]/v[1]>1,1,v[2]/v[1]))^(1/3),
#                   l.x[6]*v[1]*(v[2]/v[1])^(1/3),
#                   l.x[7]*v[2],
#                   #   l.x[8]*(g[3]+(1-g[3])*v[2]*(ifelse(v[3]/v[2]>1,1,v[3]/v[2]))^(1/3)),
#                   l.x[8]*(g[3]+(1-g[3])*v[2]*(v[3]/v[2])^(1/3)),
#                   #  l.x[9]*v[2]*(ifelse(v[3]/v[2]>1,1,v[3]/v[2]))^(1/3),
#                   l.x[9]*v[2]*(v[3]/v[2])^(1/3),
#                   l.x[10]*v[3])

# conditional probability of germination

g.v=c()
g.v[1] = g[1]/(1-(1-v[1]^(1/3))*(1-g[1]))
g.v[2] = g[2]/(1-(1-(v[1]*(v[2]/v[1])^(1/3)))*(1-g[2]))
g.v[3] = g[3]/(1-(1-(v[2]*(v[3]/v[2])^(1/3)))*(1-g[3]))

# Calculate age-specific survival probability
Oct_0.v = 1

Janpre_1.v = f.exp(t=t.sample[2])*(g[1]+(1-g[1])*v[1]^(1/3))
Jangerm_1.v = g.v[1]
Janpost_1.v = (1-g.v[1])
Oct_1.v = (f.exp(t=t.sample[4])*v[1])/(f.exp(t=t.sample[3])*(v[1])^(1/3))

Janpre_2.v =( f.exp(t=t.sample[5])*(g[2]+(1-g[2])*v[1]*(v[2]/v[1])^(1/3)))/(f.exp(t=t.sample[4])*v[1])
Jangerm_2.v = g.v[2]
Janpost_2.v = (1-g.v[2])
Oct_2.v = (f.exp(t=t.sample[7])*v[2])/(f.exp(t=t.sample[6])*v[1]*(v[2]/v[1])^(1/3))

Janpre_3.v = (f.exp(t=t.sample[8])*(g[3]+(1-g[3])*v[2]*(v[3]/v[2])^(1/3)))/(f.exp(t=t.sample[7])*v[2])
Jangerm_3.v = g.v[3]
Janpost_3.v = (1-g.v[3])
Oct_3.v = (f.exp(t=t.sample[10])*v[3])/(f.exp(t=t.sample[9])*v[2]*(v[3]/v[2])^(1/3))

# Calculate survivorship schedule
l.x.viability=c()
l.x.viability[1]=1
l.x.viability[2]=Janpre_1.v
l.x.viability[3]=Janpost_1.v*Janpre_1.v
l.x.viability[4]=Oct_1.v*Janpost_1.v*Janpre_1.v
l.x.viability[5]=Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v
l.x.viability[6]=Janpost_2.v*Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v
l.x.viability[7]=Oct_2.v*Janpost_2.v*Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v
l.x.viability[8]=Janpre_3.v*Oct_2.v*Janpost_2.v*Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v
l.x.viability[9]=Janpost_3.v*Janpre_3.v*Oct_2.v*Janpost_2.v*Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v
l.x.viability[10]=Oct_3.v*Janpost_3.v*Janpre_3.v*Oct_2.v*Janpost_2.v*Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v
l.x.viability[11]=NA


survival = l.x.viability[1:10]
pmf = c(1,l.x.viability[1:9]-l.x.viability[2:10])
hazard = c(NA,pmf[2:10]/survival[1:9])
points(x_i,hazard,
       ylim=c(0,1),
       pch=16,col='orange')

viability.df = data.frame(x = x_i,
                            survival = survival,
                            pmf = pmf,
                            hazard = hazard)

s = 1:6
hazard = pmf[c(2,4,5,7,8,10)]/survival[c(1,3,4,6,7,9)]
plot(s,hazard)

viability.model.df = data.frame(s = s, p.s = 1 - hazard)

# calculate conditional germination probabilities (gs)
g.cond1 = c(g[1]/(1-(1-v[1]^(1/3))*(1-g[1])),
           g[2]/(1-(1-v[1]*(v[2]/v[1])^(1/3))*(1-g[2])),
           g[3]/(1-(1-v[2]*(v[3]/v[2])^(1/3))*(1-g[3])))

g.cond2 = c(g.uncond[1]/l.x.viability[2],
            g.uncond[2]/l.x.viability[5],
            g.uncond[3]/l.x.viability[8])

g.v;g.cond1;g.cond2;g.uncond

g.uncond = c(
  Jangerm_1*Janpre_1,
  Jangerm_2*Janpre_2*Oct_1*Janpost_1*Janpre_1,
  Jangerm_3*Janpre_3*Oct_2*Janpost_2*Janpre_2*Oct_1*Janpost_1*Janpre_1)
g.uncond.v = c(
  Jangerm_1.v*Janpre_1.v,
  Jangerm_2.v*Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v,
  Jangerm_3.v*Janpre_3.v*Oct_2.v*Janpost_2.v*Janpre_2.v*Oct_1.v*Janpost_1.v*Janpre_1.v)
  
g.uncond;g.uncond.v



######################### 
## Figure version 1
#########################

par(mfrow=c(1,3))

# time sequence
t = seq(0,1,by=.01)

## Panel A
## Survivorship schedule and unconditional germination 
# plot(t*36,f.exp(t=t),
#      xlim=c(0,40),ylim=c(0,100),
#      type='n',frame=FALSE,
#      xlab = "Time (months)",
#      ylab = "Seed counts")
# 
# segments(x0=t.sample[c(2,5,8)]*36,
#          y0=c(0,0,0),
#          y1=l.x[c(3,6,9)]*100+l.x[c(2,5,8)]*c(Jangerm_1,Jangerm_2,Jangerm_3)*100,
#          lty='dotted')
# 
# points(t.sample[c(3,4,6,7,9,10)]*36,
#        l.x[c(3,4,6,7,9,10)]*100,
#        col='black',cex=1,
#        pch=19)
# 
# points(t.sample[c(2,5,8)]*36,
#        c(l.x[2]*Jangerm_1,l.x[5]*Jangerm_2,l.x[8]*Jangerm_3)*100,
#        bg='white',cex=1,col='black',
#        pch=21)
# 
# legend(x = 25, y = 100,
#        pch = c(19,21,NA),
#        col = c('black','black'),
#        lty = c(FALSE,FALSE,3),
#        legend = c("Intact seeds", "Seedlings","Total seeds"),
#        cex=.75,
#        box.lty=0)


## Panel A
## Survivorship schedule and conditional germination 
plot(t*36,f.exp(t=t),
     xlim=c(0,40),ylim=c(0,100),
     type='n',frame=FALSE,
     xlab = "Time (months)",
     ylab = "Seed counts")

l.x.s=l.x[c(3,6,9)]*100+l.x[c(2,5,8)]*c(Jangerm_1,Jangerm_2,Jangerm_3)*100

# plot survivorship points on to curve

segments(x0=c(t.sample[c(3,6,9)]*36),
         y0=c(l.x.s),y1=c(l.x[c(3,6,9)]*100),
         lty=c('dotted'))

points(c(t.sample[1:4],t.sample[4],t.sample[5:7],t.sample[7],t.sample[8:10])*36,
       c(l.x[1:4],l.x[4],l.x[5:7],l.x[7],l.x[8:10])*100,
       col=c('black','gray',"black","gray","black","gray","black","gray","black","gray"),
       cex=c(1,1,1,2,1,1,1,2,1,1,1,1),
       pch=c(19,19,19,19,20,19,19,19,20,19,19,19))

legend(x = 10, y = 100,
       pch = c(19,19,NA,21),
       col = c('black','gray','black','orange'),
       lty = c(FALSE,FALSE,3,FALSE),
       legend = c("Start of interval", "End of interval",'Germination','Inferred viable seeds'),
       cex=.75,
       box.lty=0)

#points(t.sample[1:10]*36,l.x[1:10]*100,col='red')
points(t.sample[1:10]*36,l.x.viability[1:10]*100,col='orange',pch=21)


## Panel B
## Viability and inferred viability
plot(t*36,f.exp(t=t),
     xlim=c(0,40),ylim=c(0,1),
     type='n',frame=FALSE,
     xlab = "Time (months)",
     ylab = "Viability probability")



x0=seq(0,1,by=1/3)
lines(x=x0*12,v[1]^(x0),
      lty='dotdash',col='orange')

x1=seq(1,2,by=.01)
lines(x=x1*12,v[1]*(v[2]/v[1])^(x1-1),
      lty='dotdash',col='orange')

x2=seq(2,3,by=.01)
lines(x=x2*12,v[2]*(v[3]/v[2])^(x2-2),
      lty='dotdash',col='orange')

points(c(0,12,24,36),
       c(1,v),
       pch=19,
       col='orange')
points(c(4,16,28), 
       c(v[1]^(1/3), v[1]*(v[2]/v[1])^(1/3), v[2]*(v[3]/v[2])^(1/3)),
       pch=21,bg='white',
       col='orange')

legend(x = 0, y = .1,
       pch = c(19,21),
       col = c('orange','orange'),
       bg = c(NA,'white'),
       legend = c("Estimated viability", "Inferred viability"),
       cex=.75,
       box.lty=0)


## Panel C
## Conditional and unconditional germination
plot(c(1,2,3),g,
     ylim=c(0,1),
     pch=16, type = 'n', frame = FALSE,
     xlab = "Seed age (years)",
     ylab = "Germination probability")

# conditional germination probabilities in gray
points( x = c(1,2,3),
        y = g,
        pch = 19, col ='gray')

# unconditional germination probabilities in black
points( x = c(1,2,3),
        g.uncond,
        col='black', pch = 21,bg="white")

# germination probabilities conditional on viability in orange
points( x = c(1,2,3),
        g.cond1,
        col='orange', pch = 21,bg="white")

# check that plot is the same
# points(c(1,2,3),
#        c(l.x[2]*Jangerm_1,l.x[5]*Jangerm_2,l.x[8]*Jangerm_3),xlim=c(0,1),ylim=c(0,1),
#        col='red',cex=1,
#        pch=1)

legend(x = 1, y = 1,
       pch = c(21,19,21),
       col = c('black','gray','orange'),
       legend = c("Unconditional","Conditional on being intact","Conditional on being viable"),
       cex=.75,
       box.lty=0)

######################### 
## Figure version 2
#########################

par(mfrow=c(1,3))
## Panel A
## Survivorship schedule and conditional germination with survival function 
plot(t*36,f.exp(t=t),
     xlim=c(0,40),ylim=c(0,100),
     type='n',frame=FALSE,
     xlab = "Time (months)",
     ylab = "Seed counts")

l.x.s=l.x[c(3,6,9)]*100+l.x[c(2,5,8)]*c(Jangerm_1,Jangerm_2,Jangerm_3)*100

# plot survivorship points on to curve

segments(x0=c(t.sample[c(3,6,9)]*36),
         y0=c(l.x.s),y1=c(l.x[c(3,6,9)]*100),
         lty=c('dotted'))

points(c(t.sample[1:4],t.sample[4],t.sample[5:7],t.sample[7],t.sample[8:10])*36,
       c(l.x[1:4],l.x[4],l.x[5:7],l.x[7],l.x[8:10])*100,
       col=c('black','gray',"black","gray","black","gray","black","gray","black","gray"),
       cex=c(1,1,1,2,1,1,1,2,1,1,1,1),
       pch=c(19,19,19,19,20,19,19,19,20,19,19,19))

t1 = seq(0/36,t.sample[3]*36/36,by=.01)
lines(t1*36, f.exp(t=t1)*100)

t2 = seq(t.sample[3]*36/36,12/36,by=.01)
lines(t2*36, (1-g[1])*f.exp(t=t2)*100)

t3 = seq(12/36,t.sample[6]*36/36,by=.01)
lines(t3*36, (1-g[1])*f.exp(t=t3)*100)

t4 = seq(t.sample[6]*36/36,24/36,by=.01)
lines(t4*36, (1-g[1])*(1-g[2])*f.exp(t=t4)*100)

t5 = seq(24/36,t.sample[9]*36/36,by=.01)
lines(t5*36, (1-g[1])*(1-g[2])*f.exp(t=t5)*100)

t6 = seq(t.sample[9]*36/36,36/36,by=.01)
lines(t6*36, (1-g[1])*(1-g[2])*(1-g[3])*f.exp(t=t6)*100)

legend(x = 20, y = 100,
       pch = c(19,19,NA),
       col = c('black','gray','black'),
       lty = c(FALSE,FALSE,3),
       legend = c("Start of interval", "End of interval",'Germination'),
       cex=.75,
       box.lty=0)

## Panel B
## Survivorship schedule and conditional germination with survival function 
## Adding viability decay curves

plot(t*36,f.exp(t=t),
     xlim=c(0,36),ylim=c(0,100),
     type='n',frame=FALSE,
     xlab = "Time (months)",
     ylab = "Seed counts")

l.x.s=l.x[c(3,6,9)]*100+l.x[c(2,5,8)]*c(Jangerm_1,Jangerm_2,Jangerm_3)*100

# plot survivorship points on to curve

segments(x0=c(t.sample[c(3,6,9)]*36),
         y0=c(l.x.s),y1=c(l.x[c(3,6,9)]*100),
         lty=c('dotted'))

points(c(t.sample[1:4],t.sample[4],t.sample[5:7],t.sample[7],t.sample[8:10])*36,
       c(l.x[1:4],l.x[4],l.x[5:7],l.x[7],l.x[8:10])*100,
       col=c('black','gray',"black","gray","black","gray","black","gray","black","gray"),
       cex=c(1,1,1,2,1,1,1,2,1,1,1,1),
       pch=c(19,19,19,19,20,19,19,19,20,19,19,19))

t1 = seq(0/36,t.sample[3]*36/36,by=.001)
lines(t1*36, f.exp(t=t1)*100)

t2 = seq(t.sample[3]*36/36,12/36,by=.001)
lines(t2*36, (1-g[1])*f.exp(t=t2)*100)

t3 = seq(12/36,t.sample[6]*36/36,by=.001)
lines(t3*36, (1-g[1])*f.exp(t=t3)*100)

t4 = seq(t.sample[6]*36/36,24/36,by=.001)
lines(t4*36, (1-g[1])*(1-g[2])*f.exp(t=t4)*100)

t5 = seq(24/36,t.sample[9]*36/36,by=.001)
lines(t5*36, (1-g[1])*(1-g[2])*f.exp(t=t5)*100)

t6 = seq(t.sample[9]*36/36,36/36,by=.001)
lines(t6*36, (1-g[1])*(1-g[2])*(1-g[3])*f.exp(t=t6)*100)

# constant decay assumption

# viability at sampling in October
# multiply probability of survival in October by probability of viability
# points(t.sample[c(4,7,10)]*36,
#        v[1]*l.x[c(4,7,10)]*100,
#        col='red',pch=19)

# AGE 1
# viability in October of the first year
points(t.sample[c(4,7,10)]*36,
       l.x.viability[c(4,7,10)]*100,
       col='blue',pch=19,cex=1.5)

points(t.sample[c(3,6,9)]*36,
       l.x.viability[c(3,6,9)]*100,
       col='blue',pch=19,cex=1.5)

x0=seq(0,1/3,by=.01)
time.varying.lx = f.exp(t=x0/3)
lines(x=x0*12,time.varying.lx*v[1]^(x0)*100,
      lty='dotdash',col='orange')

x0=seq(1/3,1,by=.01)
time.varying.lx = f.exp(t=4/36)*(1-g[1])*(f.exp(t=x0/3)/f.exp(t=4/36))
lines(x=x0*12,time.varying.lx*v[1]^(x0)*100,
      lty='dotdash',col='orange')

x1=seq(1,1+1/3,by=.01)
time.varying.lx2 = (1-g[1])*f.exp(t=x1/3)

lines(x=x1*12,time.varying.lx2*v[1]*(v[2]/v[1])^(x1-1)*100,
      lty='dotdash',col='orange')

x1=seq(1+1/3,2,by=.01)
time.varying.lx3 = (1-g[1])*(1-g[2])*f.exp(t=x1/3)

lines(x=x1*12,time.varying.lx3*v[1]*(v[2]/v[1])^(x1-1)*100,
      lty='dotdash',col='orange')

x1=seq(2,2+1/3,by=.01)
time.varying.lx4 = (1-g[1])*(1-g[2])*f.exp(t=x1/3)

lines(x=x1*12,time.varying.lx4*v[2]*(v[3]/v[2])^(x1-2)*100,
      lty='dotdash',col='orange')

x1=seq(2+1/3,3,by=.01)
time.varying.lx5 = (1-g[1])*(1-g[2])*(1-g[3])*f.exp(t=x1/3)

lines(x=x1*12,time.varying.lx5*v[2]*(v[3]/v[2])^(x1-2)*100,
      lty='dotdash',col='orange')





points(t.sample[c(4)]*36,
       v[1]*(1-g[1])*f.exp(t=12/36)*100,
       col='red',pch=19)

# 1 full season
t.v1=seq(4/12,1,.01)

# probability of viability assuming constant decay
# over course of season
pv = v[1]^t.v1

# point in season from start to end
to.v1=seq(4/36,12/36,length.out=length(t.v1))

lines(t.v1*12,
      (v[1]^t.v1)*(1-g[1])*f.exp(t=to.v1)*100,
      col='red',lty='dashed')

# viability in January of the first year
points(t.sample[c(3)]*36,
       (v[1]^(1/3))*(1-g[1])*f.exp(t=4/36)*100,
       col='red',pch=21,bg='white')

v.jan=c()
v.jan[1] = v[1]^(1/3)

# AGE 2
# viability in October of the second year
points(t.sample[c(7)]*36,
       v[2]*(1-g[2])*(1-g[1])*f.exp(t=24/36)*100,
       col='red',pch=19)

# 1 full season
t.v1=seq(4/12,1,.01)

v.jan[2] = ifelse(v[2]<v[1],
       (v[1]^(2/3))*(v[2]^(1/3)),
       v[2]^(1/3))

# probability of viability assuming constant decay
# over course of season
pv = v[2]^t.v1

# point in season from start to end
to.v1=seq(16/36,24/36,length.out=length(t.v1))

lines(12+t.v1*12,
      (v[2]^t.v1)*(1-g[2])*(1-g[1])*f.exp(t=to.v1)*100,
      col='red',lty='dashed')



# viability in January of the second year
points(t.sample[c(6)]*36,
       v[2]^(1/3)*(1-g[2])*(1-g[1])*f.exp(t=16/36)*100,
       col='red',pch=21,bg='white')

# AGE 3
# viability in October of the third year
points(t.sample[c(10)]*36,
       v[3]*(1-g[3])*(1-g[2])*(1-g[1])*f.exp(t=36/36)*100,
       col='red',pch=19)

# 1 full season
t.v1=seq(4/12,1,.01)

# probability of viability assuming constant decay
# over course of season
pv = v[3]^t.v1

# point in season from start to end
to.v1=seq(28/36,36/36,length.out=length(t.v1))

lines(24+t.v1*12,
      (v[3]^t.v1)*(1-g[3])*(1-g[2])*(1-g[1])*f.exp(t=to.v1)*100,
      col='red',lty='dashed')

# viability in January of the third year
points(t.sample[c(9)]*36,
       v[3]^(1/3)*(1-g[3])*(1-g[2])*(1-g[1])*f.exp(t=28/36)*100,
       col='red',pch=21,bg='white')

## WHAT TO DO WITH INTERPOLATION

# points(t.sample[c(6)]*36,
#        v.jan[2]*(1-g[2])*(1-g[1])*f.exp(t=16/36)*100,
#        col='blue',pch=21,bg='white')



# points(c(t.sample[1:4],t.sample[4],t.sample[5:7],t.sample[7],t.sample[8:10])*36,
#        c(l.x[1:4],l.x[4],l.x[5:7],l.x[7],l.x[8:10])*100,
#        col=c('black','gray',"black","gray","black","gray","black","gray","black","gray"),
#        cex=c(1,1,1,2,1,1,1,2,1,1,1,1),
#        pch=c(19,19,19,19,20,19,19,19,20,19,19,19))

# 
# # conditional survival 
# l.x.s=l.x[c(3,6,9)]*100+l.x[c(2,5,8)]*c(Jangerm_1,Jangerm_2,Jangerm_3)*100
# 
# # plot survivorship points on to curve
# 
# segments(x0=c(3,15,27),
#          y0=c(l.x.s),y1=c(l.x[c(3,6,9)]*100),
#          lty=c('dotted'))
# 
# points(c(t.sample[1:4],t.sample[4],t.sample[5:7],t.sample[7],t.sample[8:10])*36,
#        c(l.x[1:4],l.x[4],l.x[5:7],l.x[7],l.x[8:10])*100,
#        col=c('black','gray',"black","gray","black","gray","black","gray","black","gray"),
#        cex=c(1,1,1,2,1,1,1,2,1,1,1,1),
#        pch=c(19,19,19,19,20,19,19,19,20,19,19,19))
# 
plot(t.sample[2]*36,NA,
     xlim=c(0,36),ylim=c(0,1),
     pch=16, type = 'n', frame = FALSE,
     xlab = "Time (months)",
     ylab = "Survival probability",
     main = "Unconditional probabilities")

intact.seeds.jan = l.x[c(3,6,9)]*100
germ.jan = c(g[1]*f.exp(t=4/36)*100,
             g[2]*(1-g[1])*f.exp(t=16/36)*100,
             g[3]*(1-g[2])*(1-g[1])*f.exp(t=28/36)*100)
viable.seeds.jan = c((v[1]^(1/3))*(1-g[1])*f.exp(t=4/36)*100,
                     v[2]^(1/3)*(1-g[2])*(1-g[1])*f.exp(t=16/36)*100,
                     v[3]^(1/3)*(1-g[3])*(1-g[2])*(1-g[1])*f.exp(t=28/36)*100)
total.seeds.jan = c(f.exp(t=4/36)*100,
                    (1-g[1])*f.exp(t=16/36)*100,
                    (1-g[2])*(1-g[1])*f.exp(t=28/36)*100)

intact.seeds.oct = l.x[c(4,7,10)]*100
viable.seeds.oct = c((v[1])*(1-g[1])*f.exp(t=12/36)*100,
                     v[2]*(1-g[2])*(1-g[1])*f.exp(t=24/36)*100,
                     v[3]*(1-g[3])*(1-g[2])*(1-g[1])*f.exp(t=36/36)*100)

# January
points(t.sample[c(3,6,9)]*36,total.seeds.jan/100,
       pch=21,col='black',bg='white')
points(t.sample[c(3,6,9)]*36,intact.seeds.jan/100,
       pch=19,col='gray')
points(t.sample[c(3,6,9)]*36,viable.seeds.jan/100,
       pch=19,col='black')

# October
points(t.sample[c(4,7,10)]*36,intact.seeds.oct/100,
       pch=19,col='gray')
points(t.sample[c(4,7,10)]*36,viable.seeds.oct/100,
       pch=19,col='black')

legend(x = 0, y = 1,
       pch = c(19,19,19,10),
       col = c('gray','black',"orange",'red'),
       legend = c("Unconditional intact seeds Jan", 
                  "Unconditional viable seeds Jan", 
                  "Unconditional intact seeds Oct",
                  "Unconditional viable seeds Oct"),
       cex=.75,
       box.lty=0)

plot(t.sample[2]*36,NA,
     xlim=c(0,36),ylim=c(0,1),
     pch=16, type = 'n', frame = FALSE,
     xlab = "Time (months)",
     ylab = "Survival probability",
     main = "Conditional probabilities")


p.intact = (intact.seeds.jan+germ.jan)/100
p.germ = germ.jan/(intact.seeds.jan+germ.jan)
p.oct = intact.seeds.oct/intact.seeds.jan

# January: s1
p.intact*(p.germ+(1-p.germ)*v^(1/3))
((viable.seeds.jan+germ.jan)/100)[1]

# January: g1
p.germ/(1-(1-v^(1/3))*(1-p.germ))
((germ.jan)/(intact.seeds.jan*v^(1/3)+germ.jan))[1]

# October: s2
(intact.seeds.oct*v)/(intact.seeds.jan*v^(1/3))
p.oct*v^(2/3)




s.cond1 = (intact.seeds.jan/100)*((germ.jan/total.seeds.jan)+(1-(germ.jan/total.seeds.jan))*v^(1/3) )
points(t.sample[c(3,6,9)]*36,s.cond1,
       pch=21,col='black',bg='white')
points(t.sample[c(4,7,10)]*36,(intact.seeds.oct/intact.seeds.jan)*v^(2/3),
       pch=19,col='gray')


legend(x = 0, y = 1,
       pch = c(19,19,19,10),
       col = c('gray','black',"orange",'red'),
       legend = c("Unconditional intact seeds Jan", 
                  "Unconditional viable seeds Jan", 
                  "Unconditional intact seeds Oct",
                  "Unconditional viable seeds Oct"),
       cex=.75,
       box.lty=0)

# 
# s.c = s.i*(g+(1-g)*v^1/3)
# s2.c = s2.i*v^2/3
# 
# plot(c(1,2,3),g,
#      ylim=c(0,1),
#      pch=16, type = 'n', frame = FALSE,
#      xlab = "Seed age (years)",
#      ylab = "Germination probability")
# 
# # conditional germination probabilities in gray
# points( x = c(1,2,3),
#         y = g,
#         pch = 19, col ='gray')
# 
# # unconditional germination probabilities in black
# g.uncond = c(
#   Jangerm_1*Janpre_1,
#   Jangerm_2*Janpre_2*Oct_1*Janpost_1*Janpre_1,
#   Jangerm_3*Janpre_3*Oct_2*Janpost_2*Janpre_2*Oct_1*Janpost_1*Janpre_1)
# 
# points( x = c(1,2,3),
#         g.uncond,
#         col='black', pch = 21,bg="white")
# 
# # conditional on viability germination probabilities
# g.cond = c(
#   Jangerm_1/(1-(1-v[1]^(1/3))*(1-Jangerm_1)),
#   Jangerm_2/(1-(1-v[2]^(1/3))*(1-Jangerm_2)),
#   Jangerm_3/(1-(1-v[2]^(1/3))*(1-Jangerm_3)))
#   points( x = c(1,2,3),
#         g.cond,
#         col='black', pch = 21,bg="red")
# 
# # points(c(1,2,3),
# #        c(l.x[2]*Jangerm_1,l.x[5]*Jangerm_2,l.x[8]*Jangerm_3),xlim=c(0,1),ylim=c(0,1),
# #        col='red',cex=1,
# #        pch=1)
# 
# legend(x = 2.25, y = 1,
#        pch = c(19,21),
#        col = c('gray','black'),
#        legend = c("Conditional", "Unconditional"),
#        cex=.75,
#        box.lty=0)
# 
# #########################
# ## TWO PANEL PLOT
# #########################
# 
# par(mfrow=c(1,1))
# t.sample = c(0,3,3,12,15,15,24,28,28,36)/36
# 
# # plot exponential survival function
# plot(t*36,f.exp(t=t),
#      xlim=c(0,40),ylim=c(0,210),
#      type='n',frame=FALSE,
#      xlab = "Time (months)",
#      ylab = "",
#      axes=FALSE)
# axis(1, c(0,10,20,30,40), col.ticks = 1)
# axis(2, c(0,20,40,60,80,100), col = NA, col.ticks = 1)
# segments(x0=-1.5,y0=0,y1=100)
# mtext('Seed counts',side=2,line=3,adj=.25,col='black',cex=1)
# 
# l.x.s=l.x[c(3,6,9)]*100+l.x[c(2,5,8)]*c(Jangerm_1,Jangerm_2,Jangerm_3)*100
# 
# # plot survivorship points on to curve
# 
# segments(x0=c(3,15,27),
#          y0=c(l.x.s),y1=c(l.x[c(3,6,9)]*100),
#          lty=c('dotted'))
# 
# points(c(t.sample[1:4],t.sample[4],t.sample[5:7],t.sample[7],t.sample[8:10])*36,
#        c(l.x[1:4],l.x[4],l.x[5:7],l.x[7],l.x[8:10])*100,
#        col=c('black','gray',"black","gray","black","gray","black","gray","black","gray"),
#        cex=c(1,1,1,2,1,1,1,2,1,1,1,1),
#        pch=c(19,19,19,19,20,19,19,19,20,19,19,19))
# 
# t1 = seq(0/36,3/36,by=.01)
# lines(t1*36, f.exp(t=t1)*100)
# 
# t2 = seq(3/36,12/36,by=.01)
# lines(t2*36, (1-g[1])*f.exp(t=t2)*100)
# 
# t3 = seq(12/36,15/36,by=.01)
# lines(t3*36, (1-g[1])*f.exp(t=t3)*100)
# 
# t4 = seq(15/36,24/36,by=.01)
# lines(t4*36, (1-g[1])*(1-g[2])*f.exp(t=t4)*100)
# 
# t5 = seq(24/36,27/36,by=.01)
# lines(t5*36, (1-g[1])*(1-g[2])*f.exp(t=t5)*100)
# 
# t6 = seq(27/36,36/36,by=.01)
# lines(t6*36, (1-g[1])*(1-g[2])*(1-g[3])*f.exp(t=t6)*100)
# 
# legend(x = 30, y = 100,
#        pch = c(19,19,NA),
#        col = c('black','gray','black'),
#        lty = c(FALSE,FALSE,3),
#        legend = c("Start of interval", "End of interval",'Germination'),
#        cex=.75,
#        box.lty=0)
# 
# polygon(x=c(-1,-1,40,40),y=c(110,210,210,110),col = 'grey95')
# 
# segments(x0=c(0,3,0,0,15,0,0,28,0),
#          x1=c(3,3,12,15,15,24,28,28,36),
#          y0=rev(seq(120,200,by=10)),
#          lwd=1.5)
# segments(x0=c(3,15,27),
#          y0=c(200,170,140),
#          y1=c(190,160,130),
#          lwd=1,lty='dotted')
# points(x=c(3,3,12,15,15,24,28,28,36),
#        y=rev(seq(120,200,by=10)),
#        pch=c(19,21,19,19,21,19,19,21,19),
#        bg=c("black","white","black","black","white","black","black","white","black"))
# 
# mtext('Individual seed bag trials',side=2,line=0,adj=.9,col='black',cex=1)
# 
# legend(x = 27.5, y = 210,
#        pch = c(19,21,NA),
#        col = c('black','black','black'),
#        lty = c(FALSE,FALSE,3),
#        legend = c("Seed count", "Seedling count",'Germination'),
#        cex=.75,
#        box.lty=0,
#        bg=NA)

#########################
## TWO PANEL PLOT
#########################

pdf("~/Dropbox/clarkiaSeedBanks/products/miscellaneous/conceptual/persistence.pdf",width=6,height=8)

par(mfrow=c(1,1))
t.sample = c(0,4,4,12,16,16,24,28,28,36)/36

# plot exponential survival function
plot(t*36,f.exp(t=t),
     xlim=c(0,40),ylim=c(0,210),
     type='n',frame=FALSE,
     xlab = "Time (months)",
     ylab = "",
     axes=FALSE)
axis(1, c(0,10,20,30,40), col.ticks = 1)
axis(2, c(0,20,40,60,80,100,195,165,135),
     labels = c(0,20,40,60,80,100,"Age 1","Age 2","Age 3"),
     col = NA, col.ticks = c(1,1,1,1,1,1,0))
segments(x0=-1.5,y0=0,y1=100)
mtext('Seed counts',side=2,line=3,adj=.25,col='black',cex=1)

l.x.s=l.x[c(3,6,9)]*100+l.x[c(2,5,8)]*c(Jangerm_1,Jangerm_2,Jangerm_3)*100

# plot survivorship points on to curve

segments(x0=c(4,16,28),
         y0=c(l.x.s),y1=c(l.x[c(3,6,9)]*100),
         lty=c('dotted'))

points(c(t.sample[1:4],t.sample[4],t.sample[5:7],t.sample[7],t.sample[8:10])*36,
       c(l.x[1:4],l.x[4],l.x[5:7],l.x[7],l.x[8:10])*100,
       col=c('black','gray',"black","gray","black","gray","black","gray","black","gray"),
       cex=c(1,1,1,2,1,1,1,2,1,1,1,1),
       pch=c(19,19,19,19,20,19,19,19,20,19,19,19))

t1 = seq(0/36,4/36,by=.01)
lines(t1*36, f.exp(t=t1)*100)

t2 = seq(4/36,12/36,by=.01)
lines(t2*36, (1-g[1])*f.exp(t=t2)*100)

t3 = seq(12/36,16/36,by=.01)
lines(t3*36, (1-g[1])*f.exp(t=t3)*100)

t4 = seq(16/36,24/36,by=.01)
lines(t4*36, (1-g[1])*(1-g[2])*f.exp(t=t4)*100)

t5 = seq(24/36,28/36,by=.01)
lines(t5*36, (1-g[1])*(1-g[2])*f.exp(t=t5)*100)

t6 = seq(28/36,36/36,by=.01)
lines(t6*36, (1-g[1])*(1-g[2])*(1-g[3])*f.exp(t=t6)*100)

legend(x = 25, y = 100,
       pch = c(19,19,NA),
       col = c('black','gray','black'),
       lty = c(FALSE,FALSE,3),
       legend = c("Start of interval", "End of interval",'Germination'),
       cex=.75,
       box.lty=0)

polygon(x=c(-1,-1,40,40),y=c(110,210,210,110),col = 'grey95',border=FALSE)

segments(x0=c(0,4,0,0,16,0,0,28,0),
         x1=c(4,4,12,16,16,24,28,28,36),
         y0=c(200,190,200,170,160,170,140,130,140)-5,
         lwd=1.5)
segments(x0=c(0,0,0),
         y0=c(187.5,157.5,127.5)+5,
         y1=c(192.5,162.5,132.5)+5)
segments(x0=c(4,16,28),
         y0=c(200,170,140)-20,
         y1=c(190,160,130)+5,
         lwd=1,lty='dotted')
points(x=c(4,4,12,16,16,24,28,28,36),
       y=c(205,190,205,175,160,175,145,130,145)-10,
       pch=c(19,21,19,19,21,19,19,21,19),
       bg=c("black","white","black","black","white","black","black","white","black"))

text(x=4,201,expression(paste(y)[11]))
text(x=12,201,expression(paste(y)[12]))
text(x=6,180,expression(paste(y["g,1"])))

text(x=(16),171,expression(paste(y)[13]))
text(x=(24),171,expression(paste(y)[14]))
text(x=18,150,expression(paste(y["g,2"])))

text(x=(28),141,expression(paste(y)[15]))
text(x=(36),141,expression(paste(y)[16]))
text(x=30,120,expression(paste(y["g,3"])))

mtext('Seed bag trials',side=3,line=-.5,adj=.01,col='black',cex=1)

legend(x = 22.5, y = 210,
       pch = c(19,21,NA),
       col = c('black','black','black'),
       lty = c(FALSE,FALSE,3),
       legend = c("Seed count", "Seedling count",'Germination'),
       cex=.75,
       box.lty=0,
       bg=NA)

rect(c(0,4,12,16,24,28,36),-6,c(4,12,16,24,28,36,40),-2,
     col=c('gray50','gray90'),lwd=0)
text(c(0,4,12,16,24,28,36,40),1.5,
     c("O","J"),cex=.75)

dev.off()

#########################
## TWO PANEL PLOT FOR VIABILITY
#########################

pdf("~/Dropbox/clarkiaSeedBanks/products/miscellaneous/conceptual/viability.pdf",width=4,height=8)
par(mfrow=c(1,1))
t.sample = c(0,4,4,12,16,16,24,28,28,36)/36

## Panel B
## Viability and inferred viability
plot(t*36,f.exp(t=t),
     xlim=c(0,40),ylim=c(0,2.1),
     type='n',frame=FALSE,
     xlab = "Time (months)",
     ylab = "",
     axes=FALSE)

axis(1, c(0,10,20,30,40), col.ticks = 1)
axis(2, c(0,.2,.4,.6,.8,1),
     labels = c(0,.2,.4,.6,.8,1),
     col = NA, col.ticks = 1)
segments(x0=-1.5,y0=0,y1=1)
mtext('Viability probability',side=2,line=3,adj=.15,col='black',cex=1)

segments(x0=12,y0=v[1],x1=5,y1=1.3,lty='dotted')
segments(x0=12,y0=v[1],x1=35,y1=1.3,lty='dotted')


x0=seq(0,1,by=1/3)
lines(x=x0*12,v[1]^(x0),
      lty='dotdash',col='orange')

x1=seq(1,2,by=.01)
lines(x=x1*12,v[1]*(v[2]/v[1])^(x1-1),
      lty='dotdash',col='orange')

x2=seq(2,3,by=.01)
lines(x=x2*12,v[2]*(v[3]/v[2])^(x2-2),
      lty='dotdash',col='orange')

points(c(0,12,24,36),
       c(1,v),
       pch=19,
       col='orange')
points(c(4,16,28), 
       c(v[1]^(1/3), v[1]*(v[2]/v[1])^(1/3), v[2]*(v[3]/v[2])^(1/3)),
       pch=21,bg='white',
       col='orange')

legend(x = 20 , y = 1.05,
       pch = c(19,21),
       col = c('orange','orange'),
       bg = c(NA,'white'),
       legend = c("Estimated viability", "Inferred viability"),
       cex=.75,
       box.lty=0)


polygon(x=c(5,5,35,35),y=c(1.3,2,2,1.3),col = 'grey95',border=FALSE)


segments(x0=c(10),
         x1=c(30),
         y0=c(1.8))
points(x=c(10,30),y=c(1.8,1.8),pch=19)

segments(x0=c(10),
         x1=c(30),
         y0=c(1.4))
points(x=c(10,30),y=c(1.4,1.4),pch=19)

segments(x0=30,x1=11,y0=1.8,y1=1.425,lty='dotted')
shape::Arrows(x0=30,x1=11,y0=1.8,y1=1.425, arr.type="triangle", 
              arr.width=.2,arr.length=.2,lty=0)

text(x=9.5,1.9,expression(paste(n)["g,1"]^"viab"))
text(x=30,1.9,expression(paste(y)["g,1"]^"viab"))

text(x=9.5,1.5,expression(paste(n)["v,1"]^"viab"))
text(x=30,1.5,expression(paste(y)["v,1"]^"viab"))

mtext('Viability trials',side=3,line=-1.75,adj=.2,col='black',cex=1)

legend(x = 27.5, y = 210,
       pch = c(19,21,NA),
       col = c('black','black','black'),
       lty = c(FALSE,FALSE,3),
       legend = c("Seed count", "Seedling count",'Germination'),
       cex=.75,
       box.lty=0,
       bg=NA)

rect(c(0,4,12,16,24,28,36),-.05,c(4,12,16,24,28,36,40),-.01,
     col=c('gray50','gray90'),lwd=0)
text(c(0,4,12,16,24,28,36,40),.025,
     c("O","J"),cex=.75)

dev.off()
# 

#########################
## GERMINATION PROBABILITIES
#########################

pdf("~/Dropbox/clarkiaSeedBanks/products/miscellaneous/conceptual/germination.pdf",width=4,height=4)

par(mfrow=c(1,1))

## Panel C
## Conditional and unconditional germination
plot(c(1,2,3),g,
     ylim=c(0,1),
     pch=16, type = 'n', frame = FALSE,
     xlab = "Seed age (years)",
     ylab = "Germination probability")

# unconditional germination probabilities in black
points( x = c(1,2,3),
        g.uncond,
        col='gray', pch = 19)

# conditional germination probabilities in gray
points( x = c(1,2,3),
        y = g,
        pch = 19, col ='black')

col.orange=as.vector((col2rgb("orange",alpha=TRUE)/255)[,1])
col.orange[4]=.75

# germination probabilities conditional on viability in orange
points( x = c(1,2,3),
        g.cond1,
        col=rgb(1,.647,0,.75), pch = 19)

legend(x = 1, y = 1,
       pch = c(19,19,19),
       col = c(rgb(1,.647,0,.75),c('black','gray')),
       legend = c("Conditional on being viable","Conditional on being intact","Unconditional"),
       cex=.75,
       box.lty=0)

dev.off()

#########################
## SURVIVAL PROBABILITIES
#########################


pdf("~/Dropbox/clarkiaSeedBanks/products/miscellaneous/conceptual/survival.pdf",width=4,height=4)

par(mfrow=c(1,1))
t.sample = c(0,4,4,12,16,16,24,28,28,36)/36

# plot exponential survival function
plot(t*36,f.exp(t=t),
     xlim=c(0,40),ylim=c(0,1),
     type='n',frame=FALSE,
     xlab = "Time (months)",
     ylab = "",
     axes=FALSE)
axis(1, c(0,10,20,30,40), col.ticks = 1)
axis(2, c(0,.2,.4,.6,.8,1),
     labels = c(0,.2,.4,.6,.8,1),
     col = NA, col.ticks = 1)
segments(x0=-1.5,y0=0,y1=1)
mtext('Survival probability',side=2,line=3,adj=.5,col='black',cex=1)

l.x.s=l.x[c(3,6,9)]+l.x[c(2,5,8)]*c(Jangerm_1,Jangerm_2,Jangerm_3)

# plot survivorship points on to curve

segments(x0=c(4,16,28),
         y0=c(l.x.s),y1=c(l.x[c(3,6,9)]),
         lty=c('dotted'), col = 'black')

# points(c(t.sample[1:4],t.sample[4],t.sample[5:7],t.sample[7],t.sample[8:10])*36,
#        c(l.x[1:4],l.x[4],l.x[5:7],l.x[7],l.x[8:10]),
#        col=c('black','black',"black","black","black","black","black","black","black","black"),
#        cex=c(1,1,1,2,1,1,1,2,1,1,1,1),
#        pch=c(19,19,19,19,20,19,19,19,20,19,19,19))

## PERSISTENCE CURVE
t1 = seq(0/36,4/36,by=.01)
lines(t1*36, f.exp(t=t1), col = 'black')

t2 = seq(4/36,12/36,by=.01)
lines(t2*36, (1-g[1])*f.exp(t=t2), col = 'black')

t3 = seq(12/36,16/36,by=.01)
lines(t3*36, (1-g[1])*f.exp(t=t3), col = 'black')

t4 = seq(16/36,24/36,by=.01)
lines(t4*36, (1-g[1])*(1-g[2])*f.exp(t=t4), col = 'black')

t5 = seq(24/36,28/36,by=.01)
lines(t5*36, (1-g[1])*(1-g[2])*f.exp(t=t5), col = 'black')

t6 = seq(28/36,36/36,by=.01)
lines(t6*36, (1-g[1])*(1-g[2])*(1-g[3])*f.exp(t=t6), col = 'black')

## VIABILITY CURVE
l.x.s=l.x.viability[c(2,5,8)]

points(t.sample[1:10]*36,l.x.viability[1:10],col='orange',pch=19)

# plot survivorship points on to curve
# 
# segments(x0=c(4,16,28),
#          y0=c(f.exp(t=4/36)*(g[1]+(1-g[1])*v[1]^(1/3)),
#               (1-g[1])*f.exp(t=16/36)*(g[2]+(1-g[2])*v[1]*(v[2]/v[1])^(1/3)),
#               (1-g[2])*(1-g[1])*f.exp(t=28/36)*(g[3]+(1-g[3])*v[2]*(v[3]/v[2])^(1/3))),
#          y1=c((1-g[1])*(f.exp(t=4/36))*v[1]^(1/3),
#               (1-g[2])*(1-g[1])*f.exp(t=16/36)*v[1]*(v[2]/v[1])^(1/3),
#               (1-g[3])*(1-g[2])*(1-g[1])*f.exp(t=28/36)*v[2]*(v[3]/v[2])^(1/3)),
#          lty=c('solid'), col = 'orange')


 

# 
# x0=seq(0,1/3,by=.01)
# time.varying.lx = f.exp(t=x0/3)
# lines(x=x0*12,time.varying.lx*(g[1]+(1-g[1])*v[1]^(x0)),
#       lty='solid',col='orange')
# 
# x0=seq(1/3,1,length.out=25)
# time.varying.lx = (1-g[1])*(f.exp(t=x0/3))
# lines(x=x0*12,time.varying.lx*v[1]^(x0),
#       lty='solid',col='orange')
# 
# x1=seq(1,1+1/3,length.out=25)
# time.varying.lx2 = (1-g[1])*f.exp(t=x1/3)
# lines(x=x1*12,time.varying.lx2*(g[1]+(1-g[2])*v[1]*(v[2]/v[1])^(x1-1)),
#       lty='solid',col='orange')
# 
# x1=seq(1+1/3,2,by=.01)
# time.varying.lx3 = (1-g[1])*(1-g[2])*f.exp(t=x1/3)
# 
# lines(x=x1*12,time.varying.lx3*v[1]*(v[2]/v[1])^(x1-1),
#       lty='solid',col='orange')
# 
# x1=seq(2,2+1/3,by=.01)
# time.varying.lx4 = (1-g[1])*(1-g[2])*f.exp(t=x1/3)
# 
# lines(x=x1*12,time.varying.lx4*v[2]*(v[3]/v[2])^(x1-2),
#       lty='solid',col='orange')
# 
# x1=seq(2+1/3,3,by=.01)
# time.varying.lx5 = (1-g.cond1[1])*(1-g.cond1[2])*(1-g.cond1[3])*f.exp(t=x1/3)
# 
# lines(x=x1*12,time.varying.lx5*v[2]*(v[3]/v[2])^(x1-2),
#       lty='solid',col='orange')
# 
# 



legend(x = 10, y = 1,
       col = c('black','orange'),
       lty = c(1,NA),
       pch = c(NA,19),
       legend = c("Intact only","Intact and viable"),
       cex=.75,
       box.lty=0)



dev.off()




#########################
## CUMULATIVE SEED LOSS
#########################

# 
# par(mfrow=c(1,1))
# t.sample = c(0,4,4,12,15,15,24,28,28,36)/36
# 
# # plot exponential survival function
# plot(t*36,f.exp(t=t),
#      xlim=c(0,40),ylim=c(0,210),
#      type='n',frame=FALSE,
#      xlab = "Time (months)",
#      ylab = "",
#      axes=FALSE)
# axis(1, c(0,10,20,30,40), col.ticks = 1)
# axis(2, c(0,20,40,60,80,100,195,165,135),
#      labels = c(0,20,40,60,80,100,"Age 1","Age 2","Age 3"),
#      col = NA, col.ticks = c(1,1,1,1,1,1,0))
# segments(x0=-1.5,y0=0,y1=100)
# mtext('Seed counts',side=2,line=3,adj=.25,col='black',cex=1)
# 
# l.x.s=l.x[c(3,6,9)]*100+l.x[c(2,5,8)]*c(Jangerm_1,Jangerm_2,Jangerm_3)*100
# 
# # plot survivorship points on to curve
# 
# # segments(x0=c(4,15,27),
# #          y0=c(l.x.s),y1=c(l.x[c(3,6,9)]*100),
# #          lty=c('dotted'))
# 
# # points(c(t.sample[1:4],t.sample[4],t.sample[5:7],t.sample[7],t.sample[8:10])*36,
# #        c(l.x[1:4],l.x[4],l.x[5:7],l.x[7],l.x[8:10])*100,
# #        col=c('black','gray',"black","gray","black","gray","black","gray","black","gray"),
# #        cex=c(1,1,1,2,1,1,1,2,1,1,1,1),
# #        pch=c(19,19,19,19,20,19,19,19,20,19,19,19))
# 
# t1 = seq(0/36,3/36,by=.0001)
# #lines(t1*36, f.exp(t=t1)*100)
# polygon(x=c(rev(t1*36),t1*36),
#         y=c(rep(0,length(t1)),100-f.exp(t=t1)*100),
#         col='gray',border=NA)
# polygon(x=c(3,3,36,36),
#         y=c(0,100-f.exp(t=c(3/36,3/36))*100,0),
#         col='gray',border=NA)
# 
# t2 = seq(3/36,12/36,by=.0001)
# polygon(x=c(3,3,36,36),
#         y=100-c(f.exp(t=c(3/36))*100,(1-g[1])*f.exp(t=c(3/36))*100,(1-g[1])*f.exp(t=c(3/36))*100,f.exp(t=c(3/36))*100),
#         angle=30,density=10,col='gray',border=NA)
# polygon(x=c(rev(t2*36),t2*36),
#         y=100-c((1-g[1])*f.exp(t=rep(3/36,length(t2)))*100,(1-g[1])*f.exp(t=t2)*100),
#         col='gray',border=NA)
# polygon(x=c(12,12,36,36),
#         y=100-c((1-g[1])*f.exp(t=3/36)*100,(1-g[1])*f.exp(t=12/36)*100,(1-g[1])*f.exp(t=12/36)*100,(1-g[1])*f.exp(t=3/36)*100),
#         col='gray',border=NA)
# 
# 
# t3 = seq(12/36,15/36,by=.0001)
# polygon(x=c(rev(t3*36),t3*36),
#         y=100-c((1-g[1])*f.exp(t=rep(12/36,length(t3)))*100,(1-g[1])*f.exp(t=t3)*100),
#         col='gray',border=NA)
# polygon(x=c(15,15,36,36),
#         y=100-c((1-g[1])*f.exp(t=12/36)*100,(1-g[1])*f.exp(t=15/36)*100,(1-g[1])*f.exp(t=15/36)*100,(1-g[1])*f.exp(t=12/36)*100),
#         col='gray',border=NA)
# 
# 
# t4 = seq(15/36,24/36,by=.0001)
# polygon(x=c(15,15,36,36),
#         y=100-c((1-g[1])*f.exp(t=c(15/36))*100,(1-g[1])*(1-g[2])*f.exp(t=c(15/36))*100,(1-g[1])*(1-g[2])*f.exp(t=c(15/36))*100,(1-g[1])*f.exp(t=c(15/36))*100),
#         angle=30,density=10,col='gray',border=NA)
# polygon(x=c(rev(t4*36),t4*36),
#         y=100-c((1-g[1])*(1-g[2])*f.exp(t=rep(15/36,length(t4)))*100,(1-g[1])*(1-g[2])*f.exp(t=t4)*100),
#         col='gray',border=NA)
# polygon(x=c(24,24,36,36),
#         y=100-c((1-g[1])*(1-g[2])*f.exp(t=15/36)*100,(1-g[1])*(1-g[2])*f.exp(t=24/36)*100,(1-g[1])*(1-g[2])*f.exp(t=24/36)*100,(1-g[1])*(1-g[2])*f.exp(t=15/36)*100),
#         col='gray',border=NA)
# 
# t5 = seq(24/36,27/36,by=.0001)
# polygon(x=c(rev(t5*36),t5*36),
#         y=100-c((1-g[1])*(1-g[2])*f.exp(t=rep(24/36,length(t5)))*100,(1-g[1])*(1-g[2])*f.exp(t=t5)*100),
#         col='gray',border=NA)
# polygon(x=c(28,28,36,36),
#         y=100-c((1-g[1])*(1-g[2])*f.exp(t=24/36)*100,(1-g[1])*(1-g[2])*f.exp(t=27/36)*100,(1-g[1])*(1-g[2])*f.exp(t=27/36)*100,(1-g[1])*(1-g[2])*f.exp(t=24/36)*100),
#         col='gray',border=NA)
# 
# t6 = seq(27/36,36/36,by=.0001)
# polygon(x=c(28,28,36,36),
#         y=100-c((1-g[1])*(1-g[2])*f.exp(t=c(27/36))*100,(1-g[1])*(1-g[2])*(1-g[3])*f.exp(t=c(27/36))*100,(1-g[1])*(1-g[2])*(1-g[3])*f.exp(t=c(27/36))*100,(1-g[1])*(1-g[2])*f.exp(t=c(27/36))*100),
#         angle=30,density=10,col='gray',border=NA)
# polygon(x=c(rev(t6*36),t6*36),
#         y=100-c((1-g[1])*(1-g[2])*(1-g[3])*f.exp(t=rep(27/36,length(t6)))*100,(1-g[1])*(1-g[2])*(1-g[3])*f.exp(t=t6)*100),
#         col='gray',border=NA)

