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


# ######################### 
# ## Figure version 1
# #########################


# two one panel plots
pdf("~/Dropbox/clarkiaSeedBanks/products/manuscript/figures-overview/seed-bag-trials.pdf",width=8,height=4)

par(mfrow=c(1,1),mar=c(0,.25,.5,0),
    oma=c(4,4,1,1))
t.sample = c(0,4,4,12,16,16,24,28,28,36)/36

# plot exponential survival function
plot(t*36,f.exp(t=t),
     xlim=c(0,40),ylim=c(0,40),
     type='n',frame=FALSE,
     xlab = "Time (months)",
     ylab = "",
     axes=FALSE)
axis(1, c(0,10,20,30,40), col.ticks = 1,cex.axis = 1.25)
axis(2, c(35,22.5,10),
     labels = c("Age 1","Age 2","Age 3"),
     col = "black", col.ticks = c(1,1,1,1,1,1,0), las =2 ,cex.axis=1.25)

arrows(0, 35, 12, 35, length=0.05, angle=90, code=1,lwd=2)
arrows(0, 22.5, 24, 22.5, length=0.05, angle=90, code=1,lwd=2)
arrows(0, 10, 36, 10, length=0.05, angle=90, code=1,lwd=2)

segments(x0=c(4,16,28),
         y0=c(35,22.5,10),
         y1=c(31,18.5,6),
         lwd=2,lty='dotted')
points(x=c(4,4,12,16,16,24,28,28,36),
       y=c(35,31,35,22.5,18.5,22.5,10,6,10),
       pch=c(19,21,19,19,21,19,19,21,19),cex=2,
       bg=c("#0077bb","white","#0077bb","#0077bb","white","#0077bb","#0077bb","white","#0077bb"),
       col=c("#0077bb","#e7298a","#0077bb","#0077bb","#e7298a","#0077bb","#0077bb","#e7298a","#0077bb"))

text(x=4,38,expression(paste(y)[11]),cex=1.5)
text(x=12,38,expression(paste(y)[12]),cex=1.5)
text(x=6,31,expression(paste(y["g,1"])),cex=1.5)

text(x=(16),25.5,expression(paste(y)[13]),cex=1.5)
text(x=(24),25.5,expression(paste(y)[14]),cex=1.5)
text(x=18,18.5,expression(paste(y["g,2"])),cex=1.5)

text(x=(28),13,expression(paste(y)[15]),cex=1.5)
text(x=(36),13,expression(paste(y)[16]),cex=1.5)
text(x=30,6,expression(paste(y["g,3"])),cex=1.5)


legend(x = 25, y = 36,
       pch = c(19,21),
       col = c('#0077bb','#e7298a'),#,'black','black'),
       lty = c(FALSE,FALSE,1,3),
       legend = c("Seed count", "Seedling count"),#,"Survival",'Germination'),
       cex=1.5,
       box.lty=0,
       bg=NA)

rect(c(0,4,12,16,24,28,36),-1,c(4,12,16,24,28,36,40),0,
     col=c('gray50','gray90'),border=0,lwd=0)
text(c(0,4,12,16,24,28,36,40),2,
     c("O","J"),cex=1.25)
mtext("Time (months)",side=1,line=2.1,cex=1.75)

dev.off()


# two one panel plots
pdf("~/Dropbox/clarkiaSeedBanks/products/manuscript/figures-overview/lab-trials.pdf",width=4,height=3)

par(mfrow=c(1,1),mar=c(0,.25,.5,0),
    oma=c(4,4,1,1))
t.sample = c(0,4,4,12,16,16,24,28,28,36)/36

## Panel B
## Viability and inferred viability
plot(t*36,f.exp(t=t),
     xlim=c(9,21),ylim=c(.4,.7),
     type='n',frame=FALSE,
     xlab = "Time (months)",
     ylab = "",
     axes=FALSE)

segments(x0=c(10),
         x1=c(20),
         y0=c(.6))
points(x=c(10,20),y=c(.6,.6),pch=19)

segments(x0=c(10),
         x1=c(20),
         y0=c(.4))
points(x=c(10,20),y=c(.4,.4),pch=19)

segments(x0=20,x1=10.25,y0=.6,y1=.41,lty='dotted')
shape::Arrows(x0=20,x1=10.25,y0=.6,y1=.41, arr.type="triangle", 
              arr.width=.2,arr.length=.2,lty=0)

text(x=9.5,.625,expression(paste(n)["g,1"]^"viab"))
text(x=20.5,.625,expression(paste(y)["g,1"]^"viab"))

text(x=9.5,.425,expression(paste(n)["v,1"]^"viab"))
text(x=20.5,.425,expression(paste(y)["v,1"]^"viab"))

dev.off()

pdf("~/Dropbox/clarkiaSeedBanks/products/manuscript/figures-overview/empty.pdf",width=6,height=6)
plot(NA,NA)
dev.off()

pdf("~/Dropbox/clarkiaSeedBanks/products/manuscript/figures-overview/viability-data.pdf",width=4,height=4)
par(mfrow=c(1,1),mar=c(0,.25,.5,0),
    oma=c(4,4,1,1))
t.sample = c(0,4,4,12,16,16,24,28,28,36)/36

## Panel B
## Viability and inferred viability
plot(t*36,f.exp(t=t),
     xlim=c(0,40),ylim=c(.5,1),
     type='n',frame=FALSE,
     xlab = "Time (months)",
     ylab = "",
     axes=FALSE)

axis(1, c(0,10,20,30,40), col.ticks = 1)
axis(2, seq(.5,1,by=.1),
     labels = seq(.5,1,by=.1),
     col = 'black', col.ticks = 1)
mtext('Probability seed is viable',side=2,line=3,adj=.5,col='black',cex=1)

x0=seq(0,1,by=1/3)
lines(x=x0*12,v[1]^(x0),
      lty='dotdash',col='#1b9e77')

x1=seq(1,2,by=.01)
lines(x=x1*12,v[1]*(v[2]/v[1])^(x1-1),
      lty='dotdash',col='#1b9e77')

x2=seq(2,3,by=.01)
lines(x=x2*12,v[2]*(v[3]/v[2])^(x2-2),
      lty='dotdash',col='#1b9e77')

points(c(0,12,24,36),
       c(1,v),
       pch=19,
       col='#1b9e77')
points(c(4,16,28), 
       c(v[1]^(1/3), v[1]*(v[2]/v[1])^(1/3), v[2]*(v[3]/v[2])^(1/3)),
       pch=21,bg='white',
       col='#1b9e77')

legend(x = 22.5 , y = 1,
       pch = c(19,21),
       col = c('#1b9e77','#1b9e77'),
       bg = c(NA,'white'),
       legend = c("Estimated viability", "Inferred viability"),
       cex=.75,
       box.lty=0)

rect(c(0,4,12,16,24,28,36),.5,c(4,12,16,24,28,36,40),.525,
     col=c('gray50','gray90'),lwd=0,border=0)
text(c(0,4,12,16,24,28,36,40),.54,
     c("O","J"),cex=.75)

dev.off()

#1b9e77 purple
#e7298a green
pdf("~/Dropbox/clarkiaSeedBanks/products/manuscript/figures-overview/survival.pdf",width=8,height=8)

par(mfrow=c(1,1),mar=c(0,.25,.5,0),
    oma=c(4,4,1,1))
t.sample = c(0,4,4,12,16,16,24,28,28,36)/36

# plot exponential survival function
plot(t*36,f.exp(t=t),
     xlim=c(0,40),ylim=c(0,1),
     type='n',frame=FALSE,
     xlab = "Time (months)",
     ylab = "",
     axes=FALSE)
axis(1, c(0,10,20,30,40), col.ticks = 1,cex.axis=1.5,las=1)
axis(2, c(0,.2,.4,.6,.8,1),
     labels = c(0,.2,.4,.6,.8,1),
     col = 'black', col.ticks = 1,cex.axis=1.5,las=2)
mtext('Probability',side=2,line=3,adj=.5,col='black',cex=1.5)

l.x.s=l.x[c(3,6,9)]+l.x[c(2,5,8)]*c(Jangerm_1,Jangerm_2,Jangerm_3)

# plot survivorship points on to curve

segments(x0=c(4,16,28),
         y0=c(l.x.s),y1=c(l.x[c(3,6,9)]),
         lty=c('dotted'), col = '#e7298a',lwd=2)

# points(c(t.sample[1:4],t.sample[4],t.sample[5:7],t.sample[7],t.sample[8:10])*36,
#        c(l.x[1:4],l.x[4],l.x[5:7],l.x[7],l.x[8:10]),
#        col=c('black','black',"black","black","black","black","black","black","black","black"),
#        cex=c(1,1,1,2,1,1,1,2,1,1,1,1),
#        pch=c(19,19,19,19,20,19,19,19,20,19,19,19))

## PERSISTENCE CURVE
t1 = seq(0/36,4/36,by=.01)
lines(t1*36, f.exp(t=t1), col = '#0077bb',lwd=2)

t2 = seq(4/36,12/36,by=.01)
lines(t2*36, (1-g[1])*f.exp(t=t2), col = '#0077bb',lwd=2)

t3 = seq(12/36,16/36,by=.01)
lines(t3*36, (1-g[1])*f.exp(t=t3), col = '#0077bb',lwd=2)

t4 = seq(16/36,24/36,by=.01)
lines(t4*36, (1-g[1])*(1-g[2])*f.exp(t=t4), col = '#0077bb',lwd=2)

t5 = seq(24/36,28/36,by=.01)
lines(t5*36, (1-g[1])*(1-g[2])*f.exp(t=t5), col = '#0077bb',lwd=2)

t6 = seq(28/36,36/36,by=.01)
lines(t6*36, (1-g[1])*(1-g[2])*(1-g[3])*f.exp(t=t6), col = '#0077bb',lwd=2)

## VIABILITY CURVE
l.x.s=l.x.viability[c(2,5,8)]

points(t.sample[1:10]*36,l.x.viability[1:10],col='#1b9e77',pch=19,cex=2)


legend(x = 17.5, y = 1,
       col = c('#0077bb','#1b9e77'),
       lty = c(1,NA),
       pch = c(NA,19),
       legend = c("Estimated from seed bag trial","Adjusted with viability assays"),
       cex=1.1,
       box.lty=0)

mtext("Time (months)",side=1,line=2.1,cex=1.5)


dev.off()



