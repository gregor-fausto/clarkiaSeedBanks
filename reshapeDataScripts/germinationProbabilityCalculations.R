intact = 60
germinant = 10


prob.germination = function( n = intact, y =germinant, p.v =v){
  p = y/(n*(p.v^(1/3))+y)
  return(p)
}

par(mfrow=c(1,2))
p = seq(0,1,by=0.01)
plot(p, prob.germination(n = intact, y = germinant, p.v = p),
     xlim=c(0,1), ylim=c(0,1),
     cex = .5, pch = 16,
     xlab="Probability of viability",
     ylab="Probability of germination")
abline(h=germinant/(intact+germinant),lty='dashed')

# if all the seeds are viable, 

prob.germination2 = function( n = intact, y =germinant, p.v =v){
  g = y/(y+n)
  p = g/(1-((1-(p.v^(1/3)))*(1-g)))
  return(p)
}

p = seq(0,1,by=0.01)
plot(p, prob.germination2(n = intact, y = germinant, p.v = p),
     xlim=c(0,1), ylim=c(0,1),
     cex = .5, pch = 16,
     xlab="Probability of viability",
     ylab="Probability of germination")
abline(h=germinant/(intact+germinant),lty='dashed')


