
load(file="/Users/Gregor/Dropbox/clarkiaSeedBanks/modelBuild/output/seedbagfit.rds")

getParameters <- function(codaObject,x) { codaObject[,stringr::str_detect(colnames(codaObject),x)] }
  
aS1 <- getParameters(zc[[1]],"alphaS1")
ps<-apply(aS1,2,boot::inv.logit)

for(i in 1:length(data$n)){
rbinom(n=1,size=data$n[1],p=ps[,1])
rbinom(data$n
}

# yv.sim[i] ~ dbinom(p[i], nv[i]) 


ps[i] = ilogit(alphaS1[site[i]])

# yt.sim[i] ~ dbinom(ps[i], n[i]) 
# yg.sim[i] ~ dbinom(pg[i]*(p[i])^(1/3), yt[i]) 
# yo.sim[i] ~ dbinom(pr[i], yt[i]-yg[i]) 


# yt2.sim[i] ~ dbinom(ps2[i]*(1-pg2[i])*pr2[i]*ps3[i], n2[i])

load(file="/Users/Gregor/Dropbox/modelsF2019/output/seedbagfit")
zc[[1]][,stringr::str_detect(colnames(zc[[1]]),pattern="yt.sim")]
