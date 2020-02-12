
    model { 
    
    # theta ~ dunif(0,1)
    # kappa ~ dpar(alpha=shape,c=scale)
    # kappa ~ dpar(1.5,1)
    
    for(j in 1:nbags){
        p[j] ~ dbeta(1, 1)
        p2[j] ~ dbeta(1, 1)
    }

    for(i in 1:N){

    # v viability
    yv[i] ~ dbinom( p[bag[i]] , nv[i] )
    yv2[i] ~ dbinom( p2[bag[i]] , nv2[i] )
    
    yv.sim[i] ~ dbinom( p[bag[i]], nv[i] )
    yv2.sim[i] ~ dbinom( p2[bag[i]], nv2[i] )

    # # code for deviance from Lunn 2013
    prop[i] <- yv[i]/nv[i]
    prop2[i] <- yv2[i]/nv2[i]
    
    Ds[i] <- 2*nv[i]*(prop[i])*log((prop[i]+0.00001)/p[bag[i]]) + (1-prop[i])*log((1-prop[i]+0.00001)/(1-p[bag[i]]))
    Ds2[i] <- 2*nv2[i]*(prop2[i])*log((prop2[i]+0.00001)/p2[bag[i]]) + (1-prop2[i])*log((1-prop2[i]+0.00001)/(1-p2[bag[i]]))

   # sign[i] <- 2*step(prop[i] - p[bag[i]]) - 1
   # dev.res[i] <- sign[i]*sqrt(Ds[i])
    
    }
    
    # calculate saturated deviance
    dev.sat <- sum(Ds[])
    dev.sat2 <- sum(Ds2[])
    
    # derived quantity
    for(i in 1:nbags){
        vJoint[i] = p[i] + p2[i]*(1-p[i])
    }
    
    mean.data <- mean(yv)
    mean.sim <- mean(yv.sim)
    p.mean <- step(mean.sim - mean.data)
    
    }
    
