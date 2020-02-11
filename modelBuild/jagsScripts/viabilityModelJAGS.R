
    model { 
    
    theta ~ dunif(0,1)
    #kappa ~ dpar(alpha=shape,c=scale)
    # kappa ~ dpar(1.5,1)
    
    for(j in 1:nbags){
        p[j] ~ dbeta(1, 1)
        p2[j] ~ dbeta(1,1)
    }

    for(i in 1:N){

    # v viability
    yv[i] ~ dbinom( p[bag[i]] , nv[i] )
    yv2[i] ~ dbinom( p2[bag[i]] , nv2[i] )
    
    yv.sim[i] ~ dbinom( p[bag[i]], nv[i] )
    yv2.sim[i] ~ dbinom( p2[bag[i]], nv2[i] )

    prop[i] <- yv[i]/nv[i]
    prop2[i] <- yv2[i]/nv2[i]
    
    Ds[i] <- 2*nv[i]*(prop[i])*log((prop[i]+0.00001)/p[bag[i]]) + (1-prop[i])*log((1-prop[i]+0.00001)/(1-p[bag[i]]))
    
    }
    dev.sat <- sum(Ds[])
    
    # derived quantity
    for(i in 1:nbags){
        vJoint[i] = p[i] + p2[i]*(1-p[i])
    }
    
    mean.data <- mean(yv)
    mean.sim <- mean(yv.sim)
    p.mean <- step(mean.sim - mean.data)
    
    }
    
