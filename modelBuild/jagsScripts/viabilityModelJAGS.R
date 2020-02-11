
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

    }
    
    # derived quantity
    for(i in 1:nbags){
        vJoint[i] = p[i] + p2[i]*(1-p[i])
    }
    
    }
    
