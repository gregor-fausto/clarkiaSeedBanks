
    model { 
    
    for(j in 1:nbags){
        p[j] ~ dbeta(1, 1)
    }

    for(i in 1:N){

    # v viability
    yv[i] ~ dbinom( p[bag[i]] , nv[i])
    yv.sim[i] ~ dbinom(p[bag[i]], nv[i]) 
    
    }
    }
    
