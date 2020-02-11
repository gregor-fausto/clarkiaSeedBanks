
    model { 
    
    ##############
    ## hyperpriors
    ##############

    # viability trials
    for(j in 1:nsites){
    for(k in 1:nyears){ 
    mu.b[j,k] ~ dnorm(0,.0001)
    
    sigma.b[j,k] ~ dunif(0,100)
    tau.b[j,k] <- 1/(sigma.b[j,k] * sigma.b[j,k])
    }
    }


    ##############
    ## priors
    ##############
    
    # viability trials 1 prior for each trial
    for(i in 1:N){
    alpha[i] ~ dnorm(mu.b[site[i],year[i]],  tau.b[site[i],year[i]])
    }

    ##############
    ## likelihoods
    ##############

    for(i in 1:N){

    # v viability
    p[i] <- ilogit(alpha[i])
    yv[i] ~ dbin( p[i] , nv[i])
     yv.sim[i] ~ dbinom(p[i], nv[i]) 
    
    }
    }
    
