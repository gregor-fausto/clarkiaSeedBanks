
model { 
 ##############
    ## hyperpriors
    ##############
    
    mu.i ~ dnorm( 0, 0.0001 )
    sigma.i ~ dunif(0,100)
    tau.i <- 1/(sigma.i * sigma.i)

        mu.s ~ dnorm( 0, 0.0001 )
    sigma.s ~ dunif(0,100)
    tau.s <- 1/(sigma.s * sigma.s)
    
    ##############
    ## priors
    ##############
    
    for(j in 1:nbags){
        # lab experiments
        pg[j] ~ dbeta( 1 , 1 )
        pv[j] ~ dbeta( 1 , 1 )
    }
    
    for(j in 1:nbags){
                                        # site intercepts
        alpha.i[j] ~ dnorm(mu.i, tau.i)
        alpha.s[j] ~ dnorm(mu.s, tau.s)

    }
    
    ##############
    ## likelihoods
    ##############
    
    # germination and viability experiments
    for(i in 1:N){
        # g germination
        yg[i] ~ dbinom( pg[bag[i]] , ng[i] )
        
        # v viability
        yv[i] ~ dbinom( pv[bag[i]] , nv[i] )
        
        ygSim[i] ~ dbinom( pg[bag[i]] , ng[i] )
        yvSim[i] ~ dbinom( pv[bag[i]] , nv[i] )
        
    }

    # calculate viability    
    for(j in 1:nbags){
        viability[j] = pg[j] + pv[j]*(1-pg[j])
    }

    # seed burial experiments
    for(i in 1:N_burial){
        
                                        # s1 seed survival
        pi[i] <- ilogit(alpha.i[bag_burial[i]])
        y_total[i] ~ dbin(pi[i], n_buried[i])
        
                                        # g1 seed germination
                ps[i] <- ilogit(alpha.s[bag_burial[i]])

        y_seedlings[i] ~ dbin(ps[i]*pi[i]*(viability[bag_burial[i]]^(1/3)), n_buried[i])
        
        ySeedlingsSim[i] ~ dbinom(ps[i]*pi[i]*(viability[bag_burial[i]]^(1/3)), n_buried[i])
        yTotalSim[i] ~ dbin(pi[i], n_buried[i])
    }
    
}

