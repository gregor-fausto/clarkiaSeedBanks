
model { 
    ##############
    ## hyperpriors
    ##############
    
    for(i in 1:nsites){
    mu.i[i] ~ dnorm( 0, 0.0001 )
    sigma.i[i] ~ dunif(0,100)
    tau.i[i] <- 1/(sigma.i[i] * sigma.i[i])

    mu.s[i] ~ dnorm( 0, 0.0001 )
    sigma.s[i] ~ dunif(0,100)
    tau.s[i] <- 1/(sigma.s[i] * sigma.s[i])
    }

    for(i in 1:nsiteyears){
    mu.b.i[i] ~ dnorm( 0, 0.0001 )
    sigma.b.i[i] ~ dunif(0,100)
    tau.b.i[i] <- 1/(sigma.b.i[i] * sigma.b.i[i])

    mu.b.s[i] ~ dnorm( 0, 0.0001 )
    sigma.b.s[i] ~ dunif(0,100)
    tau.b.s[i] <- 1/(sigma.b.s[i] * sigma.b.s[i])
    }
    
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
        alpha.i[j] ~ dnorm(mu.i[site[j]], tau.i[site[j]])
        alpha.s[j] ~ dnorm(mu.s[site[j]], tau.s[site[j]])
	beta.i[j] ~ dnorm(mu.b.i[siteyear[j]], tau.b.i[siteyear[j]])
	beta.s[j] ~ dnorm(mu.b.s[siteyear[j]], tau.b.s[siteyear[j]])
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
                pi[i] <- ilogit(alpha.i[bag_burial[i]]+beta.i[bag_burial[i]])
        y_total[i] ~ dbin(pi[i], n_buried[i])
        
                                        # g1 seed germination
                        ps[i] <- ilogit(alpha.s[bag_burial[i]]+beta.s[bag_burial[i]])

        y_seedlings[i] ~ dbin(ps[i]*pi[i]*(viability[bag_burial[i]]^(1/3)), n_buried[i])
        
        ySeedlingsSim[i] ~ dbinom(ps[i]*pi[i]*(viability[bag_burial[i]]^(1/3)), n_buried[i])
        yTotalSim[i] ~ dbin(pi[i], n_buried[i])
    }
    
}

