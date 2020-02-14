
model { 
    
    ##############
    ## priors
    ##############
    
    for(j in 1:nbags){
        # lab experiments
        pg[j] ~ dbeta( 1 , 1 )
        pv[j] ~ dbeta( 1 , 1 )
    }
    
    for(j in 1:nsites){
        # site intercepts
        pi[j] ~ dbeta( 1 , 1 )
        ps[j] ~ dbeta( 1 , 1 )
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
        y_total[i] ~ dbin(pi[site[i]], n_buried[i])
        
        # g1 seed germination
        y_seedlings[i] ~ dbin(ps[site[i]]*pi[site[i]]*(viability[bag_burial[i]]^(1/3)), n_buried[i])
        
        ySeedlingsSim[i] ~ dbinom(ps[site[i]]*pi[site[i]]*(viability[bag_burial[i]]^(1/3)), n_buried[i])
        yTotalSim[i] ~ dbin(pi[site[i]], n_buried[i])
    }
    
}

