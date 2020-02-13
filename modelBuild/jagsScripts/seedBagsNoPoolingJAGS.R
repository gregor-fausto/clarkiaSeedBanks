
    model { 

    ##############
    ## priors
    ##############

        for(j in 1:nbags_burial){
        # lab experiments
        pg[j] ~ dbeta( 1 , 1 )
        pv[j] ~ dbeta( 1 , 1 )
        }

        for(i in 1:N_burial){
            # site intercepts
            pi[i] ~ dunif( 0 , 1 )
            ps[i] ~ dunif( 0 , 1 )
            }

    ##############
    ## likelihoods
    ##############

        for(i in 1:N){
    # g germination
    yg[i] ~ dbinom( pg[bag[i]] , ng[i] )

    # v viability
    yv[i] ~ dbinom( pv[i] , nv[i] )

    ygSim[i] ~ dbinom( pg[i] , ng[i] )
    yvSim[i] ~ dbinom( pv[i] , nv[i] )

    viability[i] = pg[i] + pv[i]*(1-pg[i])

        }

        for(i in 1:N_burial){

            # s1 seed survival
            y_total[i] ~ dbin(pi[i], n_buried[i])

            # g1 seed germination
            y_seedlings[i] ~ dbin(ps[i]*(viability[bag_burial[i]]^(1/3)), y_total[i])
            
            ySeedlingsSim[i] ~ dbinom(ps[i]*(viability[bag_burial[i]]^(1/3)), y_total[i]) 
            yTotalSim[i] ~ dbin(pi[i], n_buried[i])
    } 

    }
    
