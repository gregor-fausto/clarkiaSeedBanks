
    model { 

    ##############
    ## priors
    ##############

        # lab experiments
        pg ~ dunif( 0 , 1 )
        pv ~ dunif( 0 , 1 )

    # site intercepts
    pi ~ dunif( 0 , 1 )
    ps ~ dunif( 0 , 1 )

    ##############
    ## likelihoods
    ##############

        for(i in 1:N){
    # g germination
    yg[i] ~ dbinom( pg , ng[i] )

    # v viability
    yv[i] ~ dbinom( pv , nv[i] )

    ygSim[i] ~ dbinom( pg , ng[i] )
    yvSim[i] ~ dbinom( pv , nv[i] )

        }

    # derived quantities block

    viability = pg + pv*(1-pg)

        for(i in 1:N_burial){

            # s1 seed survival
            y_total[i] ~ dbin(pi, n_buried[i])

            # g1 seed germination
            y_seedlings[i] ~ dbin(ps*(viability^(1/3)), y_total[i])
            
            ySeedlingsSim[i] ~ dbinom(ps*(viability^(1/3)), y_total[i]) 
            yTotalSim[i] ~ dbin(pi, n_buried[i])
    } 

    }
    
