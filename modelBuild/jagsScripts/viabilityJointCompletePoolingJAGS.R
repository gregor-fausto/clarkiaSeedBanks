
    model { 

    pg ~ dunif(0, 1) 
    pv ~ dunif(0, 1)

    for(i in 1:N){

    # g germination
    yg[i] ~ dbinom( pg , ng[i] )

    # v viability
    yv[i] ~ dbinom( pv , nv[i] )

    yg.sim[i] ~ dbinom( pg , ng[i] )
    yv.sim[i] ~ dbinom( pv , nv[i] )

    }

    # derived quantities block

    viability = pg + pv*(1-pg)
    
    }
