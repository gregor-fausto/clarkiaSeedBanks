
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

    # # seed bags
    # for(j in 1:nsites){
    # sigmaS1[j] ~ dunif(0,100)
    # tauS1[j] <- 1/(sigmaS1[j] * sigmaS1[j])
    # 
    # sigmaG1[j] ~ dunif(0,100)
    # tauG1[j] <- 1/(sigmaG1[j] * sigmaG1[j])
    # 
    # sigmaS2[j] ~ dunif(0,100)
    # tauS2[j] <- 1/(sigmaS2[j] * sigmaS2[j])
    # }

    ##############
    ## priors
    ##############
    
    # site intercepts
    for(j in 1:nsites){
    alphaS1[j] ~ dnorm(0, .001)
    alphaG1[j] ~ dnorm(0, .001)
    alphaS2[j] ~ dnorm(0, .001)
    alphaS3[j] ~ dnorm(0, .001)
    }

    # year intercepts
    # for(j in 1:nsites){
    # for(k in 1:nyears){
    # betaS1[j,k] ~ dnorm(0, tauS1[j])
    # betaG1[j,k] ~ dnorm(0, tauG1[j])
    # betaS2[j,k] ~ dnorm(0, tauS2[j])
    # }
    # }
    # 
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
    # yv.sim[i] ~ dbinom(p[i], nv[i]) 

    # s1 seed survival
    ps[i] = ilogit(alphaS1[site[i]])
    yt[i] ~ dbin(ps[i], n[i])
    # yt.sim[i] ~ dbinom(ps[i], n[i]) 

    # g1 seed germination
    pg[i] = ilogit(alphaG1[site[i]])
    yg[i] ~ dbin(pg[i]*(p[i])^(1/3), yt[i])
    # yg.sim[i] ~ dbinom(pg[i]*(p[i])^(1/3), yt[i]) 

    # s2 seed survival
    pr[i] = ilogit(alphaS2[site[i]])
    yo[i] ~ dbin(pr[i], yt[i]-yg[i])
    # yo.sim[i] ~ dbinom(pr[i], yt[i]-yg[i]) 

    } 

    for(i in 1:N2){
    
    ps2[i] = ilogit(alphaS1[site[i]])
    pg2[i] = ilogit(alphaG1[site[i]])
    pr2[i] = ilogit(alphaS2[site[i]])

    # s3 seed survival
    ps3[i] = ilogit(alphaS3[site2[i]])
    yt2[i] ~ dbin(ps2[i]*(1-pg2[i])*pr2[i]*ps3[i], n2[i])
    # yt2.sim[i] ~ dbinom(ps2[i]*(1-pg2[i])*pr2[i]*ps3[i], n2[i])
    
    }
    }
    
