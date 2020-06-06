model {
    
    mu0_f[i] ~  dnorm(0, 0.001)
    sigma0_f[i] ~ dunif(0,1.5)
    tau0_f[i] <- 1/(sigma0_f[i]*sigma0_f[i])        

}
