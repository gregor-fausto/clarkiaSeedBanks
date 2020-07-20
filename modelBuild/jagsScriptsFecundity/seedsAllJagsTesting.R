
model {

  # # TOTAL FRUIT EQUIVALENTS
  # for(i in 1:n_site){
  # 
  #   # theta 1
  #   mu0_tfe[i] ~  dnorm(0, 0.001)
  #   sigma0_tfe[i] ~ dunif(0,10)
  #   tau0_tfe[i] <- 1/(sigma0_tfe[i]*sigma0_tfe[i])
  # 
  #   for(j in 1:n_year){
  # 
  #     # theta 1
  #     gamma_tfe[i,j] ~ dnorm(mu0_tfe[i], tau0_tfe[i])
  #     r_tfe[i,j] ~ dgamma(.001,.001)
  # 
  #   }
  # }

  for(i in 1:n_site){

    # UNDAMAGED FRUITS
    # theta 1
    mu0_und[i] ~  dnorm(0, 1)
    sigma0_und[i] ~ dnorm(0,1) T(0,)
    tau0_und[i] <- 1/(sigma0_und[i]*sigma0_und[i])

    # DAMAGED FRUITS
    # theta 1
    # mu0_dam[i] ~  dnorm(0, 0.001)
    # sigma0_dam[i] ~ dunif(0,10)
    # tau0_dam[i] <- 1/(sigma0_dam[i]*sigma0_dam[i])

    for(j in 1:n_year2){

      # theta 1
      gamma_und[i,j] ~ dnorm(mu0_und[i], tau0_und[i])
      #r_und[i,j] ~ dgamma(.001,.001)
      r_und.inv[i,j] ~ dnorm(0,1) T(0,)
      r_und[i,j] = 1/r_und.inv[i,j]
      
      # theta 1
      # gamma_dam[i,j] ~ dnorm(mu0_dam[i], tau0_dam[i])
      # r_dam[i,j] ~ dgamma(.001,.001)

    }
  }

  # # UNDAMAGED SEEDS
  # for(i in 1:n_site){
  # 
  #   # theta 1
  #   mu0_seeds[i] ~  dnorm(0, 0.001)
  #   sigma0_seeds[i] ~ dunif(0,10)
  #   tau0_seeds[i] <- 1/(sigma0_seeds[i]*sigma0_seeds[i])
  # 
  #   for(j in 1:n_year3){
  # 
  #     # theta 1
  #     gamma_seeds[i,j] ~ dnorm(mu0_seeds[i], tau0_seeds[i])
  #     r_seeds[i,j] ~ dgamma(.001,.001)
  # 
  #   }
  # }
  # 
  # # DAMAGED SEEDS
  # for(i in 1:n_site){
  # 
  #   # theta 1
  #   mu0_dam_seeds[i] ~  dnorm(0, 0.001)
  #   sigma0_dam_seeds[i] ~ dunif(0,10)
  #   tau0_dam_seeds[i] <- 1/(sigma0_dam_seeds[i]*sigma0_dam_seeds[i])
  # 
  #   for(j in 1:n_year4){
  # 
  #     # theta 1
  #     gamma_dam_seeds[i,j] ~ dnorm(mu0_dam_seeds[i], tau0_dam_seeds[i])
  #     r_dam_seeds[i,j] ~ dgamma(.001,.001)
  # 
  #   }
  # }

  # LIKELIHOODS
  # Total fruit equivalents
  # for (i in 1:n){
  #   lambda_tfe[i] = exp(gamma_tfe[site[i],year[i]])
  #   y_tfe[i] ~ dnegbin(r_tfe[site[i],year[i]]/(r_tfe[site[i],year[i]]+lambda_tfe[i]),r_tfe[site[i],year[i]])
  # }

  # Undamaged fruits
  for (i in 1:n2){
    log(lambda_und[i]) = gamma_und[site[i],year2[i]]
    y_und[i] ~ dnegbin(r_und[site[i],year2[i]]/(r_und[site[i],year2[i]]+lambda_und[i]),r_und[site[i],year2[i]])
    y_und_sim[i] ~ dnegbin(r_und[site[i],year2[i]]/(r_und[site[i],year2[i]]+lambda_und[i]),r_und[site[i],year2[i]])
    }

  # Damaged fruits
  # for (i in 1:n2){
  #   lambda_dam[i] = exp(gamma_dam[site[i],year2[i]])
  #   y_dam[i] ~ dnegbin(r_dam[site[i],year2[i]]/(r_dam[site[i],year2[i]]+lambda_dam[i]),r_dam[site[i],year2[i]])
  # }

  # # Seeds per undamaged fruit
  # for (i in 1:n3){
  #   lambda_seeds[i] = exp(gamma_seeds[site[i],year3[i]])
  #   sdno[i] ~ dnegbin(r_seeds[site[i],year3[i]]/(r_seeds[site[i],year3[i]]+lambda_seeds[i]),r_seeds[site[i],year3[i]])
  # }
  # 
  # # Seeds per damaged fruits
  # for (i in 1:n4){
  #   lambda_dam_seeds[i] = exp(gamma_dam_seeds[site[i],year4[i]])
  #   sdno_dam[i] ~ dnegbin(r_dam_seeds[site[i],year4[i]]/(r_dam_seeds[site[i],year4[i]]+lambda_dam_seeds[i]),r_dam_seeds[site[i],year4[i]])
  # }

  # for other use
  # exp(gamma_dam[site[i],year[i]])/exp(gamma[site[i],year[i]])

  ## DERIVED QUANTITIES
  # calculate the fraction of seeds in a damaged versus undamaged fruit (only for 2013-2018)
  # for(i in 1:n_site){
  #   for(j in 1:n_year4){
  #     ratio[i,j] = exp(gamma_dam_seeds[site[i],year[i]] - gamma_seeds[site[i],year[i]])
  #   }
  # }

  # use the ratio to calculate the total fruit equivalents
  # add the ratio*number of damaged fruits to the number of undamaged fruits to get a TFE

  for(i in 1:n_site){
# 
#     # p0_inv_tfe[i] ~ dnorm(mu0_tfe[i],tau0_tfe[i])
#     # p0_tfe[i] <- exp(p0_inv_tfe[i])
# 
    p0_inv_und[i] ~ dnorm(mu0_und[i],tau0_und[i])
    p0_und[i] <- exp(p0_inv_und[i])
# 
#     p0_inv_dam[i] ~ dnorm(mu0_dam[i],tau0_dam[i])
#     p0_dam[i] <- exp(p0_inv_dam[i])
# 
#     # p0_inv_und_seeds[i] ~ dnorm(mu0_seeds[i],tau0_seeds[i])
#     # p0_seeds[i] <- exp(p0_inv_und_seeds[i])
#     # 
#     # p0_inv_dam_seeds[i] ~ dnorm(mu0_dam_seeds[i],tau0_dam_seeds[i])
#     # p0_dam_seeds[i] <- exp(p0_inv_dam_seeds[i])
# 
#     # for 2006-2012
#     for(j in 1:n_year){
#       
#       p_inv_und[i] ~ dnorm(mu0_und[i],tau0_und[i])
#       p_und[i] <- exp(p_inv_und[i])
# 
#       tfe[i,j] = exp(gamma_tfe[i,j])
# 
#     }
# 
  # for 2013-2018
  for(j in 1:n_year2){
    
      log(lambda_und[i]) = gamma_und[site[i],year2[i]]
# 
#     p_und[i,j] = exp(gamma_und[i,j])
#     p_dam[i,j] = exp(gamma_dam[i,j])
# 
#     tfe_comp[i,j] = p_und[i,j] + p_dam[i,j]*ratio[i,j]
# 
  }
  }
# }

#   for (i in 1:n){
#     lambda[i] = exp(gamma[site[i],year[i]])
#
#
#     ifelse(year[i]<2013,
#      countFruitsPerPlant[i] ~ dnegbin(r[site[i],year[i]]/(r[site[i],year[i]]+lambda[i]),r[site[i],year[i]]),
#
#
#     # for 2006-2012, fruits per plant can be estimated directly from the data
#     for (year[i] in 2006:2012) {
#     }
#     # for 2013-2018, fruits per plant needs to be transformed into total fruit equivalents before fitting to the data
#     for (year[i] in 2013:2018) {
#       # calculate the fraction of seeds in a damaged versus undamaged fruit
#       ratio[site[i],year[i]] = exp(gamma_dam[site[i],year[i]] - gamma[site[i],year[i]])
#       countTFE[i] = n_undamaged[i] + round( ratio[j,k] * n_damaged[i] )
#
#       countTFE[i] ~ dnegbin(r[site[i],year[i]]/(r[site[i],year[i]]+lambda[i]),r[site[i],year[i]])
#
#     }
#
#
#     lambda[i] = exp(gamma[site[i],year[i]])
#
#     ratio[site[i],year[i]] exp(gamma_dam[site[i],year[i]] - gamma[site[i],year[i]])
#
#     n_undamaged + round( ratio[j,k] * n_damaged ) ~ dnegbin(r[site[i],year[i]]/(r[site[i],year[i]]+lambda[i]),r[site[i],year[i]])
#     countFruitsPerPlant[i] ~ dnegbin(r[site[i],year[i]]/(r[site[i],year[i]]+lambda[i]),r[site[i],year[i]])
#
#
# }
 }
