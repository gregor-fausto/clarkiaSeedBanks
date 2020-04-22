
### Model development

Here, I document progress developing models for the Clarkia demography
project.

  - [Table of datasets](#table-of-datasets)
  - [Belowground vital rates](#belowground-vital-rates)
  - [Seedling survival to fruiting](#seedling-survival-to-fruiting)
  - [Fruits per plant](#fruits-per-plant)
  - [Seeds per fruit](#seeds-per-fruit)
  - [Transition to seed bank](#transition-to-seed-bank)

I used pandoc to convert the latex file to a md file:

`pandoc -s parameter-table.tex -o parameter-table.md`

### Packages

I am using the packages `tidyverse` and `readxl` (documentation:
<https://readxl.tidyverse.org/>).

Load the libraries for data processing (see
<https://github.com/r-lib/rlang/issues/669> for the overwrite message I
am suppressing).

``` r
library(tidyverse)
library(knitr)
library(ggplot2)
```

For model building I am using:

``` r
library(rjags)
library(tidybayes)
library(MCMCvis)
```

### Table of datasets

|                                                  |                  |                          |           |
| :----------------------------------------------- | :--------------- | :----------------------: | :-------: |
| <span class="smallcaps">Seed vital rates</span>  | —                |            —             |     —     |
| Seed survival and germination                    | Seed bag burial  |  \(\bm{\mathrm{Y}}_1\)   | 2006-2009 |
| Seed viability                                   | Viability trials |  \(\bm{\mathrm{Y}}_2\)   | 2006-2009 |
| Seed survival and germination                    | Seed pots        |  \(\bm{\mathrm{Y}}_3\)   | 2013-2019 |
| <span class="smallcaps">Seedling survival</span> | —                |            —             |     —     |
| Seedling survival to fruiting                    | Field surveys    |  \(\bm{\mathrm{Y}}_4\)   | 2006-2019 |
| <span class="smallcaps">Fruits per plant</span>  | —                |            —             |     —     |
| Total fruit equivalents per plant                | Field surveys    |  \(\bm{\mathrm{Y}}_5\)   | 2006-2012 |
| Undamaged and damaged fruits per plant           | Field surveys    |  \(\bm{\mathrm{Y}}_6\)   | 2013-2019 |
| Total fruit equivalents per plant                | Extra plots      |  \(\bm{\mathrm{Y}}_7\)   | 2006-2012 |
| Undamaged and damaged fruits per plant           | Extra plots      |  \(\bm{\mathrm{Y}}_8\)   | 2013-2019 |
| <span class="smallcaps">Seeds per fruit</span>   | —                |            —             |     —     |
| Seeds per undamaged fruit                        | Lab counts       |  \(\bm{\mathrm{Y}}_9\)   | 2006-2019 |
| Seeds per damaged fruit                          | Lab counts       | \(\bm{\mathrm{Y}}_{10}\) | 2013-2019 |

### Belowground vital rates

Use \(Y_1\) and \(Y_2\) to estimate belowground vital rates. Refer to
the [appendix on using conditional
probability](https://github.com/gregor-fausto/clarkiaSeedBanks/blob/master/products/manuscript/appendix-x-conditional-probability.pdf).

Basic process. Use the data from seed bag burial experiment, germination
and viability trials to calculate \(\theta\). Each experiment is a trial
with a binomial likelihood, priors and hyperpriors that are centered and
logit transformed. I then take these probabilities use them to derive
the logit-transformed probabilities, which I then transform into the
probabilities of viability, survival, and germination.

    model { 
    
    ### Hyperpriors for germination and viability trials
    for(k in 1:n_siteViab){
    
      # germination
      mu0_g[k] ~  dnorm(0, 0.001)
      sigma0_g[k] ~ dunif(0,1.5)
      tau0_g[k] <- 1/(sigma0_g[k]*sigma0_g[k])
    
      # viability
      mu0_v[k] ~  dnorm(0, 0.001)
      sigma0_v[k] ~ dunif(0,1.5)
      tau0_v[k] <- 1/(sigma0_v[k]*sigma0_v[k])
        
        ## Priors for germination and viability trials
        for(i in 1:n_yearViab){
    
        # germination
        mu_g[k,i] ~ dnorm(mu0_g[k], tau0_g[k])
        sigma_g[k,i] ~ dunif(0,1.5)
        tau_g[k,i] <- 1/(sigma_g[k,i]*sigma_g[k,i])
    
        # viability
        mu_v[k,i] ~ dnorm(mu0_v[k], tau0_v[k])
        sigma_v[k,i] ~ dunif(0,1.5)
        tau_v[k,i] <- 1/(sigma_v[k,i]*sigma_v[k,i])
        }
      }
    
    ### Hyperpriors for seed bag burial
    for(k in 1:n_siteBags){
    
      # theta 1
      mu0_1[k] ~  dnorm(0, 0.001)
      sigma0_1[k] ~ dunif(0,1.5)
      tau0_1[k] <- 1/(sigma0_1[k]*sigma0_1[k])
    
      # theta 2
      mu0_2[k] ~  dnorm(0, 0.001)
      sigma0_2[k] ~ dunif(0,1.5)
      tau0_2[k] <- 1/(sigma0_2[k]*sigma0_2[k])
    
      # theta 3
      mu0_3[k] ~  dnorm(0, 0.001)
      sigma0_3[k] ~ dunif(0,1.5)
      tau0_3[k] <- 1/(sigma0_3[k]*sigma0_3[k])
    
        for(i in 1:n_yearBags){
    
        # theta 1
        mu_1[k,i] ~ dnorm(mu0_1[k], tau0_1[k])
        sigma_1[k,i] ~ dunif(0,1.5)
        tau_1[k,i] <- 1/(sigma_1[k,i]*sigma_1[k,i])
    
        # theta 2
        mu_2[k,i] ~ dnorm(mu0_2[k], tau0_2[k])
        sigma_2[k,i] ~ dunif(0,1.5)
        tau_2[k,i] <- 1/(sigma_2[k,i]*sigma_2[k,i])
    
        # theta 3
        mu_3[k,i] ~ dnorm(mu0_3[k], tau0_3[k])
        sigma_3[k,i] ~ dunif(0,1.5)
        tau_3[k,i] <- 1/(sigma_3[k,i]*sigma_3[k,i])
        }
      }
    
    ### Likelihood
    
    ## viability trials
    for(i in 1:n_bag){
    
      # alpha
      alpha_g[i] ~ dnorm(mu_g[siteViab[i],yearViab[i]],tau_g[siteViab[i],yearViab[i]])
      alpha_v[i] ~ dnorm(mu_v[siteViab[i],yearViab[i]],tau_v[siteViab[i],yearViab[i]])
    
      # logit 
      logit(theta_g[i]) <- alpha_g[i]
      logit(theta_v[i]) <- alpha_v[i]
    
      germCount[i] ~ dbinom( theta_g[i] , germStart[i] )
      viabStain[i] ~ dbinom( theta_v[i] , viabStart[i] )
    
    }
    
    ## seed burial experiment
    for(i in 1:n){
    
      # alpha
      alpha_1[i] ~ dnorm(mu_1[siteBags[i],yearBags[i]],tau_1[siteBags[i],yearBags[i]])
      alpha_2[i] ~ dnorm(mu_2[siteBags[i],yearBags[i]],tau_2[siteBags[i],yearBags[i]])
      alpha_3[i] ~ dnorm(mu_3[siteBags[i],yearBags[i]],tau_3[siteBags[i],yearBags[i]])
    
      # logit 
      logit(theta_1[i]) <- alpha_1[i]
      logit(theta_2[i]) <- alpha_2[i]
      logit(theta_3[i]) <- alpha_3[i]
    
      # use this as the number of trials for theta_3
      intactJan[i] = totalJan[i]-seedlingJan[i]
    
      totalJan[i] ~ dbinom(theta_1[i], seedStart[i])
      seedlingJan[i] ~ dbinom(theta_2[i], totalJan[i])
      intactOct[i] ~ dbinom(theta_3[i], intactJan[i])
    
    }
    
    ### Derived quantities
    for(i in 1:n_siteBags){
    
      p.i_g[i] ~ dnorm(mu0_g[i],tau0_g[i])
      logit(p_g[i]) <- p.i_g[i]
    
      p.i_v[i] ~ dnorm(mu0_v[i],tau0_v[i])
      logit(p_v[i]) <- p.i_v[i]
    
      # Probability of viability
      nu_1[i] = p_g[i] + p_v[i]*(1-p_g[i])
    
      p.i_1[i] ~ dnorm(mu0_1[i],tau0_1[i])
      logit(p_1[i]) <- p.i_1[i]
    
      p.i_2[i] ~ dnorm(mu0_2[i],tau0_2[i])
      logit(p_2[i]) <- p.i_2[i]
    
      p.i_3[i] ~ dnorm(mu0_3[i],tau0_3[i])
      logit(p_3[i]) <- p.i_3[i]
    
      # Probability of winter survival
      s1[i] = p_1[i]*(p_2[i] + (1-p_2[i])*(nu_1[i])^(1/3))
    
      # Probability of germination
      g1[i] = p_2[i]/(1-(1-(nu_1[i]^(1/3)))*(1-p_2[i]))
    
      # Probability of spring and summer survival
      s2[i] = p_3[i]*(nu_1[i]^(2/3))
    }
    
    }

### Seedling survival to fruiting

Use \(Y_4\) to estimate the probability of seedling survival to
fruiting.

Explore results using different ways of dealing with undercounting.

This is how I calculated the probability of survival to fruiting in fall
2019. I think this data is amenable to the same approach as the one
above. Rather than using hyperpriors at the site level, I will only get
per-year estimates for all datasets.

    model { 
    # hyperpriors
    for(i in 1:nsites){
    alphaS[i] ~ dnorm(0, .001)
    }
    
    for(i in 1:nyears){
    betaS[i] ~ dnorm(0, .001)
    }
    
    for(j in 1:nsites){
    for(k in 1:nyears){
    gammaS[j,k] ~ dnorm(0, 0.001)
    }
    }
    
    
    # likelihoods
    for (i in 1:N){ 
    phi[i] = ilogit(alphaS[siteSigma[i]] + betaS[yearSigma[i]] + gammaS[siteSigma[i],yearSigma[i]])
    y[i] ~ dbinom(phi[i], n[i])
    
    # simulated data for posterior predictive checks
    y.sim[i] ~ dbinom(phi[i], n[i]) 
    } 
    
    }

### Fruits per plant

Use \(Y_5\) and \(Y_7\) to estimate the number of fruits per plant.

Explore how to combine data from permanent plots vs. all plots at the
site. Perhaps for the fruits per plant estimate use data from all plots
at the site. Perhaps use a mixture model?

Explore how to combine the total fruit equivalents with
undamaged/damaged fruit counts.

    model { 
    # hyperpriors
    for(i in 1:nsites){
    alphaF[i] ~ dnorm(0, .001)
    }
    
    for(i in 1:nyears){
    betaF[i] ~ dnorm(0, .001)
    }
    
    for(j in 1:nsites){
    for(k in 1:nyears){
    gammaF[j,k] ~ dnorm(0, 0.001)
    rF[j,k] ~ dgamma(.001,.001)
    }
    }
    
    # likelihoods
    for (i in 1:N1){
    lambdaF[i] = exp(alphaF[siteFec[i]] + betaF[yearFec[i]] + gammaF[siteFec[i],yearFec[i]])
    y1[i] ~ dnegbin(rF[siteFec[i],yearFec[i]]/(rF[siteFec[i],yearFec[i]]+lambdaF[i]),rF[siteFec[i],yearFec[i]])
    
    # simulated data for posterior predictive checks
    y1.sim[i] ~ dnegbin(rF[siteFec[i],yearFec[i]]/(rF[siteFec[i],yearFec[i]]+lambdaF[i]),rF[siteFec[i],yearFec[i]])
    }
    
    }

### Seeds per fruit

Use \(Y_9\) to estimate seeds per fruit.

    model { 
    # hyperpriors
    for(i in 1:nsites){
    alphaP[i] ~ dnorm(0, .001)
    }
    
    for(i in 1:nyears){
    betaP[i] ~ dnorm(0, .001)
    }
    
    for(j in 1:nsites){
    for(k in 1:nyears){
    gammaP[j,k] ~ dnorm(0, 0.001)
    rP[j,k] ~ dgamma(.001,.001)
    }
    }
    
    # likelihoods
    for (i in 1:N2){
    log(lambdaP[i]) = alphaP[sitePhi[i]] + betaP[yearPhi[i]] + gammaP[sitePhi[i],yearPhi[i]]
    y2[i] ~ dnegbin(rP[sitePhi[i],yearPhi[i]]/(rP[sitePhi[i],yearPhi[i]]+lambdaP[i]),rP[sitePhi[i],yearPhi[i]])
    
    # simulated data for posterior predictive checks
    y2.sim[i] ~ dnegbin(rP[sitePhi[i],yearPhi[i]]/(rP[sitePhi[i],yearPhi[i]]+lambdaP[i]),rP[sitePhi[i],yearPhi[i]])
    }
    
    }

### Transition to seed bank

Combine estimates from belowground vital rates, seedling counts, fruits
per plot, fruits per plant (from plots), and seeds per fruit to estimate
this.

Explore how to do this with/without a seed bank.