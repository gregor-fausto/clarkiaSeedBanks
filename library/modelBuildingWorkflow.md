
# Model development

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

# Packages

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

# Table of datasets

The table below summarizes the data generated in this study (also
available [as a pdf at this
link](https://github.com/gregor-fausto/clarkiaSeedBanks/blob/master/products/tables/data-summary.pdf)).

|                                        |                  |               |           |
| :------------------------------------- | :--------------: | :-----------: | :-------: |
| **Seed vital rates**                   |        —         |       —       |     —     |
| Seed survival and germination          | Seed bag burial  | Y<sub>1</sub> | 2006-2009 |
| Seed viability                         | Viability trials | Y<sub>2</sub> | 2006-2009 |
| Seed survival and germination          |    Seed pots     | Y<sub>3</sub> | 2013-2019 |
| **Seedling survival** -                |        –         |       —       |     —     |
| Seedling survival to fruiting          |  Field surveys   | Y<sub>4</sub> | 2006-2019 |
| **Fruits per plant**                   |        —         |       —       |     —     |
| Total fruit equivalents per plant      |  Field surveys   | Y<sub>5</sub> | 2006-2012 |
| Undamaged and damaged fruits per plant |  Field surveys   | Y<sub>6</sub> | 2013-2019 |
| Total fruit equivalents per plant      |   Extra plots    | Y<sub>7</sub> | 2006-2012 |
| Undamaged and damaged fruits per plant |   Extra plots    | Y<sub>8</sub> | 2013-2019 |
| **Seeds per fruit**                    |        —         |       —       |     —     |
| Seeds per undamaged fruit              |    Lab counts    | Y<sub>9</sub> | 2006-2019 |
| Seeds per damaged fruit                |    Lab counts    | Y<sub>0</sub> | 2013-2019 |

Links to all the datasets follow:

  - [Seed bag
    burial](https://github.com/gregor-fausto/clarkiaSeedBanks/blob/master/library/dataProcessingWorkflow.md#seed-bag-data)
  - [Viability
    trials](https://github.com/gregor-fausto/clarkiaSeedBanks/blob/master/library/dataProcessingWorkflow.md#viability-trial-data)
  - [Seed
    pots](https://github.com/gregor-fausto/clarkiaSeedBanks/blob/master/library/dataProcessingWorkflow.md#seed-pot-data)
    - needs to be updated
  - [Field surveys for seedling survival to
    fruiting](https://github.com/gregor-fausto/clarkiaSeedBanks/blob/master/library/dataProcessingWorkflow.md#seedlings-and-fruiting-plant-data)
  - [Field surveys for total fruit equivalents per plant, and
    undamaged/damaged fruits per
    plant](https://github.com/gregor-fausto/clarkiaSeedBanks/blob/master/library/dataProcessingWorkflow.md#fruits-per-plant-data-extra-plots)
  - [Lab counts for seeds per fruit from undamaged and damaged
    fruits](https://github.com/gregor-fausto/clarkiaSeedBanks/blob/master/library/dataProcessingWorkflow.md#seeds-per-fruit-data)

# Analysis workflow

  - Load data
  - Rename variables for consistency
  - Filter dataset (e.g. in seedling survival model) until I develop
    measurement model
  - Create data object (list) for analysis in JAGS using for
    e.g. tidybayes \*\* this has its advantages but need to make sure
    that variable designations are consistent
  - Initialize chains \*\* currently this is done by hand but I think I
    need to improve on it \*\* Rather than writing
    list(list(…),list(…),list(…)) each with a set of parameters, I
    want to write functions to generate the appropriate initial
    conditions and then name the inits after they are in lists using a
    function (or something similar)
  - Set conditions for JAGS
  - Set parameters to sample
  - Initialize chains
  - Update chains
  - Sample chains
  - Save samples, original data, and data object for JAGS

Even just writing out these steps I think there is more organization to
be done. For example, I think it’s possible to separate the data
construction steps from the model fits. So, for example, I could use
this file to process the data for analysis in JAGS using a consistent
set of principles. I would then save the objects as data lists to a
folder, which would then be loaded in the scripts with the model.

It would also be great to write R files that are read into here rather
than copying over the JAGS code, which is liable to change.

# Belowground vital rates: age 1

## Data organization

### Seed bag data

Begin by importing the seed bag data (Y<sub>1</sub>).

``` r
seedBagsData <- readRDS("~/Dropbox/dataLibrary/postProcessingData/seedBagsData.rds")
```

The next step is to clean and organize the datasets. Begin by working
with the seed bag data.

  - I exclude rows with missing data in January. Bags that were not
    recovered in January but were recovered in October; we have no way
    of getting an estimate of how many seeds were present halfway
    through the trial.
  - I exclude the row for \``MC III 5 13 1` which has 91 intact seeds
    and 17 seedlings (this should not be possible if 100 seeds started
    the trials).
  - I exclude rows with missing data in October.
  - I exclude rows where the number of intact seeds in January is less
    than the number of intact seeds in October.

The filtering steps in the following code chunk remove 64 rows from the
dataset.

``` r
datasetSize=nrow(seedBagsData)

seedBagsData<-subset(seedBagsData,!is.na(seedBagsData$totalJan))
datasetSize <- cbind(datasetSize,nrow(seedBagsData))

seedBagsData <- seedBagsData %>%
  dplyr::filter(totalJan<=100)
datasetSize <- cbind(datasetSize,nrow(seedBagsData))

seedBagsData<-subset(seedBagsData,!is.na(seedBagsData$intactOct))
datasetSize <- cbind(datasetSize,nrow(seedBagsData))

seedBagsData<-subset(seedBagsData,!(seedBagsData$intactJan<seedBagsData$intactOct))
datasetSize <- cbind(datasetSize,nrow(seedBagsData))
```

I then create a variable for the number of seeds that start the trial.

``` r
# Create a variable for the number of seeds at the start of the trial
seedBagsData$seedStart<-as.double(100)
```

### Viability trial data

Import the viability trial data (Y<sub>2</sub>)

``` r
viabilityRawData <- readRDS("~/Dropbox/dataLibrary/postProcessingData/viabilityRawData.rds")
```

Remove summary variables from the dataset.

``` r
viabilityRawData <- viabilityRawData %>% 
  dplyr::select(-c(germPerc,germNot,viabPerc,viabPerc2,condTest))
```

Create the variable `bag` ()

``` r
viabilityRawData$bag <-as.integer(as.numeric(viabilityRawData$bagNo))
```

One of the rows in the viability dataset is coded differently: the
variable `viabStart` is `NA` instead of `0` as in the rest of the
dataset. Correct this here by recoding that value.

``` r
viabilityRawData[is.na(viabilityRawData$viabStart),]$viabStart = 0
```

I also remove data with any NAs in the viability stain column.

Check the data. Are there any rows where there are more germinants and
seeds that start the viability trials, than seeds that start the
germination trials? The answer is yes; these rows need to be checked in
the data binders.

``` r
viabilityRawData %>% dplyr::filter(germStart - germCount - viabStart<0) 
```

    ## # A tibble: 23 x 11
    ##    site  bagNo round   age block germStart germCount viabStart viabStain
    ##    <chr> <dbl> <dbl> <dbl> <dbl>     <dbl>     <dbl>     <dbl>     <dbl>
    ##  1 CF        8     1     1     5        15         7         9         8
    ##  2 CF        7     1     3    11         3         2         3         0
    ##  3 DLW      29     1     1     5        15         6        10         9
    ##  4 DLW      29     1     1    10        14         5        10         4
    ##  5 FR        1     1     1     9        15        10         6         4
    ##  6 FR       27     1     1    10        15         8         8         1
    ##  7 FR        2     1     2     7        15         7        10         9
    ##  8 GCN      11     1     2     5        15        11         9         8
    ##  9 GCN      47     2     2     9        15         6        10         6
    ## 10 KYE       3     1     1    13        15        10         6         5
    ## # … with 13 more rows, and 2 more variables: viabTotal <dbl>, bag <int>

I exclude the rows with this issue from the analysis.

``` r
viabilityRawData<-viabilityRawData %>% 
  dplyr::filter(germStart - germCount - viabStart >= 0)
```

## Indexing for JAGS

First, I filter the dataset to include only 1-year old seeds

``` r
filterData<-function(x, ageNumber) {
  x %>%
    dplyr::filter(age==ageNumber)
}

seedBagsData1<-filterData(seedBagsData,ageNumber=1)
viabilityRawData1<-filterData(viabilityRawData,ageNumber=1)

seedBagsData2<-filterData(seedBagsData,ageNumber=2)
viabilityRawData2<-filterData(viabilityRawData,ageNumber=2)
```

I assign variables that combine the site, bag number, experimental
round, and age. This creates a unique ID for each bag, for each dataset.
That unique identifier links the bags across the seed bag and viability
trial dataset.

``` r
seedBagsData1<-seedBagsData1 %>%
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="-", remove=FALSE) %>%
  tidyr::unite(col='siteBag', c(site,bagNo), sep="-", remove=FALSE) %>%
  dplyr::mutate(siteBag = as.factor(siteBag)) 

viabilityRawData1<-viabilityRawData1 %>%
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="-", remove=FALSE) %>%
  tidyr::unite(col='siteBag', c(site,bagNo), sep="-", remove=FALSE) %>%
  dplyr::mutate(siteBag = as.factor(siteBag)) 
```

I create a reference table that joins the seed bag and viability trial
data via a unique identifier.

``` r
referenceTable<-data.frame(id=union(seedBagsData1$id, viabilityRawData1$id)) %>%
  dplyr::mutate(idNo = 1:length(id)) 
```

I then append the identifiers to each dataset.

``` r
seedBagsData1<-seedBagsData1 %>%
  dplyr::left_join(referenceTable,by="id")

viabilityRawData1<-viabilityRawData1 %>%
  dplyr::left_join(referenceTable,by="id")
```

## Data for JAGS

Variable names need to be consistent in both datasets. I make the
following changes:

  - Set the `year` variable (`yearStart` in seed bags, `round` in
    viability trials) to be a factor
  - Select variables used in the model fitting
  - Rename variables to uniquely refer to each dataset (e.g. `site`
    renamed as `siteBags` in seed bag dataset)

I also summarize the viability trial dataset. Specifically, I sum each
count (e.g. `germStart`) across all replicates from a bag.

``` r
seedBagsData1 = seedBagsData1 %>%
  dplyr::mutate(year = as.factor(yearStart),
                age = as.factor(age)) %>%
  dplyr::select(site,year,age,totalJan,seedStart,seedlingJan,intactOct) %>%
  dplyr::rename(siteBags = site,
                yearBags = year,
                ageBags = age)

viabilityRawData1 = viabilityRawData1 %>%
  dplyr::mutate(year = as.factor(round),
                age = as.factor(age)) %>%
  dplyr::select(site, year, age, germStart, germCount, viabStart, viabStain, idNo) %>%
  dplyr::rename(siteViab = site,
                yearViab = year,
                ageViab = age,
                bag = idNo) %>%
 # dplyr::mutate(bag = as.factor(bag)) %>%
  dplyr::mutate(germStart = ifelse(is.na(germCount), NA, germStart),
                viabStart = ifelse(is.na(viabStain), NA, viabStart)) %>%
  dplyr::group_by(siteViab,yearViab, ageViab, bag) %>%
  # sum observations in each bag; this is ignoring some variation
  dplyr::summarise(germStart = sum(germStart),
                   germCount = sum(germCount),
                   viabStart = sum(viabStart),
                   viabStain = sum(viabStain))
```

    ## `summarise()` regrouping output by 'siteViab', 'yearViab', 'ageViab' (override with `.groups` argument)

Repeat for age 2 data:

``` r
seedBagsData2<-seedBagsData2 %>%
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="-", remove=FALSE) %>%
  tidyr::unite(col='siteBag', c(site,bagNo), sep="-", remove=FALSE) %>%
  dplyr::mutate(siteBag = as.factor(siteBag)) 

viabilityRawData2<-viabilityRawData2 %>%
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="-", remove=FALSE) %>%
  tidyr::unite(col='siteBag', c(site,bagNo), sep="-", remove=FALSE) %>%
  dplyr::mutate(siteBag = as.factor(siteBag)) 

referenceTable<-data.frame(id=union(seedBagsData2$id, viabilityRawData2$id)) %>%
  dplyr::mutate(idNo = 1:length(id)) 

seedBagsData2<-seedBagsData2 %>%
  dplyr::left_join(referenceTable,by="id")

viabilityRawData2<-viabilityRawData2 %>%
  dplyr::left_join(referenceTable,by="id")

seedBagsData2 = seedBagsData2 %>%
  dplyr::mutate(year = as.factor(yearStart),
                age = as.factor(age)) %>%
  dplyr::select(site,year,age,totalJan,seedStart,seedlingJan,intactOct) %>%
  dplyr::rename(siteBags = site,
                yearBags = year,
                ageBags = age)

viabilityRawData2 = viabilityRawData2 %>%
  dplyr::mutate(year = as.factor(round),
                age = as.factor(age)) %>%
  dplyr::select(site, year, age, germStart, germCount, viabStart, viabStain, idNo) %>%
  dplyr::rename(siteViab = site,
                yearViab = year,
                ageViab = age,
                bag = idNo) %>%
  #dplyr::mutate(bag = as.factor(bag)) %>%
  dplyr::group_by(siteViab,yearViab, ageViab, bag) %>%
    dplyr::mutate(germStart = ifelse(is.na(germCount), NA, germStart),
                viabStart = ifelse(is.na(viabStain), NA, viabStart)) %>%
  # sum observations in each bag; this is ignoring some variation
  dplyr::summarise(germStart = sum(germStart),
                   germCount = sum(germCount),
                   viabStart = sum(viabStart),
                   viabStain = sum(viabStain))
```

    ## `summarise()` regrouping output by 'siteViab', 'yearViab', 'ageViab' (override with `.groups` argument)

``` r
names(seedBagsData2)=paste(names(seedBagsData2),"2",sep="")
names(viabilityRawData2)=paste(names(viabilityRawData2),"2",sep="")
```

## Pass to JAGS

I use `tiybayes::compose_data` to create a data list that I will pass to
JAGS. I reset the number of samples `n` to match the dimension of the
seed bag dataset.

Try to only use complete rows.

``` r
viabilityRawData1<-viabilityRawData1[complete.cases(viabilityRawData1),] %>% dplyr::mutate(bag = as.factor(bag))
viabilityRawData2<-viabilityRawData2[complete.cases(viabilityRawData2),] %>% dplyr::mutate(bag2 = as.factor(bag2))

data <- tidybayes::compose_data(seedBagsData1,viabilityRawData1,
                                seedBagsData2,viabilityRawData2)

data$n = dim(seedBagsData1)[1]
data$n2 = dim(seedBagsData2)[1]
```

## Save for JAGS

``` r
fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/workflow/data/")

saveRDS(data,file=paste0(fileDirectory,"belowgroundDataAgeOneTwo.rds"))

# fileDirectory<- c("/Users/Gregor/Dropbox/dataLibrary/workflow/tidyData/")
# 
# saveRDS(seedBagsData,file=paste0(fileDirectory,"seedBagExperiment.rds"))
# saveRDS(viabilityRawData,file=paste0(fileDirectory,"viabilityExperiment.rds"))
```

Use Y<sub>1</sub> and Y<sub>2</sub> to estimate belowground vital rates.
Refer to the [appendix on using conditional
probability](https://github.com/gregor-fausto/clarkiaSeedBanks/blob/master/products/manuscript/appendix-x-conditional-probability.pdf).

These data are processed and analyzed in
`clarkiaSeedBanks/modelBuild/seedBurial/fullModelScripts.R`.

The JAGS script for this data is in
`clarkiaSeedBanks/modelBuild/jagsScriptsSeedBags/hierarchicalLogitCentered.R`.

Basic process. Use the data from seed bag burial experiment, germination
and viability trials to calculate the parameter theta. Each experiment
is a trial with a binomial likelihood, priors and hyperpriors that are
centered and logit transformed. I then take these probabilities use them
to derive the logit-transformed probabilities, which I then transform
into the probabilities of viability, survival, and germination.

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

I use Y<sub>4</sub> to estimate the probability of seedling survival to
fruiting.

These data are processed and analyzed in
`clarkiaSeedBanks/modelBuild/seedlingSurvival/fullModelScripts.R`.

The JAGS script for this data is in
`clarkiaSeedBanks/modelBuild/jagsScriptsSeedlingSurvival/hierarchicalLogitCentered.R`.

Explore results using different ways of dealing with undercounting.

This is how I calculated the probability of survival to fruiting in fall
2019. I think this data is amenable to the same approach as the one
above. Rather than using hyperpriors at the site level, I will only get
per-year estimates for all datasets.

From `seedBurial/fullModelScripts`:

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

Use Y<sub>5</sub> and Y<sub>7</sub> to estimate the number of fruits per
plant.

These data are processed and analyzed in
`clarkiaSeedBanks/modelBuild/fecundity/fitness_model.R`.

The JAGS script for this data is in
`clarkiaSeedBanks/modelBuild/jagsScriptsFecundity/fecJags.R`.

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

Use Y<sub>9</sub> to estimate seeds per fruit.

These data are processed and analyzed in
`clarkiaSeedBanks/modelBuild/fecundity/seedsPerFruitModelScripts.R`.

The JAGS script for this data is in
`clarkiaSeedBanks/modelBuild/jagsScriptsFecundity/seedsAllJags.R`.

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
