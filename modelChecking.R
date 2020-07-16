## MODEL CHECKING SCRIPTS ##

library(bayesplot)

# directory of posteriors
directory = "/Users/Gregor/Dropbox/dataLibrary/posteriors/"

# read in posterior
belowgroundSamples <- readRDS(paste0(directory, "belowgroundSamplesAllYears.RDS"))

# extract Brooks-Gelman-Rubin statistic (R-hat) 
rhats<-MCMCsummary(belowgroundSamples)$Rhat


# this is a bayesplot option
color_scheme_set("brightblue")

# plot R-hat values
mcmc_rhat_hist(rhats)

# check trace plots? MCMCtrace(belowgroundSamples,params='sigma0_1')

# Heidelberger-Welch Diagnostic

heidel.diag(belowgroundSamples)

