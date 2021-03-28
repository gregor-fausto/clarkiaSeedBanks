rm(list = setdiff(ls(all = TRUE), c("scriptConvergenceDirectory", "fileDirectory", 
  "outputDirectory")))  # if using in source(script)
options(stringsAsFactors = FALSE)

library(MCMCvis)
library(tidybayes)
library(tidyverse)
library(magrittr)
library(bayesplot)
library(rethinking)

# get fitted model directory
mcmcSampleFiles <- paste0(fileDirectory, list.files(fileDirectory))

# get MCMC samples
mcmcSamples <- readRDS(mcmcSampleFiles[[grep("seedsPerFruitSamples.rds", mcmcSampleFiles)]])

# --- Convergence diagnostics
# ------------------------------------------------------------------- ---

# --- Graphical checks ---- ---

# function to produce trace plots as JPEGs per population
f = function(x = "parm", model = "modelName", jpeg.quality = 75) {
  # get chains for parameter as a list
  
  chains <- MCMCchains(mcmcSamples, params = x, mcmc.list = TRUE)
  chains.list <- lapply(chains, as.matrix)
  
  n = (dim(chains.list[[1]])[2])/20
  counter <- 1
  parm = x
  while (counter <= n) {
    jpeg(filename = paste0(outputDirectory, "convergence-", model, "-", parm, 
      "-", counter, ".jpeg"), quality = 75)
    
    par.names = colnames(chains[[1]])
    par(mfrow = c(4, 5), oma = c(5, 4, 0, 0) + 0.1, mar = c(0, 0, 1, 1) + 0.1)
    for (i in 1:20) {
      plot(as.vector(chains[[1]][, (20 * (counter - 1) + i)]), type = "n", 
        axes = FALSE, frame = FALSE, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
      # three chains
      lines(as.vector(chains[[1]][, (20 * (counter - 1) + i)]), col = rgb(0.9882353, 
        0.5529412, 0.3490196, 0.5))
      lines(as.vector(chains[[2]][, (20 * (counter - 1) + i)]), col = rgb(1, 
        1, 0.75, 0.5))
      lines(as.vector(chains[[3]][, (20 * (counter - 1) + i)]), col = rgb(0.5686275, 
        0.7490196, 0.8588235, 0.5))
      # add parameter name to jpeg
      text(0.5 * length(as.vector(chains[[1]][, (20 * (counter - 1) + i)])), 
        0.9 * max(as.vector(chains[[1]][, (20 * (counter - 1) + i)])), par.names[(20 * 
          (counter - 1) + i)])
    }
    dev.off()
    counter <- counter + 1
  }
  
}

# print trace plots for each parameter
f(x = "nu_seeds", model = "seedsUndamaged")
f(x = "sigma0_seeds", model = "seedsUndamaged")
f(x = "sigma_seeds", model = "seedsUndamaged")

# print trace plots for each parameter
f(x = "nu_dam_seeds", model = "seedsDamaged")
f(x = "sigma0_dam_seeds", model = "seedsDamaged")
f(x = "sigma_dam_seeds", model = "seedsDamaged")

# recover chains for Undamaged seeds
all.chains = MCMCchains(mcmcSamples, params = c("nu_seeds", "sigma0_seeds", "sigma_seeds"), mcmc.list = TRUE)

# R-hat for all parameters
par(mfrow = c(1, 1))
rhat = coda::gelman.diag(all.chains, confidence = 0.95)$psrf[, 1]

# calculate Heidelberg diagnostic...
hd = coda::heidel.diag(all.chains)
p = c()
for (i in 1:3) {
  p[i] = sum(hd[[i]][, 1])/length(hd[[i]][, 1])
}
p = signif(p, 2)
dt = hist(coda::gelman.diag(all.chains, confidence = 0.95)$psrf[, 1], breaks = 25, 
  plot = FALSE)

# print histogram
jpeg(filename = paste0(outputDirectory, "rhat-heidelberg-seeds-undamaged", ".jpeg"), 
  quality = 75)
hist(rhat, col = "black", border = "white", breaks = 25, main = "Distribution of R-hat; undamaged seeds")
plot.hist = hist(rhat, breaks = 25, plot = FALSE)
prob = sum(rhat < 1.05)/length(rhat)
text(max(plot.hist$mids), max(plot.hist$counts) * 0.9, paste0("Percent of R-hat < 1.05: ", 
  signif(prob, 2)), pos = 2)
text(0.999 * max(dt$mids), 0.85 * max(dt$counts), paste0("% passing (chain 1): ", 
  p[1]), pos = 2)
text(0.999 * max(dt$mids), 0.8 * max(dt$counts), paste0("% passing (chain 2): ", 
  p[2]), pos = 2)
text(0.999 * max(dt$mids), 0.75 * max(dt$counts), paste0("% passing (chain 3): ", 
  p[3]), pos = 2)
dev.off()


# recover chains for Damaged seeds
all.chains = MCMCchains(mcmcSamples, params = c("nu_dam_seeds", "sigma0_dam_seeds", "sigma_dam_seeds"), mcmc.list = TRUE)

# R-hat for all parameters
par(mfrow = c(1, 1))
rhat = coda::gelman.diag(all.chains, confidence = 0.95)$psrf[, 1]

# calculate Heidelberg diagnostic...
hd = coda::heidel.diag(all.chains)
p = c()
for (i in 1:3) {
  p[i] = sum(hd[[i]][, 1])/length(hd[[i]][, 1])
}
p = signif(p, 2)
dt = hist(coda::gelman.diag(all.chains, confidence = 0.95)$psrf[, 1], breaks = 25, 
          plot = FALSE)

# print histogram
jpeg(filename = paste0(outputDirectory, "rhat-heidelberg-seeds-damaged", ".jpeg"), 
     quality = 75)
hist(rhat, col = "black", border = "white", breaks = 25, main = "Distribution of R-hat; damaged seeds")
plot.hist = hist(rhat, breaks = 25, plot = FALSE)
prob = sum(rhat < 1.05)/length(rhat)
text(max(plot.hist$mids), max(plot.hist$counts) * 0.9, paste0("Percent of R-hat < 1.05: ", 
                                                              signif(prob, 2)), pos = 2)
text(0.999 * max(dt$mids), 0.85 * max(dt$counts), paste0("% passing (chain 1): ", 
                                                         p[1]), pos = 2)
text(0.999 * max(dt$mids), 0.8 * max(dt$counts), paste0("% passing (chain 2): ", 
                                                        p[2]), pos = 2)
text(0.999 * max(dt$mids), 0.75 * max(dt$counts), paste0("% passing (chain 3): ", 
                                                         p[3]), pos = 2)
dev.off()
