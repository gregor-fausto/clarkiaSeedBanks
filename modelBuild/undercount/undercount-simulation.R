
setwd("~/Dropbox/projects/clarkiaScripts/data/cleanData")
load(file="fecDF.RData")

### UNDERCOUNT SEEDLINGS MODEL

set.seed(10)

# use ALL observations if modeling seedling counts as a latent state
dat <- fecDF %>% 
  dplyr::filter(year == 2007) %>%
  dplyr::select(c(year,site,transect,position,noSeedlings,noFruitingPlants)) %>%
  unique

## SIMULATION

mu = mean(dat$noSeedlings)
sigma = sd(dat$noSeedlings)
kappa = mu^2/(sigma^2-mu)
phi = 0.25
p = 0.5
prob_undercount = 0.05
w = rbinom(600, prob = prob_undercount, size = 1)

N_true <- rnbinom(600, size = kappa, mu = mu)
N_obs <- rbinom(length(N_true), size = N_true, prob = (p*w) + (1-w))
y <- rbinom(length(N_true), size = N_true, prob = phi)

dat %>% dplyr::filter(noFruitingPlants>noSeedlings)
data.frame(N_true,N_obs,y) %>% dplyr::filter(N_obs < N_true)

par(mfrow=c(3,1))
hist(dat$noSeedlings,breaks=100)
hist(N_true,breaks=100)
hist(N_obs, breaks=100)

par(mfrow=c(1,1))
plot(density(dat$noSeedlings),ylim=c(0,.4))
lines(density(N_true),col=2)

discrete.histogram(N_obs, xlab = "y", main = "Simulated data")

data = list(
  N_obs = as.double(N_obs),
  y = as.double(y),
  N = length(y)
)

inits = list(
  list( p1 = runif(1, 0, 1) , p2 = runif(1,0,1), z = rep( 1300, data$N), w = rep( 1, data$N), phi = runif(1,0,1), lambda = rep( 1, data$N), alpha = .1, upsilon = 5),
  list( p1 = runif(1, 0, 1) , p2 = runif(1,0,1), z = rep( 1300, data$N), w = rep( 1, data$N), phi = runif(1,0,1), lambda = rep( 5, data$N), alpha = 2, upsilon = .5),
  list( p1 = runif(1, 0, 1) , p2 = runif(1,0,1), z = rep( 1300, data$N), w = rep( 1, data$N), phi = runif(1,0,1), lambda = rep( 3, data$N), alpha = 10, upsilon = 2) 
)

setwd("~/Dropbox/projects/clarkiaScripts/bayesian-attempts")
sink("simulationUnderJAGS.R")
cat("
    model {
    
    # hyperpriors
    alpha ~ dgamma(0.001, 0.001)
    upsilon ~ dgamma(0.001, 0.001)

    # priors
    p1 ~ dunif (0,1)
    p2 ~ dunif (0,1)
    phi ~ dbeta(1,1)
    for(i in 1:N){ lambda[i] ~ dgamma(alpha, upsilon) } 

    
    # likelihood
    for (i in 1:N){
      w[i] ~ dbern(p1)
      z[i] ~ dpois(lambda[i])
      N_obs[i] ~ dbinom(p2 * w[i] + (1 - w[i]), z[i])
      y[i] ~ dbinom(phi, z[i])
      #y.sim[i] ~ dbinom(p * w[i] + (1 - w[i]), z[i])
    } #end of i
    
   # mean.y = mean(y)
   # mean.y.sim = mean(y.sim)
   # sd.y = sd(y)
   # sd.y.sim = sd(y.sim)
   # p.mean = step(mean.y.sim-mean.y)
   # p.sd = step(sd.y.sim-sd.y)
    
    } #end of model
    ", fill = TRUE)
sink()


# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 1000
n.update = 10000
n.iter = 100000

# Call to JAGS

# tuning (n.adapt)
jm.latent = jags.model("simulationUnderJAGS.R", data = data, inits = inits,
                       n.chains = length(inits), n.adapt = n.adapt)

# burn-in (n.update)
update(jm.latent, n.iter = n.update)

# chain (n.iter)
zc.latent = coda.samples(jm.latent, variable.names = c("p1", "p2", "z"), n.iter = n.iter)

# summary table
MCMCsummary(zc.latent, n.eff = TRUE, params = c("p1", "p2", "z"))

# trace plot
MCMCtrace(zc.latent, n.eff = TRUE, params = c("pi", "alpha", "upsilon"), pdf = FALSE)

