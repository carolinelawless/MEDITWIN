remove(list = ls())

library(rjags)
library(coda)

# Historical data
Y0 <- 12
n0 <- 20

# Current data
Y <- 10
n <- 20

# Power prior weight
a0 <- 0.5  # downweight historical data

# Data for JAGS
data_jags <- list(
  Y = Y,
  n = n,
  Y0 = Y0,
  n0 = n0,
  a0 = a0,
  zeros = 0 # needed for zeros trick
)

# JAGS model
model_string <- "
model {
  # Prior
  theta ~ dbeta(1, 1) # uniform prior on p

  # Likelihood for current data
  Y ~ dbin(theta, n)

  # Zeros trick for historical likelihood contribution
  zeros ~ dpois(phi)
  phi <- - a0 * (Y0 * log(theta) + (n0 - Y0) * log(1 - theta))
}
"

# Run model
model <- jags.model(textConnection(model_string), data = data_jags, n.chains = 3, n.adapt = 1000)
update(model, 1000)
samples <- coda.samples(model, variable.names = c("theta"), n.iter = 5000)

summary(samples)
