#########################################################################
# Bayesian modified power prior
# CL
#########################################################################

remove(list = ls())

library(rjags)
library(coda)

# Historical data
Y0 <- 12
n0 <- 20

# Current data
Y <- 10
n <- 20

# Data for JAGS
data_jags <- list(
  Y = Y,
  n = n,
  Y0 = Y0,
  n0 = n0,
  zeros = 0 # needed for zeros trick
)

# JAGS model
model_string <- "
model {
  # Prior
  theta ~ dbeta(1, 1)       # uniform prior on p
  alpha ~ dbeta(1, 1)       # power prior weight

  # Likelihood for current data
  Y ~ dbin(theta, n)

  # Zeros trick for historical likelihood contribution
  zeros ~ dpois(phi)
  phi <- - alpha * (Y0 * log(theta + 1.0E-10) + (n0 - Y0) * log(1 - theta + 1.0E-10))
}
"

# Run model
model <- jags.model(textConnection(model_string),
                    data = data_jags,
                    n.chains = 3,
                    n.adapt = 1000)

update(model, 1000)
samples <- coda.samples(model,
                        variable.names = c("theta", "alpha"),
                        n.iter = 5000)

summary(samples)
# Extract samples as a matrix
samples_mat <- as.matrix(samples)
