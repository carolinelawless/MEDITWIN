#########################################################################
# Modified power prior (PP) 
# 
# CL 08/2025
#########################################################################


remove(list = ls())
setwd("C:/Users/Caroline/Documents/MEDITWIN")

library(rjags)
library(coda)

N <- 20 #number in population


# Generate example data
set.seed(123)
Y0 <- rbinom(1, N, 0.2)
Y  <- rbinom(1, N, 0.8)


# Data for JAGS
data_jags <- list(
  Y = Y,
  n = N,
  Y0 = Y0,
  n0 = N,
  zeros = 0
)

# JAGS model
model_string <- "
model {
  # Prior
  theta ~ dbeta(1, 1)       # uniform prior on p
  alpha ~ dbeta(1, 1)

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

head(samples_mat)

mean(samples_mat[,"theta"])

mean(samples_mat[,"alpha"])

