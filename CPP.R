#########################################################################
# Calibrated power prior (CPP) 
# 
# CL
#########################################################################


remove(list = ls())
setwd("C:/Users/Caroline/Documents/MEDITWIN")

library(rjags)
library(coda)

N <- 20 #number in population

Y0 <- 2
n0 <- N

Y <- 11
n <- N


distance <- function(Y0, Y, a = 0, b = 1){
  l0 <- length(Y0)
  l <- length(Y)
  s <- max(l0, l)**0.25*as.numeric(ks.test(Y0, Y)[1])
  d <- 1/(1 + exp(a + b*log(s)))
  return(d)
}

alpha <- distance(Y0, Y)

# Data for JAGS
data_jags <- list(
  Y = Y,
  n = n,
  Y0 = Y0,
  n0 = n0,
  alpha = alpha,
  zeros = 0
)

# JAGS model
model_string <- "
model {
  # Prior
  theta ~ dbeta(1, 1)       # uniform prior on p

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
                        variable.names = c("theta"),
                        n.iter = 5000)

summary(samples)
# Extract samples as a matrix
samples_mat <- as.matrix(samples)







