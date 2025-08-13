#########################################################################
# Calibrated power prior (CPP) 
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
Y0_vec <- rbinom(N, 1, 0.2)
Y_vec  <- rbinom(N, 1, 0.8)


alpha_calibration <- function(Y0, Y, a = 0, b = 1){
  l0 <- length(Y0)
  l <- length(Y)
  s <- max(l0, l)**0.25*as.numeric(ks.test(Y0, Y)[1])
  alpha <- 1/(1 + exp(a + b*log(s)))
  return(alpha)
}


# Calibration
alpha <- alpha_calibration(Y0_vec, Y_vec)



# Data for JAGS
data_jags <- list(
  Y = sum(Y_vec),
  n = length(Y_vec),
  Y0 = sum(Y0_vec),
  n0 = length(Y0_vec),
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

head(samples_mat)

mean(samples_mat[,"theta"])





