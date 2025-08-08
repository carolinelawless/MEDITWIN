remove(list = ls())
setwd("C:/Users/Caroline/Documents/MEDITWIN")

library(rjags)
library(coda)

N <- 20 #number in population for each group

# Data
Y <- c(1, 10, 20)  # successes
n <- rep(N, length(Y))   # trials
K <- length(Y) # number of arms
phi <- c(0.1, 0.9)              # prespecified low and high response rates

data_jags <- list(Y = Y, n = n, K = K, phi = phi)

# model parameters: gamma1, gamma2, tau1 (within group variance), tau2 (how many in each group)
# latent variables: theta, I, eta, p
# observtion: Y

model_string <- "
model {
  for (j in 1:2) {
    gamma[j] <- logit(phi[j])
  }
  
  tau1 ~ dgamma(0.001, 0.001)
  #tau2 ~ dgamma(0.001, 0.001) #(I think it makes more sense to put a prior on the location rather than the variance here)
  tau2 ~ dnorm(0, 1)
  
  for (k in 1:K) {
    theta[k] ~ dnorm(tau2, 1)
    I[k] <- step(theta[k]) + 1
    eta[k] ~ dnorm(gamma[I[k]], tau1)
    logit(p[k]) <- eta[k]
    Y[k] ~ dbin(p[k], n[k])
  }
}
"
writeLines(model_string, con = "classification_model.txt")


model <- jags.model("classification_model.txt", data = data_jags, n.chains = 3, n.adapt = 1000)
update(model, 1000)  # burn-in

params <- c("gamma", "tau1", "tau2", "theta", "I", "p")
samples <- coda.samples(model, variable.names = params, n.iter = 5000)
# Extract samples as a matrix
samples_mat <- as.matrix(samples)

posterior_predictive <- function(samples_mat, N, K){
  draw <- sample(1: nrow(samples_mat), 1)
  
  # I <- samples_mat[draw , grepl("^I\\[\\d+\\]$", names(samples_mat))]
  gamma <- samples_mat[draw, grepl("^gamma\\[\\d+\\]$", colnames(samples_mat))]

  tau1 <- samples_mat[draw, "tau1"]
  tau2 <- samples_mat[draw, "tau2"]
  
  y <- vector(length = K) 
  for (k in 1:K) {
    #theta <- rnorm(1, 0, tau2)
    theta <- rnorm(1, tau2, 1)
    I <- (theta > 0) + 1
    eta <- rnorm(1, gamma[I], tau1)
    p <- plogis(eta)
    y[k] <- rbinom(1, N, p)
  }
  y <- y[order(y)]
  return(y)
}

