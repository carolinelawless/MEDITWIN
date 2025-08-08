#########################################################################
# Bayesian hierarchical model (BHM) 
# 
# CL
#########################################################################


remove(list = ls())
setwd("C:/Users/Caroline/Documents/MEDITWIN")

library(rjags)
library(coda)

N <- 20 #number in population for each group

# Data
Y <- c(1, 10, 20)  # successes
n <- rep(N, length(Y))   # trials
K <- length(Y) # number of arms

data_jags <- list(Y = Y, n = n, K = K)

model_string <- "
model {
  mu ~ dnorm(0, 0.001)
  tau ~ dgamma(0.001, 0.001)
  for (k in 1:K) {
    theta[k] ~ dnorm(mu, tau)
    logit(p[k]) <- theta[k]
    Y[k] ~ dbin(p[k], n[k])
  }
}
"

model <- jags.model(textConnection(model_string), data = data_jags, n.chains = 3, n.adapt = 1000)
update(model, 1000)
samples <- coda.samples(model, variable.names = c("mu", "tau", "theta", "p"), n.iter = 5000)
# summary(samples)

# Extract samples as a matrix
samples_mat <- as.matrix(samples)




posterior_predictive <- function(samples_mat, N, K){
  draw <- sample(1: nrow(samples_mat), 1)
  mu1 <- samples_mat[draw, "mu"]
  tau1 <- samples_mat[draw, "tau"]
  y <- vector(length = K) 
  for(k in 1:K){
    theta1 <- rnorm(1, mu1, tau1)
    p <- plogis(theta1)
    y[k] <- rbinom(1, N, p)
  }
  y <- y[order(y)]
  return(y)
}


posterior_predictive(samples_mat, N, K)

n_post_pred <- 1e2
post_pred_samples <- matrix(nrow = n_post_pred, ncol = K)

for(i in 1:n_post_pred){
  sample <- posterior_predictive(samples_mat, N, K)
  post_pred_samples[i,] <- sample
}

