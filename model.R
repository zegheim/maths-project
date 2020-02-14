###########################
# CONFIGURATION VARIABLES #
###########################

fname <- "~/Documents/diss/data/df.csv"
seed <- 17071996L
working.dir <- "~/Documents/diss/"

### DO NOT EDIT BELOW THIS LINE ###

###########
# IMPORTS #
###########

library(boot)
library(MASS)
library(Matrix)
library(optimx)
library(pscl)
library(tidyverse)

####################
# HELPER FUNCTIONS #
####################

getLaplMtrx <- function(df, res) {
  rows <- nrow(df)
  G <- Matrix(0, nrow = rows, ncol = rows, sparse = TRUE)
  for (i in seq(rows)) {
    for (j in seq(rows)) {
      if ((i != j) && (abs(df$x[i] - df$x[j]) <= res) && (abs(df$y[i] - df$y[j]) <= res)) {
        G[i, i] <- G[i, i] + 1
        G[i, j] <- -1
      }
    }
  }
  return(G)
}

getPrecMtrx <- function(theta, G) {
  #' Compute precision matrix Q defined as Q = \psi(\kappa^2 + G), where G is the Laplacian Matrix of our observed data.
  #' 
  #' @param theta tuple of (\psi, \kappa), where \psi and \kappa are real-valued.
  #' @param G Laplacian matrix of the observed data.
  #' 
  return(theta[0] * (theta[1]**2 + G))
}

gradLogLik <- function(df, par) {
  Y <- Matrix(df$count)
  Z <- cbind(
    rep(1, length(Y)), 
    as.matrix(subset(df, select = -count))
  )
  
  B_ZERO <- par[1:ncol(Z)]
  B_PLUS <- par[(ncol(Z)+1):(ncol(Z)*2)]
  X <- par[(2*ncol(Z)+1):length(par)]
  
  ZB_ZERO <- Z %*% B_ZERO
  ZB_PLUS <- Z %*% B_PLUS
  
  grad <- Matrix(0, nrow = length(par))
  
  # grads for B_ZERO
  for (i in seq(ncol(Z))) {
    grad[i] <- sum(
      ifelse(Y != 0, exp(-ZB_ZERO) / (1 + exp(-ZB_ZERO)) * Z[, i], 0) -
        ifelse(Y == 0, exp(ZB_ZERO) / (1 + exp(ZB_ZERO)) * Z[, i], 0)
    )
  }
  
  # grads for B_PLUS
  for (i in seq(ncol(Z))) {
    grad[i+ncol(Z)] <- sum(
      ifelse(Y != 0, Y * Z[, i], 0) -
        ifelse(Y != 0, exp(ZB_PLUS) * Z[, i], 0) -
        ifelse(Y != 0, exp(ZB_PLUS + X - exp(ZB_PLUS + X)) / (1 - exp(-exp(ZB_PLUS + X))) * Z[, i], 0)
    )
  }
  
  # grads for X
  for (i in seq(nrow(Z))) {
    grad[i + ncol(Z) * 2] <- sum(
      Y - 
        ifelse(Y != 0, exp(ZB_PLUS + X - exp(ZB_PLUS + X)) / (1 - exp(-exp(ZB_PLUS + X))), 0)
    )
  }
  
  return(grad)
}

gradLogPrior <- function(par, mu, Q) {
  #' Assumes parameters are follows MVNorm(mu, Q^{-1})
  #' 
  #' @param par vector containing \beta_0's, \beta_+'s, and X_i's
  #' @param mu mean vector
  #' @param Q precision matrix / inverse covariance Sigma^{-1}
  #' 
  parminmu <- Matrix(par - mu, ncol = 1)
  return(-Q %*% parminmu)
}

hessLogLik <- function(df, par) {
  Y <- Matrix(df$count)
  Z <- cbind(
    rep(1, length(Y)), 
    as.matrix(subset(df, select = -count))
  )
  
  B_ZERO <- par[1:ncol(Z)]
  B_PLUS <- par[(ncol(Z)+1):(ncol(Z)*2)]
  X <- par[(2*ncol(Z)+1):length(par)]
  
  ZB_ZERO <- Z %*% B_ZERO
  ZB_PLUS <- Z %*% B_PLUS
  ZBPX <- ZB_PLUS + X
  
  hess <- Matrix(0, nrow = length(par), ncol = length(par))
  
  # hessian for B_ZERO
  for (i in seq(ncol(Z))) {
    hess[i, i] <- sum(
      ifelse(Y != 0, exp(ZB_ZERO) / (1 + exp(ZB_ZERO))**2 * Z[, i]**2, 0) -
        ifelse(Y == 0, exp(-ZB_ZERO) / (1 + exp(-ZB_ZERO))**2 * Z[, i]**2, 0)
    )
  }
  
  # hessian for B_PLUS
  for (i in seq(ncol(Z))) {
    for (j in seq(ncol(Z))) {
      hess[ncol(Z) + i, ncol(Z) + j] <- sum(
        ifelse(Y != 0, exp(ZB_PLUS) * Z[, i]**2, 0) -
          ifelse(Y != 0, (exp(ZBPX - exp(ZBPX)) * (1 - exp(ZBPX)) - exp(-exp(ZBPX)) * Z[, i]*Z[, j]) / (1 - exp(-ZBPX))**2, 0)
      )
    }
  }
  
  # hessian for B_PLUS w/ X
  for (i in seq(ncol(Z))) {
    hess[ncol(Z) + i, (2*ncol(Z)+1):(2*ncol(Z)+nrow(Z))] <- -1 * sum(
      ifelse(
        Y != 0, 
        exp(ZBPX - exp(ZBPX)) / (1 - exp(-ZBPX))**2 * (exp(ZBPX - exp(ZBPX)) * (1 - exp(ZBPX)) - exp(-exp(ZBPX))) * Z[, i], 
        0
      )
    )
  }
  
  return(forceSymmetric(hess))
}

hessLogPrior <- function(Q) {
  return(-Q)
}

logLik <- function(par, Y, Z) {
  ZB_PLUS <- as.matrix(Z) %*% par[1:5]
  ZB_ZERO <- as.matrix(Z) %*% par[6:10]
  sum(
    ifelse(Y == 0, log(1 - inv.logit(ZB_ZERO)), 0) +
    ifelse(Y != 0, log(inv.logit(ZB_ZERO)) - exp(ZB_PLUS) + (Y*ZB_PLUS) - lgamma(Y + 1) - log(1 - exp(-exp(ZB_PLUS))), 0)
  )                                    
}

margPost <- function(df, par, theta, G) {
  Y <- Matrix(df$count)
  Z <- cbind(
    rep(1, length(Y)), 
    as.matrix(subset(df, select = -count))
  )
  
  # initial values
  Q <- getPrecMtrx(theta, G)
  grad.ll <- gradP
  
  
  while(TRUE) {
    
  }
}

########
# MAIN #
########

setwd(working.dir)
set.seed(seed)

df <- read.csv(fname)
# Normalise all covariates to have mean 0 and sd 1
df[-3] <- scale(df[-3])

Y <- df$count
# Add intercept column
Z <- cbind(rep(1, length(Y)), df[-3])
colnames(Z)[1] <- "intercept"

# Initialise parameters using the coefficients of a GLM fit
par_init <- as.numeric(c(
  glm.fit(Z, Y, family = poisson())$coefficients,
  glm.fit(Z, factor(Y > 0), family = binomial())$coefficients
))

opt <- optimx(par_init, logLik, Y = Y, Z = Z, method = "BFGS", itnmax = 10000, control=list(maximize = TRUE))

# Check with using library
fire.hurdle <- hurdle(count ~ x + y + elevation + max.temp, data = df, separate = FALSE)

# Compare final value
opt$value; fire.hurdle$optim$value

# Compare coefficients
B_PLUS_hat <- as.numeric(tail(opt, 1)[1:10])[1:5]
B_ZERO_hat <- as.numeric(tail(opt, 1)[1:10])[6:10]
as.numeric(fire.hurdle$coefficients$count); B_PLUS_hat
as.numeric(fire.hurdle$coefficients$zero); B_ZERO_hat

# Calculate Hessian
inf_matrix <- -1 * optimHess(par = as.numeric(opt[1:10]), fn = poiHuLogLik, Y = Y, Z = Z, control=list(fnscale = -1))
inf_matrix[1:5, 6:10] <- 0
inf_matrix[6:10, 1:5] <- 0
inv.inf_matrix <- solve(inf_matrix)

# generating samples of MLE
mles <- mvrnorm(n = 100, mu = as.numeric(tail(opt, 1)[1:10]), Sigma = inv.inf_matrix)

