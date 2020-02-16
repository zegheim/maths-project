###########################
# CONFIGURATION VARIABLES #
###########################

fname <- "~/Documents/diss/data/df_hires.csv"
res <- 0.1
seed <- 17071996L
working.dir <- "~/Documents/diss/"

### DO NOT EDIT BELOW THIS LINE ###

###########
# IMPORTS #
###########

library(INLA)
library(MASS)
library(Matrix)
library(optimx)
library(pscl)
library(tidyverse)

####################
# HELPER FUNCTIONS #
####################

getLaplMtrx <- function(df, res, verbose=FALSE) {
  rows <- nrow(df)
  
  if (verbose) {
    pb <- txtProgressBar(0, rows, style = 3)
  }
 
  G <- Matrix(0, nrow = rows, ncol = rows, sparse = TRUE)
  for (i in seq(rows)) {
    if (verbose) {
      setTxtProgressBar(pb, i)
    }
    
    nghbrs <- which(abs(df$x - df$x[i]) <= res & abs(df$y - df$y[i]) <= res)
    for (nghbr in nghbrs){
      if (nghbr != i) {
        G[i, nghbr] <- -1
        G[i, i] <- G[i, i] + 1
      }
    }
  }
  
  if (verbose) {
    close(pb)
  }
  
  return(G)
}

getPrecMtrx <- function(theta, G) {
  #' Compute precision matrix Q defined as Q = \psi(\kappa^2 + G), where G is the Laplacian Matrix of our observed data.
  #' 
  #' @param theta tuple of (\log(\psi), \log(\kappa)), where \log(\psi) and \log(\psi) are real-valued.
  #' @param G Laplacian matrix of the observed data.
  #' 
  return(exp(theta[1]) * (exp(theta[2])**2 * Diagonal(nrow(G), 1) + G))
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
  ZBPX <- ZB_PLUS + X
  
  # common intermediate (temporary) term
  tmp <- Y - exp(ZBPX) - exp(ZBPX) / (exp(exp(ZBPX)) - 1) 
  
  grad <- Matrix(0, nrow = length(par))
  for (i in seq(ncol(Z))) {
    # grads for B_ZERO
    grad[i] <- sum(ifelse(Y != 0, 1 / (1 + exp(-ZB_ZERO)), -1 / (1 + exp(ZB_ZERO))) * Z[,i])
    # grads for B_PLUS
    grad[i+ncol(Z)] <- sum(ifelse(Y != 0, tmp, 0) * Z[,i])
  }
  # grads for X
  grad[(2*ncol(Z)+1):length(par)] <- ifelse(Y != 0, tmp, 0)
  
  return(grad)
}

gradLogPrior <- function(par, mu, Q) {
  #' Assumes parameters are follows MVNorm(mu, Q^{-1})
  #' 
  #' @param par vector containing \beta_0's, \beta_+'s, and X_i's
  #' @param mu mean vector
  #' @param Q precision matrix / inverse covariance Sigma^{-1}
  #' 
  parminmu <- Matrix(par - mu)
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
  
  # intermediate term for B_ZERO
  tmp_b0 <- 1 / (exp(ZB_ZERO) + exp(-ZB_ZERO) + 2)
  
  # hessian for B_ZERO
  for (i in seq(ncol(Z))) {
    for (j in seq(ncol(Z))) {
      hess[i, j] <- sum(Z[,i] * Z[,j] * tmp_b0 * ((Y != 0) - (Y == 0)))
    }
  }
  
  # intermediate term for B_PLUS & B_PLUS w/ X
  tmp_bx <- -exp(ZBPX) - (exp(ZBPX) * ((exp(ZBPX) - 1) * exp(exp(ZBPX)) + 1)) / (exp(exp(ZBPX)) - 1)**2

  # hessian for B_PLUS
  for (i in seq(ncol(Z))) {
    for (j in seq(ncol(Z))) {
      hess[ncol(Z) + i, ncol(Z) + j] <- sum(ifelse(Y != 0,  Z[,i] * Z[,j] * tmp_bx, 0))
    }
  }
  
  # hessian for B_PLUS w/ X
  for (i in seq(ncol(Z))) {
    hess[ncol(Z) + i, (2*ncol(Z)+1):(2*ncol(Z)+nrow(Z))] <- ifelse(Y != 0, Z[,i] * tmp_bx, 0)
  }
  
  return(forceSymmetric(hess))
}

hessLogPrior <- function(Q) {
  return(-Q)
}

logLik <- function(df, par) {
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
  
  out <- sum(
    ifelse(Y != 0, 
      Y * ZBPX - log(1 - exp(-exp(ZBPX))) - exp(ZBPX) - log(1 + exp(-ZB_ZERO)), 
      -log(1 + exp(ZB_ZERO))
    )
  )     
  
  return(out)
}

logPrior <- function(par, mu, Q) {
  k <- length(par)
  parminmu <- Matrix(par - mu)
  L <- Matrix::expand(Cholesky(Q, perm = TRUE))$L
  halflogdetQ <- sum(log(diag(L)))
  out <- as.numeric(
    halflogdetQ - (k / 2) * log(2 * pi) - .5 * sum((parminmu) * (Q %*% parminmu))
  )
  return(out)
}

margPost <- function(df, X, par, G, verbose = FALSE, acc = 1e-7, return.X = FALSE) {
  Y <- Matrix(df$count)
  Z <- cbind(
    rep(1, length(Y)), 
    as.matrix(subset(df, select = -count))
  )
  
  mu <- rep(0, length(X))
  Q.X <- getPrecMtrx(par, G)
  Q.beta <- Diagonal(2 * ncol(Z), 1e-6)
  Q <- rbind(
    cbind(Q.beta, Matrix(0, nrow = nrow(Q.beta), ncol = ncol(Q.X))),
    cbind(Matrix(0, nrow = nrow(Q.X), ncol = ncol(Q.beta)), Q.X)
  )
  
  obj.prop <- obj.curr <- -(logLik(df, X) + logPrior(X, mu, Q))
  
  scaled.acc <- max(abs(obj.curr), 1) * acc
  tol.test <- Inf
  
  # Newton-Raphson optimisation to find X that
  # maximises log-likelihood of the full conditional
  # (i.e. find the roots of the full conditional gradient)
  while(TRUE) {
    alpha <- 1
    X.prop <- X
    obj.prop <- obj.curr
    grad.fc <- gradLogLik(df, X) + gradLogPrior(X, mu, Q)
    hess.fc <- hessLogLik(df, X) + hessLogPrior(Q)
    diff <- -Matrix::solve(hess.fc, grad.fc)
    X <- X.prop + alpha * diff
    obj.curr <- -(logLik(df, X) + logPrior(X, mu, Q))
    
    if (verbose) {
      print(paste("obj.curr = ", obj.curr))
      print(paste("obj.curr > obj.prop = ", (obj.curr > obj.prop)))
    }
    
    while (obj.curr > obj.prop) {
      alpha <- alpha / 2
      X <- X.prop + alpha * diff
      obj.curr <- -(logLik(df, X) + logPrior(X, mu, Q))
      
      if (verbose) {
        print(paste("alpha = ", alpha, ", obj.curr - obj.prop = ", round(obj.curr - obj.prop, 2)))
        print(paste("obj.curr = ", obj.curr))
      }
      
      if (obj.curr < obj.prop && verbose) {
        print("REDUCED STEP SIZE REACHED.")
      }
      
      if (alpha <= 2**(-9)) {
        X <- X.prop
        obj.curr <- -(logLik(df, X) + logPrior(X, mu, Q))
        
        if (verbose) {
          print("STOP MOVING.")
        }
        break
      }
    }
    
    scaled.acc <- max(abs(obj.curr), 1) * acc
    tol.test <- sqrt(sum(grad.fc**2))
    
    if (verbose) {
      print(paste("obj.curr = ", obj.curr, ", obj.prop = ", obj.prop))
      print(paste("tol.test = ", round(tol.test, 3), ", scaled.acc = ", round(scaled.acc, 3)))
    }
    
    nmv <- abs(sqrt(sum(X**2)) - sqrt(sum(X.prop**2)))
    
    if (tol.test < scaled.acc || nmv < acc) {
      break
    }
  }
  
  if (verbose) {
    print("EXITED WHILE LOOP.")
  }
  
  if (return.X) {
    return(X)
  } else {
    # logPrior(X, X, Q) is the Laplace approximation of
    # the full conditional evaluated at X_hat(theta). 
    # We also assume that the log prior of theta is 0.
    obj <- as.numeric(-(0 + logPrior(X, mu, Q) + logLik(df, X) - logPrior(X, X, Q)))
    counter <<- counter+1
    print(paste("counter = ", counter))
    return(obj)
  }
}

########
# MAIN #
########

setwd(working.dir)
set.seed(seed)

df <- read.csv(fname)
Y <- df$count
Z <- cbind(rep(1, length(Y)), df[-3])
colnames(Z)[1] <- "intercept"

X <- rep(0, 2 * ncol(df) + nrow(df))
# X[1:(2*ncol(df))] <- as.numeric(c(
#   glm.fit(Z, Y, family = poisson())$coefficients,
#   glm.fit(Z, factor(Y > 0), family = binomial())$coefficients
# ))
theta <- c(0.01, 0.01)
G <- getLaplMtrx(df, res, verbose = TRUE)

counter <- 0
opt <- optimx(theta, margPost, df = df, X = X, G = G, method = "Nelder-Mead", itnmax = 1000, control = list(kkt = FALSE))

X_hat <- margPost(df, X, c(opt$p1, opt$p2), G, return.X = TRUE)

# Comparison with model w/o latent spatial effects
fire.hurdle <- hurdle(count ~ x + y + elevation + avg.temp, data = df, separate = FALSE)
as.numeric(c(fire.hurdle$coefficients$zero, fire.hurdle$coefficients$count)); X_hat[1:10]



##########################

Y <- df$count
# Add intercept column
Z <- cbind(rep(1, length(Y)), df[-3])
colnames(Z)[1] <- "intercept"

# Initialise parameters using the coefficients of a GLM fit
par_init <- as.numeric(c(
  glm.fit(Z, Y, family = poisson())$coefficients,
  glm.fit(Z, factor(Y > 0), family = binomial())$coefficients
))

opt_simple <- optimx(par_init, logLik, Y = Y, Z = Z, method = "BFGS", itnmax = 10000, control=list(maximize = TRUE))

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
inf_matrix <- -1 * optimHess(par = as.numeric(opt[1:10]), fn = logLik, Y = Y, Z = Z, control=list(fnscale = -1))
inf_matrix[1:5, 6:10] <- 0
inf_matrix[6:10, 1:5] <- 0
inv.inf_matrix <- solve(inf_matrix)

# generating samples of MLE
mles <- mvrnorm(n = 100, mu = as.numeric(tail(opt, 1)[1:10]), Sigma = inv.inf_matrix)

