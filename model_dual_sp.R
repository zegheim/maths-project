##### CONFIGURATION VARIABLES #####

fname <- "df_ina_lores.csv"
res <- 0.5
working.dir <- "~/Documents/diss"

data.dir <- str_glue("{working.dir}/data")
proj.str <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
seed <- 17071996L

##### DO NOT EDIT BELOW THIS LINE ####

##### IMPORTS #####

library(extraDistr)
library(gridExtra)
library(INLA)
library(latex2exp)
library(MASS)
library(Matrix)
library(optimx)
library(pals)
library(profvis)
library(pscl)
library(raster)
library(rasterVis)
library(RColorBrewer)
library(rgeos)

##### HELPER FUNCTIONS #####

getConfInt <- function(result, alpha = 0.05) {
  Sigma_hat <- Matrix::solve(result$Q_hat, Matrix::Diagonal(nrow(result$Q_hat)))
  conf.int <- data.frame(
    lwr = as.numeric(result$X_hat) + qnorm(alpha/2) * sqrt(diag(Sigma_hat)),
    upr = as.numeric(result$X_hat) + qnorm(1 - alpha/2) * sqrt(diag(Sigma_hat))
  )
  conf.int$significant <- !(conf.int$lwr <= 0 & conf.int$upr >= 0)
  return(conf.int)
}

getLaplMtrx <- function(coords, res, verbose = FALSE) {
  rows <- nrow(coords)
  
  if (verbose) {
    pb <- txtProgressBar(0, rows, style = 3)
  }
 
  G <- Matrix(0, nrow = rows, ncol = rows, sparse = TRUE)
  for (i in seq(rows)) {
    if (verbose) {
      setTxtProgressBar(pb, i)
    }
    
    nghbrs <- which(
      round(abs(coords$x - coords$x[i]) + abs(coords$y - coords$y[i]), digits = 2) == res)
    if (length(nghbrs) > 0) {
      G[i, nghbrs] <- -1
      G[i, i] <- length(nghbrs) 
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

gradLogLik <- function(Y, Z, X) {
  idx_BZ <- 1:ncol(Z)
  idx_BP <- (ncol(Z)+1):(ncol(Z)*2)
  idx_UZ <- (2*ncol(Z)+1):(2*ncol(Z)+nrow(Z))
  idx_UP <- (2*ncol(Z)+nrow(Z)+1):length(X)
  B_ZERO <- X[idx_BZ]
  B_PLUS <- X[idx_BP]
  U_ZERO <- X[idx_UZ]
  U_PLUS <- X[idx_UP]
  ZBU_ZERO <- Z %*% B_ZERO + U_ZERO
  ZBU_PLUS <- Z %*% B_PLUS + U_PLUS
  # common intermediate (temporary) terms
  TMP_ZERO <- (Y != 0) / (1 + exp(ZBU_ZERO)) - (Y == 0) / (1 + exp(-ZBU_ZERO))
  TMP_PLUS <- Matrix(
    (Y != 0) * (Y - exp(ZBU_PLUS) - exp(ZBU_PLUS) / (exp(exp(ZBU_PLUS)) - 1)),
    sparse = FALSE
  )
  
  grad <- Matrix(0, nrow = length(X))
  grad[idx_BZ] <- t(Z) %*% TMP_ZERO
  grad[idx_BP] <- t(Z) %*% TMP_PLUS
  grad[idx_UZ] <- TMP_ZERO
  grad[idx_UP] <- TMP_PLUS
  
  return(grad)
}

gradLogPrior <- function(X, mu, Q) {
  #' Assumes parameters are follows MVNorm(mu, Q^{-1})
  #' 
  #' @param X vector containing \beta_0's, \beta_+'s, and X_i's
  #' @param mu mean vector
  #' @param Q precision matrix / inverse covariance Sigma^{-1}
  #' 
  Xminmu <- Matrix(X - mu)
  return(-Q %*% Xminmu)
}

hessLogLik <- function(Y, Z, X) {
  idx_BZ <- 1:ncol(Z)
  idx_BP <- (ncol(Z)+1):(ncol(Z)*2)
  idx_UZ <- (2*ncol(Z)+1):(2*ncol(Z)+nrow(Z))
  idx_UP <- (2*ncol(Z)+nrow(Z)+1):length(X)
  B_ZERO <- X[idx_BZ]
  B_PLUS <- X[idx_BP]
  U_ZERO <- X[idx_UZ]
  U_PLUS <- X[idx_UP]
  ZBU_ZERO <- Z %*% B_ZERO + U_ZERO
  ZBU_PLUS <- Z %*% B_PLUS + U_PLUS
  
  # common intermediate terms
  TMP_ZERO <- -1 / (exp(ZBU_ZERO) + exp(ZBU_PLUS) + 2)
  # Coerce TMP_PLUS to a dense matrix so we can perform column-wise multiplication with Z later
  TMP_PLUS <- Matrix(
    (Y != 0) * (
      (exp(ZBU_PLUS) * ((exp(ZBU_PLUS) - 1) * exp(exp(ZBU_PLUS)) + 1)) / (exp(exp(ZBU_PLUS)) - 1)**2
      - exp(ZBU_PLUS) 
    ), 
    sparse = FALSE
  )
  
  hess <- Matrix(0, nrow = length(X), ncol = length(X), sparse = TRUE)
  # hessian for B_ZERO
  hess[idx_BZ, idx_BZ] <- t(Z) %*% (Z * TMP_ZERO)
  # hessian for B_PLUS
  hess[idx_BP, idx_BP] <- t(Z) %*% (Z * TMP_PLUS)
  # hessian for B_ZERO w/ U_ZERO
  hess[idx_BZ, idx_UZ] <- t(Z * TMP_ZERO)
  # hessian for B_PLUS w/ U_PLUS
  hess[idx_BP, idx_UP] <- t(Z * TMP_PLUS)
  # hessian for U_ZERO (diagonal terms only)
  hess[idx_UZ, idx_UZ] <- Diagonal(x = as.numeric(TMP_ZERO))
  # hessian for U_PLUS (diagonal terms only)
  hess[idx_UP, idx_UP] <- Diagonal(x = as.numeric(TMP_PLUS))
  
  return(forceSymmetric(hess))
}

hessLogPrior <- function(Q) {
  return(-Q)
}

logLik <- function(Y, Z, X) {
  B_ZERO <- X[1:ncol(Z)]
  B_PLUS <- X[(ncol(Z)+1):(ncol(Z)*2)]
  U_ZERO <- X[(2*ncol(Z)+1):(2*ncol(Z)+nrow(Z))]
  U_PLUS <- X[(2*ncol(Z)+nrow(Z)+1):length(X)]
  ZBU_ZERO <- Z %*% B_ZERO + U_ZERO
  ZBU_PLUS <- Z %*% B_PLUS + U_PLUS
  
  # up to proportionality; sum(log(Y_i)) is independent of X
  lls <- (Y == 0) * (-log(1 + exp(ZBU_ZERO))) + 
    (Y != 0) * (
      Y*ZBU_PLUS - 
      log(1 + exp(-ZBU_ZERO)) - 
      exp(ZBU_PLUS) - 
      log(1 - exp(-exp(ZBU_PLUS)))
    )
  
  return(sum(lls))
}

logPrior <- function(X, mu, Q) {
  k <- length(X)
  Xminmu <- Matrix(X - mu)
  L <- Matrix::expand(Cholesky(Q, perm = TRUE))$L
  halflogdetQ <- sum(log(diag(L)))
  out <- as.numeric(
    halflogdetQ - (k / 2) * log(2 * pi) - .5 * sum((Xminmu) * (Q %*% Xminmu))
  )
  return(out)
}

margPost <- function(Y, Z, X, par, G, mu.theta = FALSE, Q.theta = FALSE, verbose = FALSE, acc = 1e-7, return.X = FALSE) {
  print(par)
  par_ZERO <- par[1:2]
  par_PLUS <- par[3:4]
 
  Q.U_ZERO <- getPrecMtrx(par_ZERO, G)
  Q.U_PLUS <- getPrecMtrx(par_PLUS, G)
  Q.BETA <- Diagonal(2 * ncol(Z), 1e-6)
  
  Q <- rbind(
    cbind(Q.BETA, Matrix(0, nrow(Q.BETA), ncol(Q.U_ZERO) + ncol(Q.U_PLUS))),
    cbind(Matrix(0, nrow(Q.U_ZERO), ncol(Q.BETA)), Q.U_ZERO, Matrix(0, nrow(Q.U_ZERO), ncol(Q.U_PLUS))),
    cbind(Matrix(0, nrow(Q.U_PLUS), ncol(Q.BETA) + ncol(Q.U_ZERO)), Q.U_PLUS)
  )
  mu <- rep(0, length(X))
  obj.prop <- obj.curr <- -(logLik(Y, Z, X) + logPrior(X, mu, Q))
  scaled.acc <- max(abs(obj.curr), 1) * acc
  tol.test <- Inf
  
  # Newton-Raphson optimisation to find X that
  # maximises log-likelihood of the full conditional
  # (i.e. find the roots of the full conditional gradient)
  while(TRUE) {
    alpha <- 1
    X.prop <- X
    obj.prop <- obj.curr
    grad.fc <- gradLogLik(Y, Z, X) + gradLogPrior(X, mu, Q)
    neg.hess.fc <- -hessLogLik(Y, Z, X) - hessLogPrior(Q)
    diff <- Matrix::solve(neg.hess.fc, grad.fc)
    X <- X.prop + alpha * diff
    
    obj.curr <- -(logLik(Y, Z, X) + logPrior(X, mu, Q))
  
    if (verbose) {
      plot(X, ylim = c(-3, 3))
      print(paste("obj.curr = ", obj.curr))
      print(paste("obj.curr > obj.prop = ", (obj.curr > obj.prop)))
    }
    
    while (obj.curr > obj.prop) {
      alpha <- alpha / 2
      X <- X.prop + alpha * diff
      obj.curr <- -(logLik(Y, Z, X) + logPrior(X, mu, Q))
      
      if (verbose) {
        print(paste("alpha = ", alpha, ", obj.curr - obj.prop = ", round(obj.curr - obj.prop, 2)))
        print(paste("obj.curr = ", obj.curr))
      }
      
      if (obj.curr < obj.prop && verbose) {
        print("REDUCED STEP SIZE REACHED.")
      }
      
      if (alpha <= 2**(-9)) {
        X <- X.prop
        obj.curr <- -(logLik(Y, Z, X) + logPrior(X, mu, Q))
        
        if (verbose) {
          print("STOP MOVING.")
        }
        break
      }
    }
    
    scaled.acc <- max(abs(obj.curr), 1) * acc
    tol.test <- sqrt(sum(grad.fc**2))
    
    if (verbose) {
      # print(paste("obj.curr = ", obj.curr, ", obj.prop = ", obj.prop))
      # print(paste("tol.test = ", round(tol.test, 3), ", scaled.acc = ", round(scaled.acc, 3)))
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
    return(list(X_hat = X, Q_hat = neg.hess.fc))
  } else {
    # logPrior(X, X, neg.hess.fc) is the Laplace approximation of
    # the full conditional evaluated at X_hat(theta). 
    # We also assume that the log prior of theta is 0.
    if (!mu.theta & !Q.theta) {
      obj <- as.numeric(
        -(0 + logPrior(X, mu, Q) + logLik(Y, Z, X) - logPrior(X, X, neg.hess.fc))
      )
    # Here we assume that the theta is a priori normally distributed with params (mu.theta, Q.theta^-1).
    } else {
      obj <- as.numeric(
        -(logPrior(par, mu.theta, Q.theta) + logPrior(X, mu, Q) + logLik(Y, Z, X) - logPrior(X, X, neg.hess.fc))
      )
    }

    counter <<- counter + 1
    print(paste("counter = ", counter))
    return(obj)
  }
}

plotMap <- function(raster, colours, plot.title, labels = NULL, bg.colour = "skyblue", padding = 0) {
  map.theme <- rasterTheme(region = colours, panel.background = list(col = bg.colour))
  map.theme$layout.heights[
    c(
      'bottom.padding', 
      'top.padding', 
      'key.sub.padding',
      'axis.xlab.padding',
      'key.axis.padding',
      'main.key.padding'
    )
  ] <- padding
  
  if (!is.null(labels)) {
    brks <- seq(0, length(labels), length.out = length(labels))
    l <- levelplot(
      raster,
      par.settings = map.theme,
      at = labels,
      colorkey = list(at = brks, labels = list(at = brks, labels = labels)),
      margin = FALSE,
    )
  } else {
    l <- levelplot(
      raster,
      par.settings = map.theme,
      margin = FALSE
    )
  }
  l$aspect.fill <- TRUE
  l$xlab <- NULL
  l$ylab <- NULL
  l$main <- plot.title
  return(l)
}

plotMapFromDataFrame <- function(coords, data, colours, plot.title, labels = NULL, bg.colour = "skyblue", padding = 0.5) {
  xyz <- coords
  xyz$z <- data
  xyz.raster <- rasterFromXYZ(xyz, crs = proj.str)
  plotMap(xyz.raster, colours, plot.title, labels = labels, bg.colour = bg.colour, padding = padding)
}

plotPanels <- function(coords, result, Y, Z) {
  covar.names <- Z@Dimnames[[2]][2:(ncol(Z) - 1)]
  B_ZERO_hat <- result$X_hat[1:ncol(Z)]
  B_PLUS_hat <- result$X_hat[(ncol(Z)+1):(2*ncol(Z))]
  U_ZERO_hat <- result$X_hat[(2*ncol(Z)+1):(2*ncol(Z)+nrow(Z))]
  U_PLUS_hat <- result$X_hat[(2*ncol(Z)+nrow(Z)+1):length(result$X_hat)]
  ZBU_ZERO_hat <- Z %*% B_ZERO_hat + U_ZERO_hat
  ZBU_PLUS_hat <- Z %*% B_PLUS_hat + U_PLUS_hat
  PR_ZERO_hat <- 1 / (1 + exp(ZBU_ZERO_hat))
  RT_PLUS_hat <- exp(ZBU_PLUS_hat)
  E_Y_hat <- (1 -  PR_ZERO_hat) * RT_PLUS_hat / (1 - exp(-RT_PLUS_hat))
  
  colours <- rev(brewer.pal(9, "RdYlGn"))
  labels <- c(0, 1, 2, 5, 10, 15, 20, 25, 31)
  
  plots <- list()
  
  # actual Y values
  plot.title <- TeX("Actual counts")
  plots[[1]] <- plotMapFromDataFrame(coords, Y, colours, plot.title, labels = labels)
  # E(Y_i)
  plot.title <- TeX("Expected counts")
  plots[[2]] <- plotMapFromDataFrame(df.coords, E_Y_hat, colours, plot.title, labels = labels)
  # U_ZERO_HAT
  plot.title <- TeX("Latent spatial effects $U_0$")
  plots[[3]] <- plotMapFromDataFrame(coords, U_ZERO_hat, colours, plot.title)
  # U_PLUS_HAT
  plot.title <- TeX("Latent spatial effects $U_+$")
  plots[[4]] <- plotMapFromDataFrame(coords, U_PLUS_hat, colours, plot.title)
  # P(Y_i > 0)
  plot.title <- TeX("$\\pi(X_i^T\\beta_0 + U_0)")
  plots[[5]] <- plotMapFromDataFrame(coords, 1 - PR_ZERO_hat, colours, plot.title)
  # Rate parameter for Y_i > 0
  plot.title <- TeX("$\\lambda(X_i^T\\beta_+ + U_+)")
  plots[[6]] <- plotMapFromDataFrame(coords, RT_PLUS_hat, colours, plot.title)
  
  for (name in covar.names) {
    plot.title <- TeX(name)
    plots[[(length(plots)+1)]] <- plotMapFromDataFrame(coords, Z[, name], colours, plot.title)
  }
  
  do.call("grid.arrange", c(plots, list(left = "Latitude", bottom = "Longitude")))
}

runOptimProcedure <- function(theta_init, X_init, Y, Z, G, acc = 1e-4, reltol = 1e-4, itnmax = 1000) {
  opt <- optimx(
    theta_init,
    margPost, 
    Y = Y,
    Z = Z, 
    X = X_init, 
    G = G,
    acc = acc,
    method = "Nelder-Mead",
    itnmax = itnmax, 
    control = list(kkt = FALSE, reltol = reltol)
  )
  theta_hat <- c(opt$p1, opt$p2, opt$p3, opt$p4)
  XQ_hat <- margPost(Y, Z, X_init, theta_hat, G, acc = acc, return.X = TRUE, verbose = TRUE)
  X_hat <- XQ_hat$X_hat
  Q_hat <- XQ_hat$Q_hat
  return(list(opt = opt, theta_hat = theta_hat, X_hat = X_hat, Q_hat = Q_hat))
}

##### MAIN #####

setwd(working.dir)

df <- read.csv(str_glue("{data.dir}/csv/{fname}"))
df.coords <- df[1:2]
df.covars <- df[c(-1,-2)]
df.covars$range.humidity <- df.covars$temp.range * df.covars$rel.humidity
df.covars$temp.humidity <- df.covars$avg.temp * df.covars$rel.humidity

G <- getLaplMtrx(df.coords, res, verbose = TRUE)
Y <- Matrix(df.covars$count)
Z.dtr <- Matrix(
  cbind(rep(1, length(Y)), as.matrix(subset(df.covars, select = -c(count, avg.temp, precipitation, temp.humidity))))
)
Z.tmp <- Matrix(
  cbind(rep(1, length(Y)), as.matrix(subset(df.covars, select = -c(count, temp.range, precipitation, range.humidity))))
)

# Initial parameters
theta <- c(0, 0, 0, 0)
X <- rep(0, 2 * (ncol(Z.dtr) + nrow(Z.dtr)))

counter <- 0
set.seed(seed)
result.dtr <- runOptimProcedure(theta, X, Y, Z.dtr, G)

counter <- 0
set.seed(seed)
result.tmp <- runOptimProcedure(theta, X, Y, Z.tmp, G)

# Should we use tmp or dtr?
log.lik.dtr <- -margPost(Y, Z.dtr, X, result.dtr$theta_hat, G, acc = 1e-4, verbose = TRUE) # 2492.772
log.lik.tmp <- -margPost(Y, Z.tmp, X, result.tmp$theta_hat, G, acc = 1e-4, verbose = TRUE) # 2493.009

# Confidence intervals
conf.int.dtr <- getConfInt(result.dtr)
conf.int.tmp <- getConfInt(result.tmp)

##### PLOTTING THE SPATIAL EFFECTS #####

plotPanels(df.coords, result.dtr, Y, Z.dtr)

plot.U_ZERO_SIG.dtr <- plotMapFromDataFrame(
  df.coords, 
  conf.int.dtr$significant[(2*ncol(Z.tmp)+1):(2*ncol(Z.tmp)+nrow(Z.tmp))], 
  c("#1A9850", "#D73027"),
  TeX("U_0 significantly different from 0? (dtr)"),
  labels = c(FALSE, 0.5, TRUE)
)

plot.U_PLUS_SIG.dtr <- plotMapFromDataFrame(
  df.coords, 
  conf.int.dtr$significant[(2*ncol(Z.tmp)+nrow(Z.tmp)+1):length(result.dtr$X_hat)], 
  c("#1A9850", "#D73027"),
  TeX("U_+ significantly different from 0? (dtr)"),
  labels = c(FALSE, 0.5, TRUE)
)
