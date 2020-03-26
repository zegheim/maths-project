# CONFIGURATION VARIABLES -------------------------------------------------

if (getOption('run.main', default = TRUE)) {
  options(run.model_checking = TRUE)
  source("~/Documents/diss/config.R")
}

# IMPORTS -----------------------------------------------------------------

library(extraDistr)
library(gridExtra)
library(latex2exp)
library(Matrix)
library(pROC)
library(raster)
library(rasterVis)
library(RColorBrewer)

# HELPER FUNCTIONS --------------------------------------------------------

getConfInt <- function(X, Q, alpha = 0.05) {
  Sigma <- Matrix::solve(Q, Matrix::Diagonal(nrow(Q)))
  conf.int <- data.frame(
    par = c(
      rep("B_ZERO", length(X$B_ZERO)),
      rep("B_PLUS", length(X$B_PLUS)),
      rep("U_ZERO", length(X$U_ZERO)),
      rep("U_PLUS", length(X$U_PLUS))
    ),
    lwr = as.numeric(unlist(X)) + qnorm(alpha/2) * sqrt(diag(Sigma)),
    upr = as.numeric(unlist(X)) + qnorm(1 - alpha/2) * sqrt(diag(Sigma))
  )
  conf.int$significant <- !(conf.int$lwr <= 0 & conf.int$upr >= 0)
  return(conf.int)
}

getProbs <- function(X, Z) {
  ZBU_ZERO <- Z %*% X$B_ZERO + X$U_ZERO
  return(1 / (1 + exp(-ZBU_ZERO)))
}

getRates <- function(X, Z) {
  ZBU_PLUS <- Z %*% X$B_PLUS + X$U_PLUS
  return(exp(ZBU_PLUS))
}

plotConfInts <- function(coords, interval) {
  colours <- c("#1A9850", "#D73027")
  labels <- c("")
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

plotPanels <- function(coords, X, Y, Z) {
  PR_FIRE <- getProbs(X, Z)
  RT_FIRE <- getRates(X, Z)
  pred <-  PR_FIRE * RT_FIRE / (1 - exp(-RT_FIRE))
  
  colours <- rev(brewer.pal(9, "RdYlGn"))
  labels <- c(0, 1, 2, 5, 10, 15, 20, 25, 31)
  
  plots <- list()
  
  # actual Y values
  plot.title <- TeX("Actual counts")
  plots[[1]] <- plotMapFromDataFrame(coords, Y, colours, plot.title, labels = labels)
  # E(Y_i)
  plot.title <- TeX("Expected counts")
  plots[[2]] <- plotMapFromDataFrame(coords, pred, colours, plot.title, labels = labels)
  # U_ZERO_HAT
  plot.title <- TeX("Latent spatial effects $U_0$")
  plots[[3]] <- plotMapFromDataFrame(coords, X$U_ZERO, colours, plot.title)
  # U_PLUS_HAT
  plot.title <- TeX("Latent spatial effects $U_+$")
  plots[[4]] <- plotMapFromDataFrame(coords, X$U_PLUS, colours, plot.title)
  # P(Y_i > 0)
  plot.title <- TeX("$\\pi(X_i^T\\beta_0 + U_0)")
  plots[[5]] <- plotMapFromDataFrame(coords, PR_FIRE, colours, plot.title)
  # Rate parameter for Y_i > 0
  plot.title <- TeX("$\\lambda(X_i^T\\beta_+ + U_+)")
  plots[[6]] <- plotMapFromDataFrame(coords, RT_FIRE, colours, plot.title)
  
  covar.names <- Z@Dimnames[[2]][2:(ncol(Z) - 1)]
  for (name in covar.names) {
    plot.title <- TeX(name)
    plots[[(length(plots)+1)]] <- plotMapFromDataFrame(coords, Z[, name], colours, plot.title)
  }
  
  do.call("grid.arrange", c(plots, list(left = "Latitude", bottom = "Longitude")))
}

plotRootogram <- function(X, Y, Z, main = "mod.spatial") {
  at <- 0:31
  probs <- as.numeric(getProbs(X, Z))
  rates <- as.numeric(getRates(X, Z))
  pred.probs <- matrix(NA, nrow = length(rates), ncol = length(at))
  pred.probs[, 1] <- 1 - as.numeric(probs)
  for (i in 2:length(at)) {
    pred.probs[, i] <- as.numeric(probs) * dtpois(at[i], as.numeric(rates), a = 0)
  }
  
  obsrvd <- as.vector(xtabs(rep(1, NROW(Y)) ~ factor(as.numeric(Y), levels = at)))
  expctd <- colSums(pred.probs)
  rootogram.default(obsrvd, expctd, breaks = -1L:31.5, xlab = "count", main = main)
}

splitParams <- function(X, Z) {
    return(list(
      B_ZERO = X[1:ncol(Z)], 
      B_PLUS = X[(ncol(Z)+1):(2*ncol(Z))], 
      U_ZERO = X[(2*ncol(Z)+1):(2*ncol(Z)+nrow(Z))], 
      U_PLUS = X[(2*ncol(Z)+nrow(Z)+1):length(X)]
    ))
}

# MAIN --------------------------------------------------------------------

if (getOption('run.model_checking', default = FALSE)) {
  load(str_glue("{data.dir}/RData/{fname}"))
  
  X <- splitParams(result$X_hat, Z)
  
  # Assess GOF for "hurdle" part using ROC curve
  obs.binary <- Y > 0
  pred.binary <- getProbs(X, Z)
  pROC.obj <- roc(
    as.numeric(obs.binary), 
    as.numeric(pred.binary),
    plot = TRUE,
    ci = TRUE, 
    auc.polygon = TRUE,
    max.auc.polygon = TRUE,
    grid = TRUE,
    print.auc = TRUE
  )
  
  # Assess GOF of "count" part using (truncated) Poisson residuals
  obs.count <- Y[Y > 0]
  rate.count <- getRates(X, Z)[Y > 0]
  pred.count <- rate.count / (1 - exp(-rate.count))
  var.pred <- (rate.count + rate.count**2) / (1 - exp(-rate.count)) - (rate.count / (1 - exp(-rate.count)))**2
  resids <- (obs.count - pred.count) / sqrt(var.pred)
  plot(pred.count, resids)
  
  # Construct confidence intervals for parameters
  conf.int <- getConfInt(X, result$Q_hat) 
}

# PLOTS -------------------------------------------------------------------

if (getOption('run.model_checking', default = FALSE)) {
  plotPanels(coords, X, Y, Z)
  plotRootogram(X, Y, Z)
  options(run.model_checking = FALSE)
}
