# CONFIGURATION VARIABLES -------------------------------------------------

if (getOption('run.main', default = TRUE)) {
  options(run.model_checking = TRUE)
  source("~/Documents/projects/diss/config.R")
}

# IMPORTS -----------------------------------------------------------------

library(countreg)
library(extraDistr)
library(gridExtra)
library(latex2exp)
library(Matrix)
library(pROC)
library(raster)
library(rasterVis)
library(RColorBrewer)

# HELPER FUNCTIONS --------------------------------------------------------

getConfInt <- function(X, Z, Q, alpha = 0.05) {
  Sigma <- Matrix::solve(Q, Matrix::Diagonal(nrow(Q)))
  conf.int <- data.frame(
    par = c(
      paste(Z@Dimnames[[2]], "_0", sep = ""),
      paste(Z@Dimnames[[2]], "_+", sep = ""),
      rep("U_ZERO", length(X$U_ZERO)),
      rep("U_PLUS", length(X$U_PLUS))
    ),
    mean = as.numeric(unlist(X)),
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

getRasterFromDataFrame <- function(coords, data, crs = proj.str) {
  xyz <- coords
  xyz$z <- data
  return(rasterFromXYZ(xyz, crs = proj.str))
}

getRates <- function(X, Z) {
  ZBU_PLUS <- Z %*% X$B_PLUS + X$U_PLUS
  return(exp(ZBU_PLUS))
}

plotConfInts <- function(coords, interval, plot.title, bg.colour = "skyblue", legend = TRUE) {
  r <- getRasterFromDataFrame(coords, interval)
  r <- ratify(r)
  rat <- data.frame(ID = c(0, 1), significance = c("FALSE", "TRUE"))
  levels(r) <- rat
  colours <- c("#1A9850", "#D73027")
  map.theme <- rasterTheme(region = colours, panel.background = list(col = bg.colour))
  levelplot(
    r, 
    par.settings = map.theme,
    att = "significance", 
    main = plot.title,
    colorkey = legend,
    xlab = "Longitude",
    ylab = "Latitude"
  )
}

plotMap <- function(raster, colours, plot.title, layers, labels = NULL, bg.colour = "skyblue", padding = 0) {
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
  return(l + layers)
}

plotMapFromDataFrame <- function(coords, data, colours, plot.title, layers, labels = NULL, bg.colour = "white", padding = 0.5) {
  xyz.raster <- getRasterFromDataFrame(coords, data, proj.str)
  plotMap(xyz.raster, colours, plot.title, layers, labels = labels, bg.colour = bg.colour, padding = padding)
}

plotPanels <- function(coords, X, Y, Z, layers) {
  PR_FIRE <- getProbs(X, Z)
  RT_FIRE <- getRates(X, Z)
  pred <-  PR_FIRE * RT_FIRE / (1 - exp(-RT_FIRE))
  
  colours <- brewer.pal(9, "Greys")
  labels <- c(0, 1, 2, 5, 10, 15, 20, 25, 31)
  
  plots <- list()
  
  # actual Y values
  plot.title <- TeX("Actual counts")
  plots[[1]] <- plotMapFromDataFrame(coords, Y, colours, plot.title, layers, labels = labels)
  # E(Y_i)
  plot.title <- TeX("Expected counts")
  plots[[2]] <- plotMapFromDataFrame(coords, pred, colours, plot.title, layers, labels = labels)
  # U_ZERO_HAT
  plot.title <- TeX("Latent spatial effects $U_0$")
  plots[[3]] <- plotMapFromDataFrame(coords, X$U_ZERO, colours, plot.title, layers)
  # U_PLUS_HAT
  plot.title <- TeX("Latent spatial effects $U_+$")
  plots[[4]] <- plotMapFromDataFrame(coords, X$U_PLUS, colours, plot.title, layers)
  # P(Y_i > 0)
  plot.title <- TeX("$\\pi(z_i^T\\beta_0 + U_0)")
  plots[[5]] <- plotMapFromDataFrame(coords, PR_FIRE, colours, plot.title, layers)
  # Rate parameter for Y_i > 0
  plot.title <- TeX("$\\lambda(z_i^T\\beta_+ + U_+)")
  plots[[6]] <- plotMapFromDataFrame(coords, RT_FIRE, colours, plot.title, layers)
  
  do.call("grid.arrange", c(plots, list(left = "Latitude", bottom = "Longitude")))
}

plotRootogram <- function(X, Y, Z, main = "Poisson hurdle model with latent spatial effects") {
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
  abline(h = c(-1, 1), lty = 2, lwd = 2, col = "blue")
  legend("topright", legend = "Tukey warning limits", lty = 2, col = "blue", box.lty = 0, cex = 0.9)
}

splitParams <- function(X, Z) {
    return(list(
      B_ZERO = X[1:ncol(Z)], 
      B_PLUS = X[(ncol(Z)+1):(2*ncol(Z))], 
      U_ZERO = X[(2*ncol(Z)+1):(2*ncol(Z)+nrow(Z))], 
      U_PLUS = X[(2*ncol(Z)+nrow(Z)+1):length(X)]
    ))
}


if (getOption('run.model_checking', default = FALSE)) {
  load(str_glue("{data.dir}/RData/{fname}"))
  Z@Dimnames[[2]][1] <- "(Intercept)"
  Z@Dimnames[[2]][6] <- "skt:rhm"
  X <- splitParams(result$X_hat, Z)
  
  # Assess GOF for "hurdle" part using ROC curve
  obs.binary <- Y > 0
  pred.binary <- getProbs(X, Z)
  pROC.obj <- roc(as.numeric(obs.binary), as.numeric(pred.binary), ci = TRUE)
  
  # Assess GOF of "count" part using (truncated) Poisson residuals
  obs.count <- Y[Y > 0]
  rate.count <- getRates(X, Z)[Y > 0]
  pred.count <- rate.count / (1 - exp(-rate.count))
  var.pred <- pred.count - exp(-rate.count) * pred.count**2
  resids <- (obs.count - pred.count) / sqrt(var.pred)
  
  # Construct confidence intervals for parameters
  conf.int <- getConfInt(X, Z, result$Q_hat)
}

# PLOTS -------------------------------------------------------------------

if (getOption('run.model_checking', default = FALSE)) {
  cpoly <- getData("GADM", path = str_glue("{data.dir}/rds"), country = cname, level = 1)
  layers <- layer(sp.lines(cpoly))
  plotPanels(coords, X, Y, Z, layers)
  plotRootogram(X, Y, Z)
  
  par(mfrow = c(1, 2), pty = "s")
  plot.roc(
    pROC.obj,
    asp = 1,
    auc.polygon = TRUE,
    max.auc.polygon = TRUE,
    grid = TRUE,
    print.auc = TRUE
  )
  par(pty = "m")
  plot(pred.count, resids, xlab = "Expected count", ylab = "Standardised residuals")
  abline(h = 0)
  par(mfrow = c(1, 1))
  
  # CI plots
  plot1 <- plotConfInts(
    coords,
    conf.int$significant[conf.int$par == "U_ZERO"],
    TeX("$U_0$ / $U_+$ significantly different from zero?"),
  )
  plot2 <- plotConfInts(
    coords,
    conf.int$significant[conf.int$par == "U_PLUS"],
    TeX("$U_+$ significantly different from zero?"),
  )
  c(plot1, plot2, y.same = TRUE, x.same = FALSE)
  options(run.model_checking = FALSE)
}