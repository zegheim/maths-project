# IMPORTS -----------------------------------------------------------------

library(gridExtra)
library(latex2exp)
library(raster)
library(rasterVis)
library(RColorBrewer)
library(stringr)

# CONFIGURATION VARIABLES -------------------------------------------------

fname <- "result.20200311163133.RData"
proj.str <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
working.dir <- "~/Documents/diss"

data.dir <- str_glue("{working.dir}/data")

# DO NOT EDIT BELOW THIS LINE ---------------------------------------------

# HELPER FUNCTIONS --------------------------------------------------------

getConfInt <- function(result, alpha = 0.05) {
  Sigma_hat <- Matrix::solve(result$Q_hat, Matrix::Diagonal(nrow(result$Q_hat)))
  conf.int <- data.frame(
    lwr = as.numeric(result$X_hat) + qnorm(alpha/2) * sqrt(diag(Sigma_hat)),
    upr = as.numeric(result$X_hat) + qnorm(1 - alpha/2) * sqrt(diag(Sigma_hat))
  )
  conf.int$significant <- !(conf.int$lwr <= 0 & conf.int$upr >= 0)
  return(conf.int)
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
  plots[[2]] <- plotMapFromDataFrame(coords, E_Y_hat, colours, plot.title, labels = labels)
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

# MAIN --------------------------------------------------------------------

setwd(working.dir)
load(str_glue("{data.dir}/RData/{fname}"))


# PLOTS -------------------------------------------------------------------

plotPanels(coords, result, Y, Z)
