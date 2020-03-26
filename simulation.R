# CONFIGURATION VARIABLES -------------------------------------------------

if (getOption("run.main", default = TRUE)) {
  options(run.simulation = TRUE)
  source("config.R")
}

# IMPORTS -----------------------------------------------------------------

library(extraDistr)
library(gridExtra)
library(latex2exp)
sourceFunctions("model_checking.R")
sourceFunctions("model_fitting.R")

# HELPER FUNCTIONS --------------------------------------------------------

# MAIN --------------------------------------------------------------------

if (getOption("run.simulation", default = FALSE)) {
  load(str_glue("{data.dir}/RData/{fname}"))
  X <- splitParams(result$X_hat, Z)
  
  # Simulate from fitted model
  set.seed(seed)
  Y.sim <- sapply(getProbs(X, Z), function(p) rbinom(1, 1, p)) * sapply(getRates(X, Z), function(r) rtpois(1, r))
  
  # Compare actual data to simulated data
  plot.actual <- plotMapFromDataFrame(
    coords, 
    Y, 
    rev(brewer.pal(9, "RdYlGn")), 
    "Actual counts",
    labels = c(0, 1, 2, 5, 10, 15, 20, 25, 31)
  )
  
  plot.simulated <- plotMapFromDataFrame(
    coords, 
    Y.sim, 
    rev(brewer.pal(9, "RdYlGn")), 
    "Simulated counts",
    labels = c(0, 1, 2, 5, 10, 15, 20, 25, 31)
  )
  
  grid.arrange(plot.actual, plot.simulated, ncol = 2)
  
  # Initialise parameters
  theta.init <- c(0, 0, 0, 0)
  X.init <- rep(0, 2 * (ncol(Z) + nrow(Z)))
  
  # Fit model to simulated data
  set.seed(seed)
  counter <- 0
  result.sim <- runOptimProcedure(theta.init, X.init, Y.sim, Z, G)
  
  # Save results
  save(coords, result.sim, Y.sim, Z, G, file = str_glue("{data.dir}/RData/{RData.name}"))
  
  # Construct confidence intervals for parameters
  X.sim <- splitParams(result.sim$X_hat, Z)
  conf.int <- getConfInt(X.sim, result.sim$Q_hat) 
  
  # CI plots
  plot1 <- plotMapFromDataFrame(
    coords,
    conf.int$significant[conf.int$par == "U_ZERO"],
    c("#1A9850", "#D73027"),
    TeX("$U_0$ significantly different from zero?"),
    labels = c(0, 0.5, 1)
  )
  plot2 <- plotMapFromDataFrame(
    coords,
    conf.int$significant[conf.int$par == "U_PLUS"],
    c("#1A9850", "#D73027"),
    TeX("$U_+$ significantly different from zero?"),
    labels = c(0, 0.5, 1)
  )
  grid.arrange(plot1, plot2)

  options(run.simulation = FALSE)
}