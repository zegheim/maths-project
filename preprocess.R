###########################
# CONFIGURATION VARIABLES #
###########################

data.dir <- "~/Documents/diss/data"
working.dir <- "~/Documents/diss/"

cname <- "Indonesia"
fname <- "/media/zegheim/Justin_SSD/nc_ina/cams_gfas_ga_1507.nc"
vname <- "frpfire"

plot.title <- "No. of wildfire occurrences in Jul 2015"

proj.str <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

### DO NOT EDIT BELOW THIS LINE ###

###########
# IMPORTS #
###########

library(boot)
library(elevatr)
library(MASS)
library(optimx)
library(pscl)
library(raster)
library(rasterVis)
library(RColorBrewer)
library(rgdal)
library(rgeos)

########################
# FUNCTION DEFINITIONS #
########################

# Flattens a RasterBrick object into a RasterLayer object that contains counts # of positive values at each grid
flattenRaster <- function(raster, map, fun) {
  args <- as.list(raster)
  # Count no. of non-zero values at each grid
  args$fun <- fun
  # flattened <- mask(do.call(mosaic, args), map)
  flattened <- do.call(mosaic, args)
  names(flattened) <- "count"
  return(flattened)
}

# Lossfully compress polygon size by removing small polygons for faster processing times
# Code from https://gis.stackexchange.com/a/62405/155373 
getSmallPolys <- function(poly, minarea=0.01) {
  # Get the areas
  areas <- lapply(poly@polygons, function(x) sapply(x@Polygons, function(y) y@area))
  # Which are the big polygons?
  bigpolys <- lapply(areas, function(x) which(x > minarea))
  # Get only the big polygons and extract them
  for(i in 1:length(bigpolys)){
    if(length(bigpolys[[i]]) >= 1 && bigpolys[[i]] >= 1){
      poly@polygons[[i]]@Polygons <- poly@polygons[[i]]@Polygons[bigpolys[[i]]]
      poly@polygons[[i]]@plotOrder <- 1:length(poly@polygons[[i]]@Polygons)
    }
  }
  return(poly)
}

# Helper function to count how many NA values in a RasterLayer
countMissing <- function(raster) {
  sum(is.na(getValues(raster)))
}

#################
# PREPROCESSING #
#################

setwd(working.dir)

cpoly <- getData("GADM", download = FALSE, path = data.dir, country = cname, level = 1)
# Simplify the geometry of the polygon shapes
cpoly.simplified <- gSimplify(getSmallPolys(cpoly), tol = 0.01, topologyPreserve = TRUE)

c.var <- raster::brick(fname, varname = vname)
c.var.flat <- flattenRaster(c.var, cpoly, function(x, na.rm) {sum(x > 0, na.rm = na.rm)})
c.var.flat.lowres <- aggregate(c.var.flat, fact = 5, fun = sum)

c.df <- as.data.frame(c.var.flat.lowres, xy = TRUE, na.rm = TRUE)

# get elevation data
elev <- get_elev_raster(c.var.flat, src = "aws", z = 6)
elev.cropped <- resample(crop(elev, c.var.flat), c.var.flat, method="bilinear")
elev.cropped.lowres <- aggregate(elev.cropped, fact = 5, fun = mean)
names(elev.cropped.lowres) <- "elevation"

# get temperature data for July
temp <- raster::brick("/media/zegheim/Justin_SSD/cpc_global_temp/tmax.2015.nc")
temp.1507 <- temp[[182:212]]
temp.1507.flattened <- crop(flattenRaster(temp.1507, cpoly, mean), c.var.flat)
temp.1507.resampled <- aggregate(resample(temp.1507.flattened, c.var.flat, method = "bilinear"), fact = 5, fun = mean)
names(temp.1507.resampled) <- "max.temp"

# coerce to data.frame
df <-  as.data.frame(mask(stack(c.var.flat.lowres, elev.cropped.lowres, temp.1507.resampled), cpoly), xy = TRUE, na.rm = TRUE)

#############
# MODELLING #
#############

set.seed(17071996L)
Y <- df$count
# Normalise all covariates to have mean 0 and sd 1
X_scaled <- scale(as.matrix(subset(df, select = -c(count))))

# Add intercept column
X <- cbind(rep(1, length(Y)), X_scaled)
colnames(X)[1] <- "intercept"

loglik <- function(par, X, Y) {
  XB_PLUS <- X %*% par[1:5]
  XB_ZERO <- X %*% par[6:10]
  sum(
    ifelse(Y == 0, log(1 - inv.logit(XB_ZERO)), 0) +
    ifelse(Y != 0, log(inv.logit(XB_ZERO)) - exp(XB_PLUS) + (Y*XB_PLUS) - lgamma(Y + 1) - log(1 - exp(-exp(XB_PLUS))), 0)
  )                                    
}

par_init <- as.numeric(c(
    glm.fit(X, Y, family = poisson())$coefficients,
    glm.fit(X, factor(Y>0), family = binomial(link = "logit"))$coefficients
))

opt <- optimx(par_init, loglik, Y = Y, X = X, method = "BFGS", itnmax = 10000, control=list(maximize = TRUE))
B_PLUS_hat <- as.numeric(tail(opt, 1)[1:10])[1:5]
B_ZERO_hat <- as.numeric(tail(opt, 1)[1:10])[6:10]

# Check with using library
fire.hurdle <- hurdle(count ~ scale(x) + scale(y) + scale(elevation) + scale(max.temp), data = df, separate = FALSE)

# Compare final value
opt$value; fire.hurdle$optim$value

# Compare coefficients
as.numeric(fire.hurdle$coefficients$count); B_PLUS_hat
as.numeric(fire.hurdle$coefficients$zero); B_ZERO_hat

inf_matrix <- -1 * optimHess(par = as.numeric(opt[1:10]), fn = loglik, Y = Y, X = X, control=list(fnscale = -1))
inf_matrix[1:5, 6:10] <- 0
inf_matrix[6:10, 1:5] <- 0
inv.inf_matrix <- solve(inf_matrix)

# generating samples of MLE
mles <- mvrnorm(n = 100, mu = as.numeric(tail(opt, 1)[1:10]), Sigma = inv.inf_matrix)

