###########################
# CONFIGURATION VARIABLES #
###########################

cname <- "Indonesia"
csv.name <- "~/Documents/diss/data/df_ina_lowres.csv"
data.dir <- "~/Documents/diss/data"
fname.fire <- "/media/zegheim/Justin_SSD/nc_ina/gfas/cams_gfas_ga_1507.nc"
fname.temp <- "/media/zegheim/Justin_SSD/nc_ina/tair/tair_2015_cropped.nc"
month <- 7
vname <- "frpfire"
working.dir <- "~/Documents/diss/"
is.lowres <- TRUE
lowres.factor <- 5

### DO NOT EDIT BELOW THIS LINE ###

###########
# IMPORTS #
###########

library(elevatr)
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

########
# MAIN #
########

setwd(working.dir)

# simplify the geometry of the polygon shapes
cpoly <- getData("GADM", path = data.dir, country = cname, level = 1)
cpoly.simplified <- gSimplify(getSmallPolys(cpoly), tol = 0.01, topologyPreserve = TRUE)

# get fire data
c.var <- raster::brick(fname.fire, varname = vname)
if (is.lowres) {
  c.var <- aggregate(c.var, fact = lowres.factor, fun = mean)
}
c.var.flat <- flattenRaster(c.var, cpoly.simplified, function(x, na.rm) {sum(x > 0, na.rm = na.rm)})

# get elevation data
elev <- get_elev_raster(c.var.flat, src = "aws", z = 6)
elev.cropped <- resample(elev, c.var.flat, method="bilinear")
names(elev.cropped) <- "elevation"

# get temperature data
temp <- raster::brick(fname.temp)[[month]]
temp.resampled <- resample(temp, c.var.flat)
names(temp.resampled) <- "avg.temp"

# coerce to data.frame
data.raster <- mask(stack(c.var.flat, elev.cropped, temp.resampled), cpoly.simplified)
df <- as.data.frame(data.raster, xy = TRUE, na.rm = TRUE)
df$avg.temp <- df$avg.temp - 273.15
df$elevation <- pmax(df$elevation, 0)
df$x <- round(df$x, 2)
df$y <- round(df$y, 2)
write.csv(df, csv.name, row.names = FALSE)
