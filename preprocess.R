# CONFIGURATION VARIABLES -------------------------------------------------

if (getOption("run.main", default = TRUE)) {
  options(run.preprocess = TRUE)
  source("~/Documents/diss/config.R")
}

# IMPORTS -----------------------------------------------------------------

library(elevatr)
library(raster)
library(rgdal)
library(rgeos)

# HELPER FUNCTIONS --------------------------------------------------------

flattenRaster <- function(raster, map, fun) {
  #'
  #'  Flattens a RasterBrick object into a RasterLayer object that contains counts 
  #'  # of positive values at each grid.
  #'  
  
  args <- as.list(raster)
  # Count no. of non-zero values at each grid
  args$fun <- fun
  # flattened <- mask(do.call(mosaic, args), map)
  flattened <- do.call(mosaic, args)
  names(flattened) <- "count"
  return(flattened)
}

getSmallPolys <- function(poly, minarea=0.01) {
  #' 
  #' Lossfully compress polygon size by removing small polygons for faster processing times
  #' Code from https://gis.stackexchange.com/a/62405/155373 
  #' 
  
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

getSVPressure <- function(temp, is.water = TRUE) {
  #'
  #' Calculate saturation vapour pressure using Teten's formula.
  #' See Eqn (7.5) in https://www.ecmwf.int/en/elibrary/16648-part-iv-physical-processes 
  #' for more details.
  #' 
  
  a.1 <- 611.21
  a.3 <- ifelse(is.water, 17.502, 22.587)
  a.4 <- ifelse(is.water, 32.19, -0.7)
  T_0 <- 273.16
  
  return(a.1 * exp(a.3 * ((temp - T_0)/(temp - a.4))))
}

loadCovarRaster <- function(fname, vname, raster) {
  covar <- raster(fname, varname = vname)
  covar <- resample(covar, raster)
  names(covar) <- vname
  return(covar)
}

# MAIN --------------------------------------------------------------------

if (getOption("run.preprocess", default = FALSE)) {
  cpoly <- getData("GADM", path = str_glue("{data.dir}/rds"), country = cname, level = 1)
  
  fire <- raster::brick(fname.fire, varname = vname.fire)
  if (is.lowres) {
    fire <- aggregate(fire, fact = lowres.factor, fun = sum)
  }
  fire <- flattenRaster(fire, cpoly, function(x, na.rm) {sum(x > 0, na.rm = na.rm)})
  
  airt <- loadCovarRaster(fname.airt, vname.airt, fire)
  dewp <- loadCovarRaster(fname.dewp, vname.dewp, fire)
  temp <- loadCovarRaster(fname.temp, vname.temp, fire)
  elev <- loadCovarRaster(fname.elev, vname.elev, fire)
  vegc <- loadCovarRaster(fname.vegc, vname.vegc, fire)
  
  # coerce to data.frame
  data.raster <- mask(stack(fire, airt, dewp, temp, elev, vegc), cpoly)
  df <- as.data.frame(data.raster, xy = TRUE, na.rm = TRUE)
  df$x <- round(df$x, 2)
  df$y <- round(df$y, 2)
  df$ptc <- df$ptc / 100
  df$skt <- df$skt - 273.15
  df$rhm <- getSVPressure(df$d2m) / getSVPressure(df$t2m)
  df <- subset(df, select = -c(d2m, t2m))
  
  write.csv(df, str_glue("{data.dir}/csv/{csv.name}"), row.names = FALSE)
  options(run.preprocess = FALSE)
}

