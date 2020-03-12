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

# MAIN --------------------------------------------------------------------

if (getOption("run.preprocess", default = FALSE)) {
  cpoly <- getData("GADM", path = str_glue("{data.dir}/rds"), country = cname, level = 1)
  
  c.var <- raster::brick(fname.fire, varname = vname.fire)
  if (is.lowres) {
    c.var <- aggregate(c.var, fact = lowres.factor, fun = mean)
  }
  c.var.flat <- flattenRaster(c.var, cpoly, function(x, na.rm) {sum(x > 0, na.rm = na.rm)})
  
  elev <- get_elev_raster(c.var.flat, src = "aws", z = 6)
  elev.resampled <- resample(elev, c.var.flat)
  names(elev.resampled) <- "elevation"
  
  dewp <- raster::brick(fname.era5, varname = vname.dewp)[[2]]
  dewp.resampled <- resample(dewp, c.var.flat)
  names(dewp.resampled) <- "dewpoint.temp"
  
  airt <- raster::brick(fname.airt, varname = vname.airt)[[2]]
  airt.resampled <- resample(airt, c.var.flat)
  names(airt.resampled) <- "air.temp.2m"
  
  dtrg <- raster::brick(fname.dtrg, varname = vname.dtrg)[[12 * (year - 1) + month]]
  dtrg.resampled <- resample(dtrg, c.var.flat)
  names(dtrg.resampled) <- "temp.range"
  
  temp <- raster::brick(fname.temp)[[month]]
  temp.resampled <- resample(temp, c.var.flat)
  names(temp.resampled) <- "avg.temp"
  
  prec <- raster::brick(fname.prec, varname = vname.prec)[[12 * (year - 1) + month]]
  prec.resampled <- resample(prec, c.var.flat)
  names(prec.resampled) <- "precipitation"
  
  # coerce to data.frame
  data.raster <- mask(
    stack(
      c.var.flat,
      elev.resampled,
      prec.resampled,
      temp.resampled,
      airt.resampled,
      dewp.resampled,
      dtrg.resampled
    ), 
    cpoly
  )
  
  df <- as.data.frame(data.raster, xy = TRUE, na.rm = TRUE)
  df$x <- round(df$x, 2)
  df$y <- round(df$y, 2)
  df$elevation <- pmax(df$elevation, 0)
  df$avg.temp <- df$avg.temp - 273.15
  df$rel.humidity <- getSVPressure(df$dewpoint.temp) / getSVPressure(df$air.temp.2m)
  df <- subset(df, select = -c(dewpoint.temp, air.temp.2m))
  
  write.csv(df, str_glue("{data.dir}/csv/{csv.name}"), row.names = FALSE)
}
