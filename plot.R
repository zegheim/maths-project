###########################
# CONFIGURATION VARIABLES #
###########################

data.dir <- "~/Documents/diss/data"
working.dir <- "~/Documents/diss/"

cname <- "Australia"
fname <- "/media/zegheim/Justin_SSD/nc_aus/cams_gfas_ga_0612.nc"
vname <- "frpfire"

plot.title <- "No. of wildfire occurrences in Dec 2006"

### DO NOT EDIT BELOW THIS LINE ###

###########
# IMPORTS #
###########

library(raster)
library(rasterVis)
library(RColorBrewer)
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

# Plot flattened maps
plotMap <- function(raster, title, lines) {
  my.at <- c(0, 1, 2, 5, 10, 15, 20, 25, 31)
  my.brks <- seq(0, 9, length.out = 9)
  my.color.key <- list(at = my.brks, labels = list(at = my.brks, labels = my.at))
  map.theme <- rasterTheme(region = c("#006400", brewer.pal(8, "Oranges")), panel.background = list(col = "skyblue"))
  levelplot(
    raster,
    par.settings = map.theme, 
    at = my.at,
    colorkey = my.color.key,
    main = title,
    margin = FALSE,
  ) + lines
}

########
# PLOT #
########

setwd(working.dir)
cpoly <- getData("GADM", download = FALSE, path = data.dir, country = cname, level = 1)
# Simplify the geometry of the polygon shapes
cpoly.simplified <- gSimplify(getSmallPolys(cpoly, minarea=1.875e-2), tol = 1e-2, topologyPreserve = TRUE)
c.var <- raster::brick(fname, varname = vname)
c.var.flat <- flattenRaster(c.var, cpoly, function(x, na.rm) {sum(x > 0, na.rm = na.rm)})
plotMap(crop(mask(c.var.flat, cpoly.simplified), cpoly.simplified), plot.title, layer(sp.lines(cpoly.simplified)))
