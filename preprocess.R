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
  flattened <- mask(do.call(mosaic, args), map)
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

###############
# MAIN SCRIPT #
###############

setwd(working.dir)

cpoly <- getData("GADM", download = FALSE, path = data.dir, country = cname, level = 1)
# Simplify the geometry of the polygon shapes
cpoly.simplified <- gSimplify(getSmallPolys(cpoly), tol = 0.01, topologyPreserve = TRUE)

c.var <- raster::brick(fname, varname = vname)
c.var.flat <- flattenRaster(c.var, cpoly.simplified, function(x, na.rm) {sum(x > 0, na.rm = na.rm)})
c.var.flat.lowres <- aggregate(c.var.flat, fact = 5, fun = sum)

# Plot flattened map
my.at <- c(0, 1, 2, 5, 10, 25, 50, 100, 200)
my.brks <- seq(0, 9, length.out = 9)
my.color.key <- list(at = my.brks, labels = list(at = my.brks, labels = my.at))
map.theme <- rasterTheme(region = c("#006400", brewer.pal(8, "Oranges")), panel.background = list(col = "skyblue"))
map.lines <- layer(sp.lines(cpoly.simplified))
levelplot(
  c.var.flat.lowres,
  par.settings = map.theme, 
  at = my.at,
  colorkey = my.color.key,
  main = plot.title,
  margin = FALSE
) + map.lines

c.df <- as.data.frame(c.var.flat.lowres, xy = TRUE, na.rm = TRUE)

# get elevation data
t <- get_elev_raster(c.var.flat.lowres, src = "aws", z = 6)
t.cropped <- resample(crop(t, c.var.flat), c.var.flat, method="bilinear")
t.cropped.lowres <- aggregate(t.cropped, fact = 5, fun = mean)
names(t.cropped.lowres) <- "elevation"

# get temperature data for July
temp <- raster::brick("/media/zegheim/Justin_SSD/cpc_global_temp/tavg.2015.nc")
temp.1507 <- dropLayer(temp, c(1:181, 213:365))
temp.1507.flattened <- flattenRaster(crop(temp.1507, c.var.flat.lowres), cpoly.simplified, mean)
temp.1507.resampled <- resample(temp.1507.flattened, c.var.flat.lowres, method = "bilinear")
names(temp.1507.resampled) <- "avg.temp"

# coerce to data.frame
test <- stack(c.var.flat.lowres, t.cropped.lowres, temp.1507.resampled)
test.df <- as.data.frame(test, xy = TRUE, na.rm = TRUE)
