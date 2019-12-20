library(raster)
library(rasterVis)
library(RColorBrewer)
library(rgeos)

setwd("~/Documents/diss")

# code from https://gis.stackexchange.com/a/62405/155373 
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

ina <- getData("GADM", download = FALSE, path = '~/Documents/diss/data/', country = "Indonesia", level = 1)
ina.simplified <- gSimplify(getSmallPolys(ina), tol = 0.01, topologyPreserve = TRUE)
frpfire.ina <- raster::brick("~/Documents/diss/data/cams_gfas_ga.nc", varname = "frpfire")
frpfire.ina.masked <- mask(frpfire.ina, ina.simplified)

# Only include values between June and October 2015
frpfire.summer <- dropLayer(frpfire.ina.masked, c(1:4534, 4688:6178))
frpfire.jun <- dropLayer(frpfire.summer, 31:153)
frpfire.jul <- dropLayer(frpfire.summer, 1:30, 62:153)
frpfire.aug <- dropLayer(frpfire.summer, 1:61, 92:153)
frpfire.sep <- dropLayer(frpfire.summer, 1:91, 123:153)
frpfire.oct <- dropLayer(frpfire.summer, 1:122)

flattenRaster <- function(raster, map) {
  args <- as.list(raster)
  # Count no. of non-zero values at each grid
  args$fun <- function(x, na.rm) {sum(x > 0, na.rm = na.rm)}
  return(mask(do.call(mosaic, args), map))
}

frpfire.summer.flattened <- flattenRaster(frpfire.summer, ina.simplified)
frpfire.jun.flattened <- flattenRaster(frpfire.jun, ina.simplified)
frpfire.jul.flattened <- flattenRaster(frpfire.jul, ina.simplified)
frpfire.aug.flattened <- flattenRaster(frpfire.aug, ina.simplified)
frpfire.sep.flattened <- flattenRaster(frpfire.sep, ina.simplified)
frpfire.oct.flattened <- flattenRaster(frpfire.oct, ina.simplified)

# Plot flattened maps
plotMap <- function(raster, title) {
  my.at <- c(0, 1, 2, 5, 10, 25, 50, 100)
  my.brks <- seq(0, 8, length.out = 8)
  my.color.key <- list(at = my.brks, labels = list(at = my.brks, labels = my.at))
  map.theme <- rasterTheme(region = c("#006400", brewer.pal(7, "Oranges")), panel.background = list(col = "skyblue"))
  ina.simplified.lines <- layer(sp.lines(ina.simplified))
  levelplot(
    raster,
    par.settings = map.theme, 
    at = my.at,
    colorkey = my.color.key,
    main = title,
    margin = FALSE
  ) + ina.simplified.lines
}

plotMap(frpfire.summer.flattened, "No. of wildfire occurrences between Jun - Oct 2015")
plotMap(frpfire.jun.flattened, "No. of wildfire occurrences in Jun 2015")
plotMap(frpfire.jul.flattened, "No. of wildfire occurrences in Jul 2015")
plotMap(frpfire.aug.flattened, "No. of wildfire occurrences in Aug 2015")
plotMap(frpfire.sep.flattened, "No. of wildfire occurrences in Sep 2015")
plotMap(frpfire.oct.flattened, "No. of wildfire occurrences in Oct 2015")