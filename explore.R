library(animation)
library(latticeExtra)
library(maps)
library(maptools)
library(raster)
library(rasterVis)
library(RColorBrewer)
library(rgeos)
library(sp)

setwd("~/Documents/diss/data")

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
ina <- getData("GADM", country = "Indonesia", level = 1)
ina.simplified <- gSimplify(getSmallPolys(ina), tol = 0.01, topologyPreserve = TRUE)

frpfire.ina <- raster::brick("cams_gfas_ga.nc", varname = "frpfire")
frpfire.ina.masked <- mask(frpfire.ina, ina.simplified)

my.at <- c(0, 1e-2, 1e-1, 2.5e-1, 5e-1, 7.5e-1, 1, 2.5, 5, 10)
my.brks <- seq(0, 10, length.out = 10)
my.color.key <- list(at = my.brks, labels = list(at = my.brks, labels = my.at))
map.theme <- rasterTheme(region = c("#006400", brewer.pal(9, "Oranges")), panel.background = list(col = "skyblue"))
ina.simplified.lines <- layer(sp.lines(ina.simplified))

dates <- strftime(strptime(names(frpfire.ina.masked), format = "X%Y.%m.%d"), format = "%d-%m-%Y")

i <- 4630
levelplot(
  frpfire.ina.masked[[i]], 
  par.settings = map.theme, 
  at = my.at, 
  colorkey = my.color.key,
  xlab.top = dates[i],
  margin = FALSE) + 
  ina.simplified.lines

saveGIF({
  for (i in 4384:4748) {
    l <- levelplot(
      frpfire.ina.masked[[i]], 
      par.settings = map.theme, 
      at = my.at, 
      colorkey = my.color.key,
      xlab.top = dates[i],
      margin = FALSE) + 
      ina.simplified.lines
    plot(l)
  }
}, 
interval = 0.1, 
movie.name = "frpfire.ina.2015.gif", 
ani.height = 540, 
ani.width = 960
)
  

