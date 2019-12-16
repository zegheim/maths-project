library(animation)
library(maps)
library(maptools)
library(raster)
library(rasterVis)
library(RColorBrewer)
library(sp)

setwd("~/Documents/diss/")

ina <- getData("GADM", country = "Indonesia", level = 1)

frpfire.1807 <- raster::brick("/media/zegheim/Justin/nc/cams_gfas_ga_1807.nc")
frpfire.1807.cropped <- crop(frpfire.1807, ina)
frpfire.1807.cropped.masked <- mask(frpfire.1807.cropped, ina)

# Create a world map SpatialLines object
countries <- map("world", plot=FALSE) 
countries <- map2SpatialLines(countries, proj4string = CRS("+proj=longlat"))

my.at <- c(0, 1e-2, 1e-1, 2.5e-1, 5e-1, 7.5e-1, 1, 2.5, 5, 10)
my.brks <- seq(0, 10, length.out = 10)
my.color.key <- list(at = my.brks, labels = list(at = my.brks, labels = my.at))
map.theme <- rasterTheme(region = c("#000000", brewer.pal(9, "Oranges")))

levelplot(frpfire.1807.cropped.masked[[17]], par.settings = map.theme, at = my.at, colorkey = my.color.key, margin = FALSE) + layer(sp.lines(ina))

saveGIF({
  for (i in 1:nlayers(frpfire.gfas.1807)) {
    l <- levelplot(
      frpfire.gfas.1807[[i]], 
      par.settings = map.theme, 
      at = my.at, 
      colorkey = my.color.key, 
      margin = FALSE) + 
      layer(sp.lines(countries))
    plot(l)
  }
}, interval = 0.2, movie.name = "frpfire.1807.gif", ani.height = 500, ani.width = 1000)


