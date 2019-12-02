library(maps)
library(maptools)
library(raster)
library(rasterVis)
library(RColorBrewer)
library(sp)

frpfire.gfas.190717 <- raster("/media/zegheim/Justin/nc/cams_gfas_ga_1907.nc", band = 17, varname = "frpfire")
frpfire.gfas.190717 <- rotate(frpfire.gfas.190717)

# Convert small values to 0
frpfire.gfas.190717[frpfire.gfas.190717 < 1e-10] <- 0

# Checking the coordinates of maximum value
xyFromCell(frpfire.gfas.190717,  which.max(frpfire.gfas.190717))

# Create a world map SpatialLines object
countries <- map("world", plot=FALSE) 
countries <- map2SpatialLines(countries, proj4string = CRS("+proj=longlat"))

my.at <- c(0, 1e-2, 1e-1, 2.5e-1, 5e-1, 7.5e-1, 1, 2.5, 5, 10)
my.brks <- seq(0, 10, length.out = 10)
my.color.key <- list(at = my.brks, labels = list(at = my.brks, labels = my.at))
map.theme <- rasterTheme(region = c("#00688B", brewer.pal(9, "Oranges")))

levelplot(frpfire.gfas.190717, par.settings = map.theme, at = my.at, colorkey = my.color.key, margin = FALSE) + layer(sp.lines(countries))