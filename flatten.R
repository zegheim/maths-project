library(raster)

setwd("~/Documents/diss")
ina <- getData("GADM", download = FALSE, path = '~/Documents/diss/data/', country = "Indonesia", level = 1)
ina.simplified <- gSimplify(getSmallPolys(ina), tol = 0.01, topologyPreserve = TRUE)
frpfire.ina <- raster::brick("~/Documents/diss/data/cams_gfas_ga.nc", varname = "frpfire")
frpfire.ina.masked <- mask(frpfire.ina, ina.simplified)
# Only include values between June and October 2015
frpfire.summer <- dropLayer(frpfire.ina.masked, c(1:4534, 4688:6178))

# Count no. of non-zero values at each grid
args <- as.list(frpfire.summer)
args$fun <- function(x, na.rm) {sum(x > 0, na.rm = na.rm)}
args$na.rm <- TRUE
frpfire.flattened <- do.call(mosaic, args)
frpfire.flattened.masked <- mask(frpfire.flattened, ina.simplified)

# Plot flattened map
my.at <- c(0, 1, 2, 5, 10, 25, 50, 100)
my.brks <- seq(0, 8, length.out = 8)
my.color.key <- list(at = my.brks, labels = list(at = my.brks, labels = my.at))
map.theme <- rasterTheme(region = c("#006400", brewer.pal(7, "Oranges")), panel.background = list(col = "skyblue"))
ina.simplified.lines <- layer(sp.lines(ina.simplified))

levelplot(
  frpfire.flattened.masked, 
  par.settings = map.theme, 
  at = my.at, 
  colorkey = my.color.key,
  main = "No. of wildfire occurences between Jun - Oct 2015",
  margin = FALSE
) + ina.simplified.lines