library(raster)
library(latex2exp)
library(rasterVis)
library(RColorBrewer)
library(gridExtra)
library(stringr)

getRasterFromDataFrame <- function(coords, data, crs = proj.str) {
  xyz <- coords
  xyz$z <- data
  return(rasterFromXYZ(xyz, crs = proj.str))
}

plotMap <- function(raster, colours, plot.title, layers, labels = NULL, bg.colour = "skyblue", padding = 0) {
  map.theme <- rasterTheme(region = colours, panel.background = list(col = bg.colour))
  map.theme$layout.heights[
    c(
      'bottom.padding', 
      'top.padding', 
      'key.sub.padding',
      'axis.xlab.padding',
      'key.axis.padding',
      'main.key.padding'
    )
    ] <- padding
  
  if (!is.null(labels)) {
    brks <- seq(0, length(labels), length.out = length(labels))
    
    l <- levelplot(
      raster,
      par.settings = map.theme,
      at = labels,
      colorkey = list(at = brks, labels = list(at = brks, labels = labels)),
      margin = FALSE,
    )
  } else {
    l <- levelplot(
      raster,
      par.settings = map.theme,
      margin = FALSE
    )
  }
  
  l$aspect.fill <- TRUE
  l$xlab <- NULL
  l$ylab <- NULL
  l$main <- plot.title
  
  return(l + layers)
}

plotMapFromDataFrame <- function(coords, data, colours, plot.title, layers, labels = NULL, bg.colour = "white", padding = 0.5, labels.gap = 1) {
  xyz.raster <- getRasterFromDataFrame(coords, data, proj.str)
  plotMap(xyz.raster, colours, plot.title, layers, labels = labels, bg.colour = bg.colour, padding = padding)
}

setLabelGaps <- function(plot, labels.gap = 5) {
  at <- plot$legend$right$args$key$labels$at
  plot$legend$right$args$key$labels$at <- at[seq(1, length(at), by = labels.gap)]
  labels <- plot$legend$right$args$key$labels$labels
  plot$legend$right$args$key$labels$labels <- labels[seq(1, length(labels), by = labels.gap)]
  return(plot)
}

proj.str <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
colours <- brewer.pal(9, "Greys")
labels <- c(0, 1, 2, 5, 10, 15, 20, 25, 31)

cname1 <- "Indonesia"
cpoly1 <- getData("GADM", path = str_glue("~/Documents/projects/diss/data/rds"), country = cname1, level = 1)
layers1 <- layer(sp.lines(cpoly1))
df1 = read.csv("~/Documents/projects/diss/data/csv/df_ina_lores.csv")
coords1 <- df1[c("x", "y")]

plot1 <- plotMapFromDataFrame(coords1, df1$count, colours, "", layers1, labels = labels)

cname2 <- "Australia"
cpoly2 <- getData("GADM", path = str_glue("~/Documents/projects/diss/data/rds"), country = cname2, level = 1)
layers2 <- layer(sp.lines(cpoly2))
df2 = read.csv("~/Documents/projects/diss/data/csv/df_aus_lores.csv")
coords2 <- df2[c("x", "y")]

plot2 <- plotMapFromDataFrame(coords2, df2$count, colours, "", layers2, labels = labels)

grid.arrange(plot1, plot2, ncol = 2)

labels.temp <- seq(10, 45)
labels.rhm <- seq(0, 1, 0.1)

climate11 <- plotMapFromDataFrame(coords1, df1$skt, colours, TeX("Surface temperature (\u00B0C)"), layers1, labels = labels.temp)
climate11 <- setLabelGaps(climate11)
climate12 <- plotMapFromDataFrame(coords1, df1$rhm, colours, TeX("Relative humidity"), layers1, labels = labels.rhm)
climate21 <- plotMapFromDataFrame(coords2, df2$skt, colours, TeX("Surface temperature (\u00B0C)"), layers2, labels = labels.temp)
climate21 <- setLabelGaps(climate21)
climate22 <- plotMapFromDataFrame(coords2, df2$rhm, colours, TeX("Relative humidity"), layers2, labels = labels.rhm)

plots.climate <- list(climate11, climate12, climate21, climate22)
do.call("grid.arrange", c(plots.climate, list(left = "Latitude", bottom = "Longitude")))

###

labels.elev <- seq(0, 2500, 100)
labels.ptc <- seq(0, 1, 0.1)

topograph11 <- plotMapFromDataFrame(coords1, df1$elev, colours, TeX("Elevation above sea level (m)"), layers1, labels = labels.elev)
topograph11 <- setLabelGaps(topograph11)
topograph12 <- plotMapFromDataFrame(coords1, df1$ptc, colours, TeX("Percent tree cover"), layers1, labels = labels.ptc)
topograph21 <- plotMapFromDataFrame(coords2, df2$elev, colours, TeX("Elevation above sea level (m)"), layers2, labels = labels.elev)
topograph21 <- setLabelGaps(topograph21)
topograph22 <- plotMapFromDataFrame(coords2, df2$ptc, colours, TeX("Percent tree cover"), layers2, labels = labels.ptc)

plots.topograph <- list(topograph11, topograph12, topograph21, topograph22)
do.call("grid.arrange", c(plots.topograph, list(left = "Latitude", bottom = "Longitude")))

