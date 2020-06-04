# IMPORTS -----------------------------------------------------------------

library(stringr)

# COMMON CONFIG -----------------------------------------------------------

sourceFunctions <- function(file) {
  options(run.main = FALSE)
  source(file)
  options(run.main = TRUE)
}

working.dir <- "~/Documents/diss"
data.dir <- str_glue("{working.dir}/data")

setwd(working.dir)

# PREPROCESSING -----------------------------------------------------------

if (getOption("run.preprocess", default = FALSE)) {
  cname <- "Indonesia"
  ccode <- "ina"
  
  is.lowres <- FALSE
  lowres.factor <- 5
  
  month <- 9
  year <- 15
  month.str <- formatC(month, digits = 1, flag = "0", format = "d")
  year.str <- formatC(year, digits = 1, flag = "0", format = "d")
  
  vname.fire <- "frpfire"
  vname.airt <- "t2m"
  vname.dewp <- "d2m"
  vname.temp <- "skt"
  vname.elev <- "elev"
  vname.vegc <- "ptc"
  
  storage.dir <- str_glue("/media/zegheim/Justin_SSD/nc_{ccode}")
  fname.fire <- str_glue("{storage.dir}/gfas/cams_gfas_ga_{year.str}{month.str}.nc")
  fname.airt <- str_glue("{storage.dir}/era5/{vname.airt}.{year.str}{month.str}.nc")
  fname.dewp <- str_glue("{storage.dir}/era5/{vname.dewp}.{year.str}{month.str}.nc")
  fname.temp <- str_glue("{storage.dir}/era5/{vname.temp}.{year.str}{month.str}.nc")
  fname.elev <- str_glue("{storage.dir}/globalmaps/elevation/gm_el_v2.nc")
  fname.vegc <- str_glue("{storage.dir}/globalmaps/ptc/gm_ve_v2.nc")
  
  csv.name <- str_glue("df_{ccode}_{ifelse(is.lowres, 'lores', 'hires')}.csv")
}

# MODEL FITTING -----------------------------------------------------------

if (getOption("run.model_fitting", default = FALSE)) {
  ccode <- "ina"
  is.lowres <- FALSE
  csv.name <- str_glue("df_{ccode}_{ifelse(is.lowres, 'lores', 'hires')}.csv")
  RData.name <- str_glue("result.{ccode}.{strftime(Sys.time(), format = '%Y%m%d%H%M%S')}.RData")
  res <- 0.1
  seed <- 17071996L
}

# MODEL CHECKING ----------------------------------------------------------

if (getOption("run.model_checking", default = FALSE)) {
  fname <- "result.aus.20200330231943.RData"
  proj.str <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
}

# SIMULATION --------------------------------------------------------------

if (getOption("run.simulation", default = FALSE)) {
  fname <- "result.ina.20200311163133.RData"
  proj.str <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  RData.name <- str_glue("sim.ina.{strftime(Sys.time(), format = '%Y%m%d%H%M%S')}.RData")
  seed <- 26031997L
}
