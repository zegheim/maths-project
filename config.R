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
  cname <- "Australia"
  ccode <- "aus"
  
  is.lowres <- TRUE
  lowres.factor <- 5
  
  month <- 12
  year <- 15
  month.str <- formatC(month, digits = 1, flag = "0", format = "d")
  year.str <- formatC(year, digits = 1, flag = "0", format = "d")
  
  storage.dir <- str_glue("/media/zegheim/Justin_SSD/nc_{ccode}")
  fname.fire <- str_glue("{storage.dir}/gfas/cams_gfas_ga_{year.str}{month.str}.nc")
  fname.era5 <- str_glue("{storage.dir}/era5/era5.land.data.nc")
  fname.airt <- str_glue("{storage.dir}/era5/t2m.nc")
  fname.temp <- str_glue("{storage.dir}/tair/tair_2015_cropped.nc")
  fname.dtrg <- str_glue("{storage.dir}/ceda/dtrg/cru_ts4.03.2001.2018.dtr.dat.nc")
  fname.prec <- str_glue("{storage.dir}/ceda/prec/cru_ts4.03.2001.2018.pre.dat.nc")
  
  vname.fire <- "frpfire"
  vname.prec <- "pre"
  vname.dewp <- "d2m"
  vname.dtrg <- "dtr"
  vname.airt <- "t2m"
  
  csv.name <- str_glue("df_{ccode}_{ifelse(is.lowres, 'lores', 'hires')}.csv")
}

# MODEL FITTING -----------------------------------------------------------

if (getOption("run.model_fitting", default = FALSE)) {
  csv.name <- "df_ina_lores.csv"
  RData.name <- str_glue("result.ina.{strftime(Sys.time(), format = '%Y%m%d%H%M%S')}.RData")
  res <- 0.5
  seed <- 17071996L
}

# MODEL CHECKING ----------------------------------------------------------

if (getOption("run.model_checking", default = FALSE)) {
  fname <- "result.ina.20200311163133.RData"
  proj.str <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
}

# SIMULATION --------------------------------------------------------------

if (getOption("run.simulation", default = FALSE)) {
  fname <- "result.ina.20200311163133.RData"
  proj.str <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  RData.name <- str_glue("sim.ina.{strftime(Sys.time(), format = '%Y%m%d%H%M%S')}.RData")
  seed <- 26031997L
}
