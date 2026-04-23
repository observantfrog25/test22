library(terra)
library(here)
library(dplyr)

#' Build a prediction grid for a single state and year.
#' Returns a SpatRaster with layers: year, median_aadt, tmin, tmean, ppt,
#' FctImp, LndCov -- ready for predict(model, grid).
#'
#' @param state Two-letter abbreviation: "CA", "OR", or "WA"
#' @param year  Prediction year (integer)
#' @param chunk_pattern  Glob pattern for the traffic flow chunk RDS files.
#'   Must contain one "\%s" for the state abbreviation.
#' @param exclude_lndcov LndCov values to set to NA (e.g. c(12, 41) for
#'   classes that cause complete separation in the model)
build_pred_grid <- function(state, year,
                            chunk_pattern = "Roads/traffic_flow_extracted/%s",
                            exclude_lndcov = c(12, 41)) {

  state_upper <- toupper(state)
  grid_file <- here::here(paste0(state_upper, "5km_raster.tif" ))
  base_grid <- terra::rast(grid_file)

  # --- Traffic flow from chunks ---
  chunk_dir <- here::here(sprintf(chunk_pattern, state_upper))
  chunk_files <- list.files(
    chunk_dir,
    pattern = paste0(state_upper, "_2024_chunk_.*\\.rds"),
    full.names = TRUE
  )

  mean_aadt_rast     <- base_grid * NA
  median_aadt_rast   <- base_grid * NA
  harmonic_aadt_rast <- base_grid * NA
  dist_rast          <- base_grid * NA

  for (f in chunk_files) {
    chunk <- readRDS(f)
    mean_aadt_rast[chunk$cells]     <- chunk$metrics$mean_aadt
    median_aadt_rast[chunk$cells]   <- chunk$metrics$median_aadt
    harmonic_aadt_rast[chunk$cells] <- chunk$metrics$harmonic_aadt
    dist_rast[chunk$cells]          <- chunk$metrics$dist_d
  }

  exposure <- c(mean_aadt_rast, median_aadt_rast, harmonic_aadt_rast, dist_rast)
  names(exposure) <- c("mean_aadt", "median_aadt", "harmonic_aadt", "dist_d")

  # --- Climate ---
  tmin  <- terra::rast(here::here("TMIN/Pre stack",  paste0(state_upper, year, ".tif")))
  tmean <- terra::rast(here::here("TMEAN/Pre stack", paste0(state_upper, year, ".tif")))
  ppt   <- terra::rast(here::here("PPT",
                                  paste0(state_upper, " pre-stack"),
                                  paste0(state_upper, year, ".tif")))
  names(tmin)  <- "tmin"
  names(tmean) <- "tmean"
  names(ppt)   <- "ppt"

  # --- NLCD (use last available layer; repeat for future years) ---
  fctimp_stack <- terra::rast(here::here("NLCD/state_FctImp", paste0(state_upper, ".tif")))
  lndcov_stack <- terra::rast(here::here("NLCD/state_LndCov", paste0(state_upper, ".tif")))

  fctimp <- fctimp_stack[[terra::nlyr(fctimp_stack)]]
  lndcov <- lndcov_stack[[terra::nlyr(lndcov_stack)]]

  names(fctimp) <- "FctImp"
  names(lndcov) <- "LndCov"

  if (length(exclude_lndcov) > 0) {
    lndcov <- terra::ifel(lndcov %in% exclude_lndcov, NA, lndcov)
  }

  # --- Year layer ---
  year_rast <- base_grid * 0 + year
  names(year_rast) <- "year"

  # --- Stack ---
  grid <- c(year_rast, exposure[["median_aadt"]], tmin, tmean, ppt, fctimp, lndcov)
  terra::time(grid) <- rep(year, terra::nlyr(grid))

  grid
}

#' Build prediction grids for multiple years and stack them into a list.
#'
#' @param state Two-letter abbreviation
#' @param years Integer vector of years
#' @param ... Additional arguments passed to build_pred_grid
#' @return Named list of SpatRasters, one per year
build_pred_stack <- function(state, years, ...) {
  grids <- lapply(years, function(yr) {
    message(paste("Building", state, yr, "grid..."))
    build_pred_grid(state, yr, ...)
  })
  names(grids) <- years
  grids
}
