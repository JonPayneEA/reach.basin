# ---- Internal helpers -------------------------------------------------------

#' @keywords internal
#' @noRd
#'
#' Resolve `files` to an ordered character vector of NetCDF paths, or return
#' a SpatRaster as-is.
#'
#' @param files  SpatRaster, character(1) directory, or character(≥1) file paths.
#' @param var_name  NetCDF sub-dataset / variable name.
#' @param expected_layers  Optional integer.  If given, warn when any file's
#'   layer count is not a multiple of this number.
#' @param context  Short string used in warning messages (e.g. "H19", "MOSES PE").
#' @param ...  Passed to terra::rast().
#' @return A terra SpatRaster concatenated across all files.
load_nc_stack <- function(files, var_name, expected_layers = NULL,
                          context = "raster", ...) {
  # --- SpatRaster: pass through unchanged ------------------------------------
  if (inherits(files, "SpatRaster")) {
    return(files)
  }

  if (!is.character(files) || length(files) == 0L) {
    stop(
      "'files' must be a SpatRaster, a directory path, or a character vector ",
      "of NetCDF file paths.",
      call. = FALSE
    )
  }

  # --- Directory path ---------------------------------------------------------
  if (length(files) == 1L && dir.exists(files)) {
    nc_files <- list.files(files, pattern = "\\.nc$", full.names = TRUE)
    if (length(nc_files) == 0L) {
      stop("No .nc files found in directory: ", files, call. = FALSE)
    }
    files <- nc_files
  }

  # Sort lexicographically so date-encoded filenames are in chronological order
  files <- sort(files)

  missing_files <- files[!file.exists(files)]
  if (length(missing_files)) {
    stop(
      "The following files do not exist:\n",
      paste(" ", missing_files, collapse = "\n"),
      call. = FALSE
    )
  }

  # --- Load each file and check layer count -----------------------------------
  rasters <- lapply(files, function(f) {
    r <- tryCatch(
      terra::rast(f, subds = var_name, ...),
      error = function(e) {
        stop(
          sprintf("Failed to open '%s' (subds = '%s'): %s", f, var_name, conditionMessage(e)),
          call. = FALSE
        )
      }
    )

    if (!is.null(expected_layers) && (terra::nlyr(r) %% expected_layers) != 0L) {
      warning(
        sprintf(
          "%s: file '%s' has %d layers; expected a multiple of %d.",
          context, basename(f), terra::nlyr(r), expected_layers
        ),
        call. = FALSE
      )
    }

    r
  })

  # Concatenate all layers
  do.call(terra::c, rasters)
}


#' @keywords internal
#' @noRd
#'
#' Shared extraction core: reproject catchment to raster CRS, run
#' exact_extract, build output data.frame.
#'
#' @param stack       SpatRaster already loaded.
#' @param catchment   sf polygon.
#' @param value_col   Name for the extracted values column in the result.
#' @return data.frame with columns timestep, datetime, <value_col>.
.extract_areal_mean <- function(stack, catchment, value_col) {
  if (!requireNamespace("exactextractr", quietly = TRUE)) {
    stop(
      "Package 'exactextractr' is required for this function.\n",
      "Install it with: install.packages('exactextractr')",
      call. = FALSE
    )
  }

  catch_reproj <- sf::st_transform(catchment, crs = terra::crs(stack))

  result <- exactextractr::exact_extract(stack, catch_reproj, "mean",
                                         progress = FALSE)

  means    <- as.numeric(result[1L, ])
  n_layers <- terra::nlyr(stack)
  times    <- tryCatch(terra::time(stack), error = function(e) NULL)

  out <- data.frame(
    timestep = seq_len(n_layers),
    datetime = if (!is.null(times) && length(times) == n_layers)
                 as.POSIXct(times)
               else
                 rep(as.POSIXct(NA), n_layers),
    stringsAsFactors = FALSE
  )
  out[[value_col]] <- means
  out
}


# ---- Public functions -------------------------------------------------------

#' Extract catchment-average rainfall from H19 (or other) radar NetCDF files
#'
#' Loads one or more gridded rainfall NetCDF files, concatenates them into a
#' single time series, and returns the area-weighted mean rainfall over a
#' catchment polygon for every timestep.
#'
#' Supports three input forms for `files`:
#' * A pre-loaded terra `SpatRaster` (backward-compatible with the original
#'   interface).
#' * A single directory path — all `.nc` files in the directory are loaded,
#'   sorted lexicographically (filename-date encoding makes this chronological).
#' * A character vector of individual file paths — sorted and loaded in order.
#'
#' Uses [exactextractr::exact_extract()] for area-weighted extraction, which
#' handles partial cells at catchment boundaries correctly.  This matters for
#' small catchments relative to the radar grid resolution (1 km for H19/NIMROD,
#' 1 km for CEH-GEAR).
#'
#' @param catchment  An sf polygon from [load_catchment()].
#' @param files  One of:
#'   * A terra `SpatRaster` (backward-compatible).
#'   * A length-1 character string giving a directory; all `.nc` files inside
#'     are loaded in lexicographic order.
#'   * A character vector of individual `.nc` file paths, loaded in
#'     lexicographic order.
#' @param product  Character. Radar product identifier.  Currently only
#'   `"H19"` is supported.  Used to validate expected file structure and will
#'   be used to select product-specific readers in future versions.
#' @param var_name  Character. NetCDF sub-dataset / variable name to extract
#'   from each file.  Default `"rainfall_amount"`.  Adjust if your files use
#'   a different variable name (e.g. `"rain"` or `"precip_rate"`).
#' @param ...  Additional arguments passed to [terra::rast()].
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{`timestep`}{Integer index (1 … total layers across all files).}
#'     \item{`datetime`}{POSIXct timestamp read from `terra::time()`, or `NA`
#'       if time metadata is absent from the files.}
#'     \item{`areal_rainfall_mm`}{Catchment-average rainfall in mm.}
#'   }
#'
#' @section H19 product notes:
#' H19 is a 15-minute NIMROD composite; each daily file contains 96 layers
#' (24 hours × 4 timesteps/hour).  A warning is issued for any file whose
#' layer count is not a multiple of 96.
#'
#' @seealso [extract_moses_pe()], [areal_reduction_factor()], [load_catchment()]
#'
#' @examples
#' \dontrun{
#' catchment <- load_catchment("wye.shp")
#'
#' # Single pre-loaded raster (original interface, still works)
#' stk   <- terra::rast("h19_20240101.nc", subds = "rainfall_amount")
#' areal <- extract_radar_rainfall(catchment, stk)
#'
#' # Directory of daily files
#' areal <- extract_radar_rainfall(catchment, "data/h19/jan2024/")
#'
#' # Explicit file list
#' files <- list.files("data/h19/", pattern = "\\.nc$", full.names = TRUE)
#' areal <- extract_radar_rainfall(catchment, files)
#' }
#'
#' @export
extract_radar_rainfall <- function(catchment,
                                   files,
                                   product  = "H19",
                                   var_name = "rainfall_amount",
                                   ...) {
  product <- match.arg(product, choices = "H19")

  expected_layers <- switch(product, H19 = 96L)

  stack <- load_nc_stack(
    files,
    var_name        = var_name,
    expected_layers = expected_layers,
    context         = product,
    ...
  )

  # CRS check: H19 should be OSGB
  rast_epsg <- terra::crs(stack, describe = TRUE)$code
  if (!isTRUE(as.integer(rast_epsg) == 27700L)) {
    warning(
      sprintf(
        "extract_radar_rainfall: raster CRS is EPSG:%s; expected EPSG:27700 (OSGB). ",
        rast_epsg %||% "unknown"
      ),
      "Catchment will be reprojected to match the raster.",
      call. = FALSE
    )
  }

  .extract_areal_mean(stack, catchment, value_col = "areal_rainfall_mm")
}


#' Extract catchment-average potential evapotranspiration from MOSES PE grids
#'
#' Loads one or more MOSES (Met Office Surface Exchange Scheme) PE NetCDF
#' files, concatenates them into a single hourly time series, and returns the
#' area-weighted mean PE over a catchment polygon for every timestep.
#'
#' MOSES PE is delivered as daily NetCDF files each containing 24 hourly
#' layers.  The function accepts the same flexible `files` input as
#' [extract_radar_rainfall()]: a pre-loaded `SpatRaster`, a directory path,
#' or a vector of file paths.
#'
#' @param catchment  An sf polygon from [load_catchment()].
#' @param files  One of:
#'   * A terra `SpatRaster` (backward-compatible).
#'   * A length-1 character string giving a directory; all `.nc` files inside
#'     are loaded in lexicographic order.
#'   * A character vector of individual `.nc` file paths.
#' @param var_name  Character. NetCDF variable name for PE.  Default `"peti"`
#'   (PE with interception correction, the variable used in PDM/FMP workflows).
#'   Other common values: `"pe"`, `"pet"`.
#' @param ...  Additional arguments passed to [terra::rast()].
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{`timestep`}{Integer index (1 … total layers across all files).}
#'     \item{`datetime`}{POSIXct timestamp from `terra::time()`, or `NA`.}
#'     \item{`areal_pe_mm`}{Catchment-average potential evapotranspiration
#'       in mm.}
#'   }
#'
#' @section MOSES PE product notes:
#' Each daily MOSES PE file is expected to contain 24 hourly layers.  A
#' warning is issued for any file whose layer count differs from 24.
#' The default variable name `"peti"` extracts PE with interception
#' correction, which is the standard input to PDM and FMP model runs.
#'
#' @seealso [extract_radar_rainfall()], [load_catchment()]
#'
#' @examples
#' \dontrun{
#' catchment <- load_catchment("wye.shp")
#'
#' # Directory of daily MOSES PE files
#' pe <- extract_moses_pe(catchment, "data/moses_pe/jan2024/")
#'
#' # Explicit file list with non-default variable name
#' pe <- extract_moses_pe(catchment, "data/moses_pe/jan2024/", var_name = "pe")
#' }
#'
#' @export
extract_moses_pe <- function(catchment,
                              files,
                              var_name = "peti",
                              ...) {
  stack <- load_nc_stack(
    files,
    var_name        = var_name,
    expected_layers = 24L,
    context         = "MOSES PE",
    ...
  )

  rast_epsg <- terra::crs(stack, describe = TRUE)$code
  if (!isTRUE(as.integer(rast_epsg) == 27700L)) {
    warning(
      sprintf(
        "extract_moses_pe: raster CRS is EPSG:%s; expected EPSG:27700 (OSGB). ",
        rast_epsg %||% "unknown"
      ),
      "Catchment will be reprojected to match the raster.",
      call. = FALSE
    )
  }

  .extract_areal_mean(stack, catchment, value_col = "areal_pe_mm")
}


# ---- ARF (unchanged) --------------------------------------------------------

#' Areal Reduction Factor (ARF) for design storm scaling
#'
#' Returns the FEH/NERC Areal Reduction Factor for a given catchment area and
#' storm duration.  ARF < 1 accounts for the spatial variability of rainfall —
#' a point estimate of rainfall is reduced when applied to a catchment area.
#'
#' Values are bilinearly interpolated from Table 2.3 of the Flood Estimation
#' Handbook (Reed 1999).  The table covers:
#' * Areas: 1, 5, 10, 30, 100, 300, 1 000, 3 000, 10 000 km²
#' * Durations: 1, 2, 6, 12, 24, 48, 96 hours
#'
#' For areas < 1 km² the ARF is extrapolated but will approach 1.0.  For areas
#' > 10 000 km² a warning is issued and the 10 000 km² row is used.
#'
#' @param area_km2    Numeric. Catchment area in km².
#' @param duration_hr Numeric. Storm duration in hours.
#'
#' @return Numeric vector of ARF values (same length as the longer of
#'   `area_km2` and `duration_hr`; both are recycled if necessary).
#'
#' @references
#' Reed, D.W. (1999) *Flood Estimation Handbook*, Vol. 1.  Centre for Ecology
#' and Hydrology.  Table 2.3.
#'
#' @examples
#' areal_reduction_factor(area_km2 = 250, duration_hr = 24)
#' areal_reduction_factor(area_km2 = c(10, 100, 500), duration_hr = 6)
#'
#' @export
areal_reduction_factor <- function(area_km2, duration_hr) {
  area_breaks     <- c(1, 5, 10, 30, 100, 300, 1000, 3000, 10000)
  duration_breaks <- c(1, 2, 6, 12, 24, 48, 96)

  # fmt: off
  arf_table <- matrix(c(
  # 1h     2h     6h     12h    24h    48h    96h
    0.978, 0.983, 0.990, 0.993, 0.996, 0.997, 0.998,  # 1 km²
    0.944, 0.953, 0.969, 0.977, 0.983, 0.988, 0.992,  # 5 km²
    0.923, 0.935, 0.955, 0.966, 0.975, 0.981, 0.987,  # 10 km²
    0.882, 0.899, 0.927, 0.943, 0.957, 0.967, 0.976,  # 30 km²
    0.837, 0.860, 0.897, 0.919, 0.938, 0.952, 0.964,  # 100 km²
    0.790, 0.818, 0.863, 0.890, 0.915, 0.933, 0.949,  # 300 km²
    0.734, 0.767, 0.820, 0.854, 0.884, 0.907, 0.928,  # 1000 km²
    0.672, 0.710, 0.772, 0.812, 0.848, 0.877, 0.903,  # 3000 km²
    0.604, 0.646, 0.717, 0.764, 0.806, 0.841, 0.874   # 10000 km²
  ), nrow = length(area_breaks), ncol = length(duration_breaks), byrow = TRUE)
  # fmt: on

  rownames(arf_table) <- area_breaks
  colnames(arf_table) <- duration_breaks

  if (any(area_km2 > 10000)) {
    warning(
      "area_km2 > 10 000 km\u00b2; ARF table only covers up to 10 000 km\u00b2. ",
      "Using 10 000 km\u00b2 row for out-of-range values.",
      call. = FALSE
    )
    area_km2 <- pmin(area_km2, 10000)
  }
  area_km2    <- pmax(area_km2, area_breaks[1L])
  duration_hr <- pmax(duration_hr, duration_breaks[1L])
  duration_hr <- pmin(duration_hr, duration_breaks[length(duration_breaks)])

  n <- max(length(area_km2), length(duration_hr))
  area_km2    <- rep_len(area_km2, n)
  duration_hr <- rep_len(duration_hr, n)

  vapply(seq_len(n), function(i) {
    .bilinear_arf(area_km2[i], duration_hr[i],
                  area_breaks, duration_breaks, arf_table)
  }, numeric(1L))
}


#' @keywords internal
#' @noRd
.bilinear_arf <- function(area, dur, area_breaks, dur_breaks, tbl) {
  ai_lo <- findInterval(area, area_breaks, all.inside = TRUE)
  ai_hi <- min(ai_lo + 1L, length(area_breaks))

  di_lo <- findInterval(dur, dur_breaks, all.inside = TRUE)
  di_hi <- min(di_lo + 1L, length(dur_breaks))

  t_a <- if (area_breaks[ai_hi] == area_breaks[ai_lo]) 0 else
    (area - area_breaks[ai_lo]) / (area_breaks[ai_hi] - area_breaks[ai_lo])
  t_d <- if (dur_breaks[di_hi] == dur_breaks[di_lo]) 0 else
    (dur - dur_breaks[di_lo]) / (dur_breaks[di_hi] - dur_breaks[di_lo])

  v00 <- tbl[ai_lo, di_lo]
  v10 <- tbl[ai_hi, di_lo]
  v01 <- tbl[ai_lo, di_hi]
  v11 <- tbl[ai_hi, di_hi]

  (1 - t_a) * (1 - t_d) * v00 +
    t_a * (1 - t_d) * v10 +
    (1 - t_a) * t_d * v01 +
    t_a * t_d * v11
}
