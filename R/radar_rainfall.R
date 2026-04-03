#' Extract areal rainfall from a radar or gridded rainfall stack
#'
#' Computes the catchment-average rainfall for each time step in a terra
#' `SpatRaster` stack by taking the area-weighted mean of all grid cells that
#' intersect the catchment polygon.  Uses `exactextractr::exact_extract()` to
#' handle partial cells at catchment boundaries correctly — this is important
#' for small catchments relative to the radar grid resolution (1 km for NIMROD,
#' 1 km for CEH-GEAR).
#'
#' @param catchment  An sf polygon from [load_catchment()].
#' @param radar_stack  A terra `SpatRaster` with one layer per time step.
#'   Time metadata (if set via `terra::time()`) is used to label the output.
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{`timestep`}{Integer index (1, 2, ..., `nlyr(radar_stack)`).}
#'     \item{`datetime`}{POSIXct timestamp (if available from `terra::time()`),
#'       otherwise `NA`.}
#'     \item{`areal_rainfall_mm`}{Catchment-average rainfall in mm.}
#'   }
#'
#' @seealso [areal_reduction_factor()]
#'
#' @examples
#' \dontrun{
#' catchment  <- load_catchment("wye.shp")
#' nimrod_stk <- terra::rast("nimrod_2024.nc")
#' areal      <- extract_radar_rainfall(catchment, nimrod_stk)
#' }
#'
#' @export
extract_radar_rainfall <- function(catchment, radar_stack) {
  if (!requireNamespace("exactextractr", quietly = TRUE)) {
    stop(
      "Package 'exactextractr' is required for extract_radar_rainfall().\n",
      "Install it with: install.packages('exactextractr')",
      call. = FALSE
    )
  }
  if (!inherits(radar_stack, "SpatRaster")) {
    stop("radar_stack must be a terra SpatRaster.", call. = FALSE)
  }

  # Reproject catchment to match raster CRS if necessary
  rast_crs <- terra::crs(radar_stack, describe = TRUE)
  catch_reproj <- sf::st_transform(catchment, crs = terra::crs(radar_stack))

  # exactextractr returns a list (one element per feature); catchment is 1 row
  result <- exactextractr::exact_extract(radar_stack, catch_reproj, "mean",
                                         progress = FALSE)

  # result is a data.frame with one row and one column per layer
  means <- as.numeric(result[1L, ])

  n_layers <- terra::nlyr(radar_stack)
  times    <- tryCatch(terra::time(radar_stack), error = function(e) NULL)

  data.frame(
    timestep          = seq_len(n_layers),
    datetime          = if (!is.null(times) && length(times) == n_layers) times
                        else rep(NA_real_, n_layers),
    areal_rainfall_mm = means,
    stringsAsFactors  = FALSE
  )
}


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
  # FEH Table 2.3 ARF values
  # Rows = area breaks (km²), Columns = duration breaks (hours)
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

  # Clamp area with warning
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

  # Vectorise over recycled inputs
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
  # Find bounding indices for area
  ai_lo <- findInterval(area, area_breaks, all.inside = TRUE)
  ai_hi <- min(ai_lo + 1L, length(area_breaks))

  # Find bounding indices for duration
  di_lo <- findInterval(dur, dur_breaks, all.inside = TRUE)
  di_hi <- min(di_lo + 1L, length(dur_breaks))

  # Fraction along each axis
  t_a <- if (area_breaks[ai_hi] == area_breaks[ai_lo]) 0 else
    (area - area_breaks[ai_lo]) / (area_breaks[ai_hi] - area_breaks[ai_lo])
  t_d <- if (dur_breaks[di_hi] == dur_breaks[di_lo]) 0 else
    (dur - dur_breaks[di_lo]) / (dur_breaks[di_hi] - dur_breaks[di_lo])

  # Bilinear interpolation
  v00 <- tbl[ai_lo, di_lo]
  v10 <- tbl[ai_hi, di_lo]
  v01 <- tbl[ai_lo, di_hi]
  v11 <- tbl[ai_hi, di_hi]

  (1 - t_a) * (1 - t_d) * v00 +
    t_a * (1 - t_d) * v10 +
    (1 - t_a) * t_d * v01 +
    t_a * t_d * v11
}
