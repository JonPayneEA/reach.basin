#' Compute elevation-band weights for rain gauges
#'
#' Uses a hypsometric curve (area–elevation relationship) to assign each gauge
#' a weight proportional to the fraction of the catchment area within that
#' gauge's elevation band.
#'
#' Elevation bands are defined by the mid-points between successive gauge
#' elevations, with the catchment minimum and maximum as the outer bounds.
#' Each band's area fraction is read from the hypsometric curve by linear
#' interpolation.
#'
#' @param gauges  A data.frame with at minimum a `notation` column and an
#'   `elevation` column (metres above OD).
#' @param hypsometric_data  A data.frame with columns `elevation` (metres
#'   above OD) and `area_fraction` (cumulative fraction of catchment area at
#'   or above that elevation, in [0, 1]).  Typically derived from a DEM using
#'   `whitebox` or equivalent.
#'
#' @return A data.frame with columns `notation` and `elevation_weight`.  Weights
#'   sum to 1.
#'
#' @seealso [gauge_to_catchment_weights()], [plot_hypsometric_curve()]
#'
#' @examples
#' \dontrun{
#' hyps <- data.frame(
#'   elevation     = c(100, 200, 300, 400, 500),
#'   area_fraction = c(1.0, 0.75, 0.50, 0.25, 0.0)
#' )
#' gauges <- data.frame(notation = c("A", "B"), elevation = c(150, 350))
#' elevation_weights(gauges, hyps)
#' }
#'
#' @export
elevation_weights <- function(gauges, hypsometric_data) {
  required_g <- c("notation", "elevation")
  required_h <- c("elevation", "area_fraction")
  miss_g <- setdiff(required_g, names(gauges))
  miss_h <- setdiff(required_h, names(hypsometric_data))
  if (length(miss_g)) stop("gauges is missing columns: ", paste(miss_g, collapse = ", "), call. = FALSE)
  if (length(miss_h)) stop("hypsometric_data is missing columns: ", paste(miss_h, collapse = ", "), call. = FALSE)

  # Sort gauges by elevation
  g    <- gauges[order(gauges$elevation), ]
  elev <- g$elevation
  n    <- nrow(g)

  hyps <- hypsometric_data[order(hypsometric_data$elevation), ]

  # Helper: interpolate area_fraction at a given elevation from the hyps curve
  frac_at <- function(elev_val) {
    stats::approx(
      x    = hyps$elevation,
      y    = hyps$area_fraction,
      xout = elev_val,
      rule = 2L   # extrapolate with boundary values
    )$y
  }

  # Band boundaries: min of hyps / midpoints / max of hyps
  hyps_min <- min(hyps$elevation)
  hyps_max <- max(hyps$elevation)

  lower_bounds <- c(hyps_min, (elev[-n] + elev[-1L]) / 2)
  upper_bounds <- c((elev[-n] + elev[-1L]) / 2, hyps_max)

  # Area fraction in each band = difference in cumulative fraction
  # (area_fraction is fraction ABOVE elevation, so higher elev → lower frac)
  band_fracs <- frac_at(lower_bounds) - frac_at(upper_bounds)
  band_fracs <- pmax(band_fracs, 0)  # guard against floating-point negatives

  total <- sum(band_fracs)
  if (total == 0) {
    stop("All elevation bands have zero area — check that gauge elevations overlap the hypsometric curve.", call. = FALSE)
  }

  data.frame(
    notation         = g$notation,
    elevation_weight = band_fracs / total,
    stringsAsFactors = FALSE
  )
}


#' Compute catchment-area weights for gauges (unified interface)
#'
#' A single entry point that dispatches to the appropriate weighting method and
#' returns a consistent named numeric vector.  Callers can swap methods without
#' changing any downstream code.
#'
#' @param gauges  An sf POINT object or data.frame with coordinate columns.
#'   Must contain a `notation` column.  For `method = "elevation"`, must also
#'   contain an `elevation` column.
#' @param catchment  An sf polygon object from [load_catchment()].
#' @param method  Character.  One of:
#'   \describe{
#'     \item{`"thiessen"`}{Voronoi area weights via [thiessen_weights()].}
#'     \item{`"elevation"`}{Hypsometric band weights via [elevation_weights()].
#'       Requires the `hypsometric_data` argument.}
#'     \item{`"distance_decay"`}{Inverse-distance-squared weights.}
#'   }
#' @param ...  Additional arguments passed to the underlying weighting function
#'   (`hypsometric_data` for `method = "elevation"`; `decay_exp` for
#'   `method = "distance_decay"`).
#'
#' @return A named numeric vector of weights (names = `notation` values).
#'   Weights sum to 1.
#'
#' @seealso [thiessen_weights()], [elevation_weights()]
#'
#' @examples
#' \dontrun{
#' catchment <- load_catchment("wye.shp")
#' gauges    <- reach.io::find_stations(names = "wye")
#'
#' # Thiessen
#' w_t <- gauge_to_catchment_weights(gauges, catchment, method = "thiessen")
#'
#' # Elevation
#' hyps <- data.frame(elevation = c(100, 500), area_fraction = c(1, 0))
#' w_e <- gauge_to_catchment_weights(gauges, catchment, method = "elevation",
#'                                   hypsometric_data = hyps)
#'
#' # Distance decay
#' w_d <- gauge_to_catchment_weights(gauges, catchment, method = "distance_decay")
#' }
#'
#' @export
gauge_to_catchment_weights <- function(gauges, catchment, method = "thiessen", ...) {
  method <- match.arg(method, c("thiessen", "elevation", "distance_decay"))

  dots <- list(...)

  w <- switch(method,
    thiessen = {
      df <- thiessen_weights(gauges, catchment)
      stats::setNames(df$thiessen_weight, df$notation)
    },
    elevation = {
      if (!"hypsometric_data" %in% names(dots)) {
        stop("method = 'elevation' requires a 'hypsometric_data' argument.", call. = FALSE)
      }
      df <- elevation_weights(gauges, dots$hypsometric_data)
      stats::setNames(df$elevation_weight, df$notation)
    },
    distance_decay = {
      decay_exp <- dots$decay_exp %||% 2
      .distance_decay_weights(gauges, catchment, decay_exp = decay_exp)
    }
  )

  w / sum(w)  # normalise (guard against floating-point drift)
}


# Internal: inverse-distance-decay weights to catchment centroid.
#' @keywords internal
#' @noRd
.distance_decay_weights <- function(gauges, catchment, decay_exp = 2) {
  pts     <- as_sf_points(gauges, crs = sf::st_crs(catchment))
  centroid <- sf::st_centroid(sf::st_union(catchment))

  dists   <- as.numeric(sf::st_distance(pts, centroid))
  dists[dists == 0] <- 1e-6  # avoid division by zero

  raw_w <- 1 / dists^decay_exp
  notation <- if ("notation" %in% names(gauges)) gauges$notation else seq_len(nrow(pts))
  stats::setNames(raw_w / sum(raw_w), notation)
}
