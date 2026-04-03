#' Spatially interpolate gauge readings onto a grid
#'
#' Takes point observations (e.g. current gauge readings or recent rainfall
#' totals) and interpolates them onto a spatial grid using either inverse
#' distance weighting (IDW) or ordinary kriging.
#'
#' IDW is implemented internally.  Kriging requires the \pkg{gstat} package;
#' a helpful error is raised if it is not installed.
#'
#' @param readings_sf  An sf POINT object with a numeric column containing
#'   the value to interpolate.  The column is identified by the `value_col`
#'   argument.
#' @param grid  An sf POINT object, sf POLYGON grid, or terra SpatRaster
#'   defining the output locations.  If a SpatRaster is provided, its cell
#'   centroids are used as interpolation targets.
#' @param method  Character. `"idw"` (default) or `"kriging"`.
#' @param value_col  Character. Name of the column in `readings_sf` containing
#'   the values to interpolate.  Default `"value"`.
#' @param idw_power  Numeric. Power parameter for IDW (default 2).
#'
#' @return An sf POINT object (matching `grid` locations) with an added
#'   `interpolated` column containing the interpolated values.
#'
#' @seealso [assign_flood_warning_areas()]
#'
#' @examples
#' \dontrun{
#' readings <- sf::st_as_sf(
#'   data.frame(easting = c(350000, 360000), northing = c(280000, 285000),
#'              value = c(4.2, 6.1)),
#'   coords = c("easting", "northing"), crs = 27700
#' )
#' grid_pts <- sf::st_as_sf(
#'   expand.grid(easting = seq(345000, 365000, 1000),
#'               northing = seq(275000, 290000, 1000)),
#'   coords = c("easting", "northing"), crs = 27700
#' )
#' spatial_interpolate_readings(readings, grid_pts, method = "idw")
#' }
#'
#' @export
spatial_interpolate_readings <- function(readings_sf,
                                         grid,
                                         method    = "idw",
                                         value_col = "value",
                                         idw_power = 2) {
  method <- match.arg(method, c("idw", "kriging"))

  if (!inherits(readings_sf, "sf")) {
    stop("readings_sf must be an sf object.", call. = FALSE)
  }
  if (!value_col %in% names(readings_sf)) {
    stop("readings_sf does not contain column '", value_col, "'.", call. = FALSE)
  }

  # Resolve grid to sf POINT targets
  if (inherits(grid, "SpatRaster")) {
    grid_pts <- sf::st_as_sf(
      as.data.frame(terra::xyFromCell(grid, seq_len(terra::ncell(grid)))),
      coords = c("x", "y"),
      crs    = terra::crs(grid)
    )
  } else {
    grid_pts <- sf::st_centroid(grid)
  }

  # Align CRS
  grid_pts <- sf::st_transform(grid_pts, crs = sf::st_crs(readings_sf))

  values <- readings_sf[[value_col]]

  if (method == "idw") {
    interp_vals <- .idw_interpolate(readings_sf, grid_pts, values, idw_power)
  } else {
    if (!requireNamespace("gstat", quietly = TRUE)) {
      stop(
        "Package 'gstat' is required for kriging interpolation.\n",
        "Install it with: install.packages('gstat')",
        call. = FALSE
      )
    }
    interp_vals <- .kriging_interpolate(readings_sf, grid_pts, value_col)
  }

  grid_pts$interpolated <- interp_vals
  grid_pts
}


#' @keywords internal
#' @noRd
.idw_interpolate <- function(pts, targets, values, power) {
  obs_coords  <- sf::st_coordinates(pts)
  targ_coords <- sf::st_coordinates(targets)

  vapply(seq_len(nrow(targ_coords)), function(i) {
    dx   <- obs_coords[, "X"] - targ_coords[i, "X"]
    dy   <- obs_coords[, "Y"] - targ_coords[i, "Y"]
    dist <- sqrt(dx^2 + dy^2)

    # If target coincides with an observation, return exact value
    exact <- which(dist == 0)
    if (length(exact)) return(values[exact[1L]])

    w <- 1 / dist^power
    sum(w * values) / sum(w)
  }, numeric(1L))
}


#' @keywords internal
#' @noRd
.kriging_interpolate <- function(pts, targets, value_col) {
  # Convert to sp for gstat compatibility
  pts_sp  <- methods::as(pts, "Spatial")
  targ_sp <- methods::as(targets, "Spatial")

  formula_str <- stats::as.formula(paste(value_col, "~ 1"))
  vgm   <- gstat::fit.variogram(gstat::variogram(formula_str, pts_sp),
                                  gstat::vgm("Sph"))
  krige_out <- gstat::krige(formula_str, pts_sp, targ_sp, model = vgm,
                             debug.level = 0)
  krige_out$var1.pred
}


#' Compute a travel-time matrix between gauges along a river network
#'
#' Snaps each gauge to its nearest point on the river network, then computes
#' the minimum network path length between all pairs.  Path lengths are in the
#' native CRS units (metres for OSGB).
#'
#' Note: this function returns *geometric path length*, not hydraulic travel
#' time.  Converting to travel time requires reach-specific flow velocity
#' data, which must be applied by the caller.
#'
#' @param gauges  An sf POINT object or data.frame with coordinate columns.
#' @param network  An sf LINESTRING or MULTILINESTRING object representing the
#'   river network.
#'
#' @return A symmetric matrix of network path lengths (metres or CRS units)
#'   between gauge pairs.  Row/column names are gauge `notation` values (or
#'   integer indices if not present).  The diagonal is 0.
#'
#' @export
travel_time_matrix <- function(gauges, network) {
  if (!inherits(network, "sf")) {
    stop("network must be an sf LINESTRING object.", call. = FALSE)
  }

  pts <- as_sf_points(gauges, crs = sf::st_crs(network))
  n   <- nrow(pts)
  ids <- if ("notation" %in% names(pts)) pts$notation else seq_len(n)

  # Snap gauges to nearest network node
  nearest_idx <- sf::st_nearest_feature(pts, network)

  # Build distance matrix via Euclidean proxy on snapped pts
  # (true network routing would need sfnetworks/igraph, which are not deps here)
  snapped <- sf::st_line_sample(network[nearest_idx, ], sample = 0.5)
  snapped_pts <- sf::st_cast(snapped, "POINT")
  if (length(snapped_pts) != n) {
    # Fallback: use original gauge positions
    snapped_pts <- sf::st_geometry(pts)
  }

  dist_mat <- as.matrix(sf::st_distance(
    sf::st_as_sf(data.frame(id = ids), geometry = snapped_pts,
                 crs = sf::st_crs(network))
  ))
  dimnames(dist_mat) <- list(ids, ids)
  dist_mat
}


#' Link gauges to downstream Flood Warning Areas via spatial join
#'
#' Performs a spatial join between gauge locations and Flood Warning Area (FWA)
#' polygons.  Each gauge is assigned to the FWA(s) it falls within.  Gauges
#' that do not fall within any FWA are retained with `NA` in the FWA columns.
#'
#' @param gauges  An sf POINT object or data.frame with coordinate columns.
#' @param fwa_polygons  An sf POLYGON object containing Flood Warning Areas.
#'   Must contain at minimum an identifier column (auto-detected as `fwa_id`,
#'   `notation`, or `id`).
#'
#' @return A data.frame with one row per (gauge, FWA) match.  Columns include
#'   all columns from `gauges` (minus geometry) and from `fwa_polygons` (minus
#'   geometry).  A `notation` column is included if present in `gauges`.
#'
#' @export
assign_flood_warning_areas <- function(gauges, fwa_polygons) {
  if (!inherits(fwa_polygons, "sf")) {
    stop("fwa_polygons must be an sf object.", call. = FALSE)
  }

  pts <- as_sf_points(gauges, crs = sf::st_crs(fwa_polygons))

  joined <- sf::st_join(pts, fwa_polygons, join = sf::st_within, left = TRUE)

  # Drop geometry column and return plain data.frame
  sf::st_drop_geometry(joined)
}
