#' Compute Thiessen (Voronoi) areal weights for rain gauges within a catchment
#'
#' Constructs Voronoi polygons from gauge locations, clips them to the
#' catchment boundary, and returns the proportion of catchment area attributed
#' to each gauge.  This replaces the three-step workflow previously done with
#' `mappER::teeSun()`, `mappER::intersectPoly()`, and `mappER::gaugeProp()`.
#'
#' Unlike the old `mappER` functions, this implementation:
#' * Handles MULTIPOLYGON catchments via [sf::st_cast()] before intersecting.
#' * **Errors** if any gauge falls outside the catchment bounding box, rather
#'   than silently assigning a zero weight, which caused hard-to-trace errors
#'   in operational workflows.
#'
#' @param coords  An sf POINT object, or a data.frame / data.table with
#'   coordinate columns (`easting`/`northing` in EPSG:27700, or `long`/`lat`
#'   in EPSG:4326).  Typically the output of `reach.io::find_stations()`.
#'   Must contain a `notation` column (station SUID) used as the gauge
#'   identifier in the output; if absent, a sequential integer id is added.
#' @param catchment  An sf polygon object, e.g. from [load_catchment()].
#'
#' @return A data.frame with one row per gauge and columns:
#'   \describe{
#'     \item{`notation`}{Station SUID (or integer index if not supplied).}
#'     \item{`label`}{Station name, if present in `coords`.}
#'     \item{`easting`, `northing`}{Gauge coordinates in EPSG:27700.}
#'     \item{`thiessen_weight`}{Area proportion [0, 1]; sums to 1.}
#'   }
#'
#' @seealso [saar_adjusted_weights()], [gauge_to_catchment_weights()]
#'
#' @examples
#' \dontrun{
#' catchment <- load_catchment("wye_catchment.shp")
#' gauges    <- reach.io::find_stations(names = "wye")
#' weights   <- thiessen_weights(gauges, catchment)
#' }
#'
#' @export
thiessen_weights <- function(coords, catchment) {
  # --- input coercion & CRS alignment ----------------------------------------
  pts <- as_sf_points(coords, crs = sf::st_crs(catchment))

  # --- bounding-box guard (error, not silent zero) ----------------------------
  bbox <- sf::st_bbox(catchment)
  coords_xy <- sf::st_coordinates(pts)

  outside <- (
    coords_xy[, "X"] < bbox["xmin"] |
    coords_xy[, "X"] > bbox["xmax"] |
    coords_xy[, "Y"] < bbox["ymin"] |
    coords_xy[, "Y"] > bbox["ymax"]
  )
  if (any(outside)) {
    bad <- which(outside)
    ids <- if ("notation" %in% names(pts)) pts$notation[bad] else bad
    stop(
      sprintf(
        "%d gauge(s) fall outside the catchment bounding box: %s\n",
        length(bad), paste(ids, collapse = ", ")
      ),
      "Remove these gauges or check coordinates before calling thiessen_weights().",
      call. = FALSE
    )
  }

  # --- Voronoi construction ---------------------------------------------------
  # st_voronoi needs a GEOMETRYCOLLECTION from st_combine; then explode back
  catchment_poly <- sf::st_cast(sf::st_union(catchment), "POLYGON")

  # Extend envelope slightly so edge gauges get full Voronoi tiles
  env <- sf::st_buffer(catchment_poly, dist = max(
    diff(bbox[c("xmin", "xmax")]),
    diff(bbox[c("ymin", "ymax")])
  ) * 0.1)

  voronoi_gc  <- sf::st_voronoi(sf::st_combine(pts), envelope = env)
  voronoi_sfc <- sf::st_cast(voronoi_gc, "POLYGON")

  # Match each Voronoi tile back to the nearest gauge point
  voronoi_sf <- sf::st_as_sf(data.frame(tile_id = seq_along(voronoi_sfc)),
                              geometry = voronoi_sfc,
                              crs = sf::st_crs(catchment))

  # Clip to catchment
  catchment_union <- sf::st_union(catchment)
  clipped <- sf::st_intersection(voronoi_sf, catchment_union)

  # st_voronoi returns tiles in the order of the input points
  # (guaranteed by GEOS ≥ 3.6), so tile_id == row index in pts.

  # --- area weights -----------------------------------------------------------
  total_area <- sum(sf::st_area(clipped))
  clipped$thiessen_weight <- as.numeric(sf::st_area(clipped) / total_area)

  # --- build output data.frame ------------------------------------------------
  coords_bng <- sf::st_transform(pts, crs = 27700)
  xy          <- sf::st_coordinates(coords_bng)

  out <- data.frame(
    notation        = if ("notation" %in% names(pts)) pts$notation
                      else seq_len(nrow(pts)),
    easting         = xy[, "X"],
    northing        = xy[, "Y"],
    thiessen_weight = clipped$thiessen_weight[order(clipped$tile_id)],
    stringsAsFactors = FALSE
  )

  if ("label" %in% names(pts)) {
    out <- cbind(data.frame(label = pts$label, stringsAsFactors = FALSE), out)
  }

  out
}
